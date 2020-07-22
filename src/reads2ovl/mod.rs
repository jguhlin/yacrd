/*
Copyright (c) 2019 Pierre Marijon <pmarijon@mpi-inf.mpg.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

/* crate use */
use anyhow::{anyhow, bail, Context, Result};

/* local mod */
pub mod fullmemory;
pub mod ondisk;

/* stuff declare in submod need to be accessible from mod level */
pub use self::fullmemory::*;
pub use self::ondisk::*;

/* std use */
pub use self::fullmemory::*;
use std::io::BufReader;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::fs::File;
use std::thread::Builder;
use std::thread;
use std::time::{Instant};
use std::sync::Arc;

/* external crates */
use twox_hash::XxHash64;
use flate2::bufread::GzDecoder;
use t1ha::{t1ha0};
use indicatif::{HumanDuration, MultiProgress, ProgressBar, ProgressStyle};
use sled;
use thincollections::thin_vec::ThinVec;
use crossbeam::utils::Backoff;
use crossbeam::queue::{ArrayQueue, PushError};
use byteorder::{BigEndian};
use serde::{Serialize, Deserialize};
use bincode;
use sled::{Batch, open};

use zerocopy::{byteorder::U64, 
    AsBytes, 
    FromBytes, 
    LayoutVerified, 
    Unaligned, 
    U16, 
    U32,};

pub const MAX_READS:  usize = 256_000_000;

#[derive(PartialEq)]
pub enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Vec<csv::StringRecord>> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> Vec<csv::StringRecord> {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<HashMap<String, ThinVec<(u32, u32)>, BuildHasherDefault<XxHash64>>> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> HashMap<String, ThinVec<(u32, u32)>, BuildHasherDefault<XxHash64>> {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}


/* local use */
use crate::error;
use crate::util;

pub trait Reads2Ovl {
    fn init(&mut self, filename: &str) -> Result<()> {
        self.sub_init(filename)
    }

    fn sub_init(&mut self, filename: &str) -> Result<()> {
        let (input, _) = util::read_file(filename)?;

        
        let reads2len_worker;
        {
            let filename = filename.to_string();
            reads2len_worker = thread::spawn(move || {
                parse_paf("output".to_string(), 64, 8 * 1024, filename)
            });
        }

        /*

        match util::get_file_type(filename) {
            Some(util::FileType::Paf) => self
                .init_paf(input)
                .with_context(|| anyhow!("Filename: {}", filename.to_string()))?,
            Some(util::FileType::M4) => self
                .init_m4(input)
                .with_context(|| anyhow!("Filename: {}", filename.to_string()))?,
            Some(util::FileType::Fasta) => bail!(error::Error::CantRunOperationOnFile {
                operation: "overlap parsing".to_string(),
                filetype: util::FileType::Fasta,
                filename: filename.to_string()
            }),
            Some(util::FileType::Fastq) => bail!(error::Error::CantRunOperationOnFile {
                operation: "overlap parsing".to_string(),
                filetype: util::FileType::Fastq,
                filename: filename.to_string()
            }),
            Some(util::FileType::Yacrd) => bail!(error::Error::CantRunOperationOnFile {
                operation: "overlap parsing".to_string(),
                filetype: util::FileType::Yacrd,
                filename: filename.to_string()
            }),
            None | Some(util::FileType::YacrdOverlap) => {
                bail!(error::Error::UnableToDetectFileFormat {
                    filename: filename.to_string()
                })
            }
        }

        println!("Joining reads2len thread...");

        */

        let reads2len = reads2len_worker.join().expect("Unable to join reads2len thread");

        Ok(())
    }

    fn init_m4(&mut self, input: Box<dyn std::io::Read>) -> Result<()> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b' ')
            .has_headers(false)
            .from_reader(input);

        for record in reader.records() {
            let result = record.with_context(|| error::Error::ReadingErrorNoFilename {
                format: util::FileType::M4,
            })?;

            if result.len() < 12 {
                bail!(error::Error::ReadingErrorNoFilename {
                    format: util::FileType::M4,
                });
            }

            let id_a = result[0].to_string();
            let id_b = result[1].to_string();

            let ovl_a = (util::str2u32(&result[5])?, util::str2u32(&result[6])?);
            let ovl_b = (util::str2u32(&result[9])?, util::str2u32(&result[10])?);

            self.add_overlap(id_a, ovl_a)?;
            self.add_overlap(id_b, ovl_b)?;
        }

        Ok(())
    }

    fn overlap(&self, id: &str) -> Result<Vec<(u32, u32)>>;
    fn length(&self, id: &str) -> usize;
    fn init_paf(&mut self, input: Box<dyn std::io::Read>) -> Result<()>;

    fn add_overlap(&mut self, id: String, ovl: (u32, u32)) -> Result<()>;

    fn get_reads(&self) -> std::collections::HashSet<String>;
}

// Aiming for single threaded here...
// HashMaps are a bit slow, here we only need to insert an item only once, and never update it
// So we do it "manually"
pub struct FastReadsIdx {
    pub idx:  Vec<Option<core::num::NonZeroU64>>,
    pub ids:  Vec<Option<String>>,
    pub lens: Vec<Option<core::num::NonZeroU32>>,
    pub entries: usize,
}

impl FastReadsIdx {
    pub fn convert_to_hashmap(&mut self) -> HashMap<String, u32, BuildHasherDefault<XxHash64>> {
        let mut reads2len: HashMap<String, u32, BuildHasherDefault<XxHash64>> = Default::default();

        reads2len.reserve(self.entries);

        reads2len = self.idx.iter().filter(|x| x.is_some()).map(|x| {
            let id = x.unwrap().get() as usize;
            let readid = self.ids[id].clone().unwrap();
            let len = self.lens[id].unwrap();
            (readid, len.get())
        }).collect();

        reads2len

    }

    pub fn new() -> FastReadsIdx {
        // Overkill, but multi-threading gives it a slight speed boost

        let idx_builder = match Builder::new()
                        .name("Idx Builder".into())
                        .spawn(|| 
                        {
                            let mut idx = Vec::with_capacity(MAX_READS);
                            idx.resize(MAX_READS, None);
                            idx
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

        let lens_builder = match Builder::new()
                        .name("Lens Builder".into())
                        .spawn(|| 
                        {
                            let mut lens = Vec::with_capacity(MAX_READS);
                            lens.resize(MAX_READS, None);
                            lens
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };

        let ids_builder = match Builder::new()
                        .name("Ids Builder".into())
                        .spawn(|| 
                        {
                            let mut ids = Vec::with_capacity(MAX_READS);
                            ids.resize(MAX_READS, None);
                            ids
                        })
                    {
                        Ok(x)  => x,
                        Err(y) => panic!("{}", y)
                    };


        let idx = idx_builder.join().unwrap();
        let lens = lens_builder.join().unwrap();
        let ids = ids_builder.join().unwrap(); 
        let entries = 0;

        FastReadsIdx { idx, lens, ids, entries }
    }

    // Code stolen from another project of mine, made for multithreaded use, so may be a few artifacts...
    #[inline]
    fn get_id(&self, readid: &str) -> usize {

        let hash = t1ha0(readid.as_bytes(), 42_988_123) as usize % MAX_READS;

        let mut id = hash;
        let mut cur_readid = &self.ids[id];
        let mut retry: usize = 0;

        while self.idx[id] != None 
            && 
            !cur_readid.is_none()
            &&
            cur_readid.as_ref().unwrap() != readid
        {
            retry = retry.saturating_add(1); // Saturating add is faster
            id = (hash.wrapping_add(retry)) % MAX_READS;
            if retry > 100_000 {
              assert!(self.entries < MAX_READS as usize, "More entries than MAX_READS!");
              println!("Error: More than 100,000 tries...");
            }
            cur_readid = &self.ids[id];
        }

        id as usize
    }

    #[inline]
    fn add_len(&mut self, readid: &str, len: &str) {
        // We pass len as &str because we don't want to do a conversion unless absolutely necessary...

        let id = self.get_id(readid);
        if let None = self.ids[id] {
            self.idx[id] = core::num::NonZeroU64::new(self.entries as u64);
            self.ids[self.entries] = Some(readid.to_string());
            self.lens[self.entries] = core::num::NonZeroU32::new(len.parse::<u32>().unwrap());

            self.entries += 1;
        }
    }
}

// JGG: TODO: Make compressed paf optional
pub fn parse_paf(prefix: String, 
                 threads: usize,
                 batch_size: usize,
                 filename: String) 
    -> HashMap<String, u32, BuildHasherDefault<XxHash64>> {

    let file = File::open(filename).expect("Unable to open file");
    let pb = ProgressBar::new(file.metadata().expect("Unable to get file metadata").len());
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})"));

    let output_channel:  Arc<ArrayQueue<ThreadCommand<HashMap<String, ThinVec<(u32, u32)>, BuildHasherDefault<XxHash64>>>>> = Arc::new(ArrayQueue::new(threads * 2));
    let process_channel: Arc<ArrayQueue<ThreadCommand<Vec<csv::StringRecord>>>> = Arc::new(ArrayQueue::new(threads * 2));

    let mut workers = Vec::with_capacity(threads);

    for i in 0..threads {
        let process_channel = Arc::clone(&process_channel);
        let output_channel = Arc::clone(&output_channel);
        let main_thread = std::thread::current();
        // let x = i.clone();
        
        let worker = thread::spawn(move || {
            let backoff = Backoff::new();
            let mut overlaps: HashMap<String, 
                                ThinVec<(u32, u32)>, 
                                BuildHasherDefault<XxHash64>> = Default::default();

            overlaps.reserve(batch_size);
            let mut count: usize = 0;

            loop {
                if let Ok(command) = process_channel.pop() {
                    if let ThreadCommand::Terminate = command {
                        // Do cleanup
                        if count > 0 {
                            let mut result = output_channel.push(ThreadCommand::Work(overlaps));
                            while let Err(PushError(overlaps)) = result {
                                backoff.snooze();
                                result = output_channel.push(overlaps);
                            }
                        }
                        return
                    }
                    let records = command.unwrap();

                    for result in records {
                        let id_a = result[0].to_string();
                        let id_b = result[5].to_string();

                        let ovl_a = (util::str2u32(&result[2]).unwrap(), util::str2u32(&result[3]).unwrap());
                        let ovl_b = (util::str2u32(&result[7]).unwrap(), util::str2u32(&result[8]).unwrap());

                        // JGG: TODO: Make into function...
                        let x = match overlaps.get_mut(&id_a) {
                            None => { overlaps.insert(id_a.clone(), ThinVec::new());
                                      overlaps.get_mut(&id_a).unwrap()
                                    },
                            Some(x) => x
                        };

                        x.push(ovl_a);

                        let x = match overlaps.get_mut(&id_b) {
                            None => { overlaps.insert(id_b.clone(), ThinVec::new());
                                      overlaps.get_mut(&id_b).unwrap()
                                    },
                            Some(x) => x
                        };

                        x.push(ovl_b);

                        count += 2;
                }
                    if count >= batch_size {
                        let mut result = output_channel.push(ThreadCommand::Work(overlaps));
                        while let Err(PushError(overlaps)) = result {
                            std::thread::park();
                            backoff.snooze();
                            result = output_channel.push(overlaps);
                        }
                        main_thread.unpark();
                        overlaps = Default::default();
                        overlaps.reserve(batch_size);
                    }
                } else {
                    main_thread.unpark();
                    std::thread::park();
                    backoff.snooze();
                }
            }});
        workers.push(worker);
    }

    let output_worker;

    {
        let main_thread = std::thread::current();
        let output_channel = Arc::clone(&output_channel);
        output_worker = thread::spawn(move || {
            let backoff = Backoff::new();
            let db = sled::Config::default()
                            .path(prefix.to_string())
                            .create_new(true)
                            .open().expect("Unable to open (or create) database!");

            loop {
                if let Ok(command) = output_channel.pop() {
                    if let ThreadCommand::Terminate = command {
                        db.flush().expect("Unable to perform final db flush");
                        return db
                    }

                    main_thread.unpark();

                    let mut output = command.unwrap();
                    let mut batch = Batch::default();

                    for (k, vs) in output.drain() {
                        let new_val: Vec<(u32, u32)> = match db.get(k.as_bytes()).expect("Unable to read from database") {
                            Some(x) => {
                                let mut orig: Vec<(u32, u32)> = bincode::deserialize(&x).expect("Unable to deserialize vector");
                                orig.extend(vs);
                                orig
                            },
                            None    => vs.to_vec(),
                        };

                        batch.insert(k.as_bytes(), bincode::serialize(&new_val).expect("Unable to serialize vector"));
                    }

                    db.apply_batch(batch).expect("Error applying batch update");
                    db.flush().expect("Unable to flush database");
                } else {
                    std::thread::park();
                    backoff.snooze();
                }
            }

        });
    }

    let now = Instant::now();
    println!("Starting to get reads2len");
    println!("R2L {}", now.elapsed().as_secs());

    let reader = BufReader::with_capacity(128 * 1024 * 1024, pb.wrap_read(file));
    let reader = BufReader::with_capacity(32 * 1024 * 1024, GzDecoder::new(reader));

    let mut readsidx = FastReadsIdx::new();

    let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(false)
            .from_reader(reader);

    let chunk_size = batch_size;
    let mut chunk = Vec::with_capacity(chunk_size);

    let backoff = Backoff::new();

    for record in reader.records() {
        let result = record.expect("Unable to read from CSV file");

        readsidx.add_len(&result[0], &result[1]);
        readsidx.add_len(&result[5], &result[6]);

        chunk.push(result);
        if chunk.len() == chunk_size {
            output_worker.thread().unpark();

            for x in &workers {
                x.thread().unpark();
            }

            let mut result = process_channel.push(ThreadCommand::Work(chunk));
            while let Err(PushError(chunk)) = result {
                pb.tick();
                pb.reset_eta();
                backoff.snooze();
                result = process_channel.push(chunk);
            }
            chunk = Vec::with_capacity(chunk_size);
        }
    }

    // Process any final pieces...
    if chunk.len() > 0 {
        let mut result = process_channel.push(ThreadCommand::Work(chunk));
        while let Err(PushError(chunk)) = result {
            backoff.snooze();
            result = process_channel.push(chunk);
        }
    }

    println!("Snoozing until no jobs left... {} currently left, {} in output", process_channel.len(), output_channel.len());
    while process_channel.len() > 0 {
        backoff.snooze();
        output_worker.thread().unpark();

        for x in &workers {
            x.thread().unpark();
        }
    }

    for _ in 0..workers.len() {
        let mut result = process_channel.push(ThreadCommand::Terminate);
        while let Err(PushError(x)) = result {
            result = process_channel.push(x);
        }
    }

    println!("Terminate sent, again snoozing until no jobs left... {} currently left, {} in output", process_channel.len(), output_channel.len());
    while process_channel.len() > 0 {
        backoff.snooze();
        output_worker.thread().unpark();

        for x in &workers {
            x.thread().unpark();
        }
    }

    for worker in workers {
        worker.join().expect("Unable to join worker");
    }

    pb.finish();

    println!("Worker threads joined...");
    println!("Process channel empty and closed. Snoozing until output channel empty. {} currently in queue", output_channel.len());
    let pb = ProgressBar::new(output_channel.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}]"));

    while output_channel.len() > 0 {
        pb.set_position(output_channel.len() as u64);
        backoff.snooze();
        pb.tick();
        pb.reset_eta();
        pb.reset_elapsed();
        println!("{} currently in queue", output_channel.len());
        output_worker.thread().unpark();
    }

    let mut result = output_channel.push(ThreadCommand::Terminate);
    while let Err(PushError(x)) = result {
        result = output_channel.push(x);
    }

    println!("Finished reads2len, converting to hashmap");
    println!("R2L {}", now.elapsed().as_secs());

    let result = readsidx.convert_to_hashmap();
    println!("Finished reads2len hashmap conversion");
    println!("R2L {}", now.elapsed().as_secs());

    let db = output_worker.join().expect("Unable to join output worker and get the db handle...");
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;

    extern crate tempfile;

    const PAF_FILE: &'static [u8] = b"1\t12000\t20\t4500\t-\t2\t10000\t5500\t10000\t4500\t4500\t255
1\t12000\t5500\t10000\t-\t3\t10000\t0\t4500\t4500\t4500\t255
";

    const M4_FILE: &'static [u8] = b"1 2 0.1 2 0 20 4500 12000 0 5500 10000 10000
1 3 0.1 2 0 5500 10000 12000 0 0 4500 10000
";

    #[test]
    fn paf() {
        let mut paf = tempfile::Builder::new()
            .suffix(".paf")
            .tempfile()
            .expect("Can't create tmpfile");

        paf.as_file_mut()
            .write_all(PAF_FILE)
            .expect("Error durring write of paf in temp file");

        let mut ovl = FullMemory::new();

        ovl.init(paf.into_temp_path().to_str().unwrap())
            .expect("Error in overlap init");

        assert_eq!(
            ["1".to_string(), "2".to_string(), "3".to_string(),]
                .iter()
                .cloned()
                .collect::<std::collections::HashSet<String>>(),
            ovl.get_reads()
        );

        assert_eq!(vec![(20, 4500), (5500, 10000)], ovl.overlap("1").unwrap());
        assert_eq!(vec![(5500, 10000)], ovl.overlap("2").unwrap());
        assert_eq!(vec![(0, 4500)], ovl.overlap("3").unwrap());
    }

    #[test]
    fn m4() {
        let mut m4 = tempfile::Builder::new()
            .suffix(".m4")
            .tempfile()
            .expect("Can't create tmpfile");

        m4.as_file_mut()
            .write_all(M4_FILE)
            .expect("Error durring write of paf in temp file");

        let mut ovl = FullMemory::new();

        ovl.init(m4.into_temp_path().to_str().unwrap())
            .expect("Error in overlap init");

        assert_eq!(
            ["1".to_string(), "2".to_string(), "3".to_string(),]
                .iter()
                .cloned()
                .collect::<std::collections::HashSet<String>>(),
            ovl.get_reads()
        );

        assert_eq!(vec![(20, 4500), (5500, 10000)], ovl.overlap("1").unwrap());
        assert_eq!(vec![(5500, 10000)], ovl.overlap("2").unwrap());
        assert_eq!(vec![(0, 4500)], ovl.overlap("3").unwrap());
    }
}
