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

use std::io::Write;

/* crate use */
use log::info;
use fs3::FileExt;
use anyhow::{anyhow, bail, Context, Result};
use crossbeam::utils::Backoff;

/* local use */
use crate::error;
use crate::reads2ovl;
use crate::util;

// Notes:
// JGG: twox hash is one of the fastest hashing algorithms out there
// Not always the fastest, but faster than the default
use twox_hash::XxHash64;
use std::hash::BuildHasherDefault;
use std::collections::HashMap;
use crossbeam::queue::{ArrayQueue, PushError};
use std::thread;
use std::sync::Arc;
use std::time::{Duration, Instant};

use thincollections::thin_vec::ThinVec;

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

impl ThreadCommand<(String, usize)> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> (String, usize) {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}


pub struct OnDisk {
    reads2ovl: HashMap<String, ThinVec<(u32, u32)>, BuildHasherDefault<XxHash64>>,
    reads2len: HashMap<String, usize, BuildHasherDefault<XxHash64>>,
    prefix: String,
    number_of_value: u64,
    buffer_size: u64,
}

impl OnDisk {
    pub fn new(prefix: String, buffer_size: u64) -> Self {
        let mut reads2ovl: HashMap<String, ThinVec<(u32, u32)>, BuildHasherDefault<XxHash64>> = Default::default();
        let mut reads2len: HashMap<String, usize, BuildHasherDefault<XxHash64>> = Default::default();

        reads2ovl.reserve(buffer_size as usize);
        reads2len.reserve(1024 * 1024 * 8);

        OnDisk {
            reads2ovl,
            reads2len,
            prefix,
            number_of_value: 0,
            buffer_size,
        }
    }

    // JGG: TODO: Enable snappy compression for temp files...
    fn clean_buffer(&mut self) -> Result<()> {
        info!(
            "Clear cache, number of value in cache is {}",
            self.number_of_value
        );

        for (key, values) in self.reads2ovl.iter_mut() {
            let prefix = self.prefix.clone();
            let mut output = std::io::BufWriter::with_capacity(1 * 1024 * 1024, OnDisk::create_yacrd_ovl_file(&prefix, key));

            for v in values.iter() {
                writeln!(output, "{},{}", v.0, v.1).with_context(|| {
                    error::Error::WritingError {
                        filename: format!("{}{}", &prefix, key),
                        format: util::FileType::YacrdOverlap,
                    }
                })?;
            }

            // JGG: Output will drop here, releasing the file lock

            values.clear();
        }

        self.number_of_value = 0;

        Ok(())
    }

    fn create_yacrd_ovl_file(prefix: &str, id: &str) -> std::fs::File {
        /* build path */
        let path = prefix_id2pathbuf(prefix, id);

        /* create parent directory if it's required */
        if let Some(parent_path) = path.parent() {
            std::fs::create_dir_all(parent_path).with_context(|| {
                error::Error::PathCreationError {
                    path: parent_path.to_path_buf(),
                }
            }).expect("Unable to create directory");
        }

        /* create file */
        let file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .with_context(|| error::Error::CantWriteFile {
                filename: path.to_string_lossy().to_string(),
            }).expect("Unable to open file!");

        // JGG: TODO: Should be a "try" with a timeout
        // But, should also be just fine and safe right now...
        file.lock_exclusive().expect("Unable to get file lock");
        file
    }
}

pub(crate) fn prefix_id2pathbuf(prefix: &str, id: &str) -> std::path::PathBuf {
    let mut path = std::path::PathBuf::from(prefix);
    path.push(id);
    path.set_extension("yovl");

    path
}

impl reads2ovl::Reads2Ovl for OnDisk {
    fn init(&mut self, filename: &str) -> Result<()> {
        self.sub_init(filename)?;

        self.clean_buffer()
            .with_context(|| anyhow!("Error during creation of tempory file"))?;
        self.number_of_value = 0;

        Ok(())
    }

    // JGG: This pull the overlap from the file, so the fact we basically
    // create a new reads2ovl instance for each thread seems to be ok...
    fn overlap(&self, id: &str) -> Result<Vec<(u32, u32)>> {
        let filename = format!("{}{}.yovl", self.prefix, id);
        if std::path::Path::new(&filename).exists() {
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(b',')
                .has_headers(false)
                .from_reader(std::io::BufReader::with_capacity(
                    4 * 1024 * 1024, // Files are usually small, so read them completely into memory
                    std::fs::File::open(&filename).with_context(|| error::Error::CantReadFile {
                        filename: filename.clone(),
                    })?,
                ));

            let mut ovls = Vec::new();
            for record in reader.records() {
                let result = record.with_context(|| error::Error::ReadingError {
                    filename: filename.clone(),
                    format: util::FileType::YacrdOverlap,
                })?;

                ovls.push((util::str2u32(&result[0])?, util::str2u32(&result[1])?));
            }

            Ok(ovls)
        } else {
            Ok(Vec::new())
        }
    }

    fn length(&self, id: &str) -> usize {
        *self.reads2len.get(&id.to_string()).unwrap_or(&0)
    }

    fn add_overlap(&mut self, id: String, ovl: (u32, u32)) -> Result<()> {
        // JGG: Rust's entry API is good, but slow...
        // self.reads2ovl.entry(id).or_insert_with(Vec::new).push(ovl);
        let x = match self.reads2ovl.get_mut(&id) {
            None => { self.reads2ovl.insert(id.clone(), ThinVec::new());
                      self.reads2ovl.get_mut(&id).unwrap()
                    },
            Some(x) => x
        };

        x.push(ovl);

        self.number_of_value += 1;

        if self.number_of_value >= self.buffer_size {
            self.clean_buffer()?;
        }
        Ok(())
    }

    // Don't think we actually need this here anymore...
    fn add_length(&mut self, id: String, length: usize) -> bool { // Return true if insert was against an empty entry, false if it already existed
        // self.reads2len.entry(id).or_insert(length);
        // Don't need to check for it, since re-adding will just overwrite the previous one...
        match self.reads2len.insert(id, length) {
            None => true,
            Some(_) => false
        }
    }

    fn get_reads(&self) -> std::collections::HashSet<String> {
        self.reads2len.keys().map(|x| x.to_string()).collect()
    }

    // JGG: TODO: Refactor to be more functional
    fn init_paf(&mut self, input: Box<dyn std::io::Read>) -> Result<()> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(false)
            .from_reader(input);

        let mut children = Vec::new();

        let channel: Arc<ArrayQueue<ThreadCommand<Vec<csv::StringRecord>>>> = Arc::new(ArrayQueue::new(512));

        let reads2len_channel: Arc<ArrayQueue<ThreadCommand<(String, usize)>>> = Arc::new(ArrayQueue::new(8192 * 12)); // Never want to block because of this...

        let now = Instant::now();

        let hashmap_child;

        {
            let reads2len_channel = Arc::clone(&reads2len_channel);
            hashmap_child = thread::spawn(move || {
                let mut reads2len: HashMap<String, usize, BuildHasherDefault<XxHash64>> = Default::default();
                reads2len.reserve(512 * 1024 * 1024); // Support 512 million reads before needing to re-allocate

                let backoff = Backoff::new();

                loop {
                    if let Ok(command) = reads2len_channel.pop() {
                        if let ThreadCommand::Terminate = command {
                            return reads2len;
                        }

                        let entry: (String, usize) = command.unwrap();

                        reads2len.insert(entry.0, entry.1);

                    } else {
                        backoff.snooze();
                    }
                }
            });
        }

        // JGG: TODO: Hardcoded 48 threads
        let threads = 48;
        for _ in 0..threads {
            let channel = Arc::clone(&channel);
            let reads2len_channel = Arc::clone(&reads2len_channel);
            // Need to make more of these
            // TODO: Maybe make individual functions or move this out of the
            // threading area?
            let mut r2o = OnDisk::new(self.prefix.clone(), self.buffer_size);
            let child = thread::spawn(move || {
                let backoff = Backoff::new();
                loop {
                    if let Ok(command) = channel.pop() {
                        if let ThreadCommand::Terminate = command {
                            // Flush out anything still in memory...
                            r2o.clean_buffer().expect("Unable to clean buffer");
                            return;
                        }

                        let records = command.unwrap();

                        for result in records {
                            
                            if result.len() < 9 {
                                panic!(error::Error::ReadingErrorNoFilename {
                                    format: util::FileType::Paf,
                                });
                            }
                
                            let id_a = result[0].to_string();
                            let id_b = result[5].to_string();
                
                            let len_a = util::str2usize(&result[1]).unwrap();
                            let len_b = util::str2usize(&result[6]).unwrap();
                
                            let ovl_a = (util::str2u32(&result[2]).unwrap(), util::str2u32(&result[3]).unwrap());
                            let ovl_b = (util::str2u32(&result[7]).unwrap(), util::str2u32(&result[8]).unwrap());
                
                            if r2o.add_length(id_a.clone(), len_a) {
                                let mut result = reads2len_channel.push(ThreadCommand::Work((id_a.clone(), len_a)));
                                while let Err(PushError(chunk)) = result {
                                    println!("Buffer full, waiting...");
                                    result = reads2len_channel.push(chunk);
                                }
                            }
                            
                            if r2o.add_length(id_b.clone(), len_b) {
                                let mut result = reads2len_channel.push(ThreadCommand::Work((id_b.clone(), len_b)));
                                while let Err(PushError(chunk)) = result {
                                    println!("Buffer full, waiting...");
                                    result = reads2len_channel.push(chunk);
                                }
                            }
                
                            r2o.add_overlap(id_a, ovl_a).unwrap();
                            r2o.add_overlap(id_b, ovl_b).unwrap();
                        }
                    } else {
                        backoff.snooze();
                    }
                }
            });

            children.push(child);
        }

        let chunk_size = 8192;
        let mut chunk = Vec::with_capacity(chunk_size);
        let backoff = Backoff::new();
        for record in reader.records() {
            let record = record.with_context(|| error::Error::ReadingErrorNoFilename {
                format: util::FileType::Paf,
            }).expect("Unable to read file properly...");
            chunk.push(record);
            if chunk.len() == chunk_size {
                let mut result = channel.push(ThreadCommand::Work(chunk));
                while let Err(PushError(chunk)) = result {
                    println!("Buffer full, waiting...");
                    backoff.spin();
                    result = channel.push(chunk);
                }
                chunk = Vec::with_capacity(chunk_size);
            }
        }

        println!("File read...");
        println!("2: {}", now.elapsed().as_secs());

        backoff.reset();
        for _ in 0..4 { // Very slight delay then issue terminate commands...
            backoff.spin();
        }

        for _ in 0..children.len() {
            let mut result = channel.push(ThreadCommand::Terminate);
            while let Err(PushError(chunk)) = result {
                result = channel.push(chunk);
            }
        }

        println!("Merging children... {} jobs left", channel.len());
        println!("2: {}", now.elapsed().as_secs());

        // JGG: Because of the multiple threads, need
        // to bring everything back into the main impl
        // TODO: Make more functional
        for child in children {
            child.join().expect("Unable to join one of the worker child threads");
        }

        println!("Children joined, now doing hashmap worker");

        self.reads2len = hashmap_child.join().expect("Unable to join hashmap child");

        println!("Hashmap retrieved succesfully");

        println!("Returning...");
        println!("2: {}", now.elapsed().as_secs());

        Ok(())
    }
}
