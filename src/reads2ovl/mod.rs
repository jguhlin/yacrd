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

/* external crates */
use twox_hash::XxHash64;
use flate2::bufread::GzDecoder;
use t1ha::{t1ha0};
use indicatif::{HumanDuration, MultiProgress, ProgressBar, ProgressStyle};


pub const MAX_READS:  usize = 256_000_000;


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
                get_lengths_from_paf(filename)
            });
        }

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
pub fn get_lengths_from_paf(filename: String) -> HashMap<String, u32, BuildHasherDefault<XxHash64>> {
    let now = Instant::now();
    println!("Starting to get reads2len");
    println!("R2L {}", now.elapsed().as_secs());

    let reader = BufReader::with_capacity(
        64 * 1024 * 1024,
        File::open(filename).expect("Unable to open file"));

    let reader = BufReader::with_capacity(32 * 1024 * 1024, GzDecoder::new(reader));

    let mut readsidx = FastReadsIdx::new();

    let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(false)
            .from_reader(reader);

    for record in reader.records() {
        let result = record.expect("Unable to read from CSV file");

        readsidx.add_len(&result[0], &result[1]);
        readsidx.add_len(&result[5], &result[6]);
    }

    println!("Finished reads2len, converting to hashmap");
    println!("R2L {}", now.elapsed().as_secs());

    let result = readsidx.convert_to_hashmap();
    println!("Finished reads2len hashmap conversion");
    println!("R2L {}", now.elapsed().as_secs());
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
