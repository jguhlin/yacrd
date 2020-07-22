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
use anyhow::{anyhow, Context, Result};
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
use std::time::{Instant};

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
    fn clean_buffer(&mut self, 
        output_channel: &Arc<ArrayQueue<ThreadCommand<Vec<(String, Vec<(u32, u32)>)>>>>) 
            -> Result<()> {
        info!(
            "Clear cache, number of value in cache is {}",
            self.number_of_value
        );

        for (key, values) in self.reads2ovl.iter_mut() {

            if values.len() > 0 {
                let prefix = self.prefix.clone();
                let mut file = OnDisk::create_yacrd_ovl_file(&prefix, key);
                let mut output = std::io::BufWriter::with_capacity(1 * 1024 * 1024, file);

                for v in values.iter() {
                    writeln!(output, "{},{}", v.0, v.1).with_context(|| {
                        error::Error::WritingError {
                            filename: format!("{}{}", &prefix, key),
                            format: util::FileType::YacrdOverlap,
                        }
                    })?;
                }
                file = output.into_inner().unwrap();
                file.unlock().expect("Unable to unlock file");
    
                values.clear();
            }
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

        let backoff = Backoff::new();

        while let Err(_) = file.try_lock_exclusive() {
            println!("{} still locked...", path.display());
            backoff.snooze();
        }

        // file.lock_exclusive().expect("Unable to get file lock")
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

        Ok(())
    }

    fn get_reads(&self) -> std::collections::HashSet<String> {
        self.reads2len.keys().map(|x| x.to_string()).collect()
    }

    fn init_paf(&mut self, input: Box<dyn std::io::Read>) -> Result<()> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(false)
            .from_reader(input);

        let mut children = Vec::new();

        let channel: Arc<ArrayQueue<ThreadCommand<Vec<csv::StringRecord>>>> = Arc::new(ArrayQueue::new(4096));
        let output_channel: Arc<ArrayQueue<ThreadCommand<Vec<(String, Vec<(u32, u32)>)>>>> = Arc::new(ArrayQueue::new(96));

        let now = Instant::now();

        // JGG: TODO: Hardcoded 48 threads
        let threads = 48;
        for i in 0..threads {
            let channel = Arc::clone(&channel);
            let output_channel = Arc::clone(&output_channel);
            // Need to make more of these
            // TODO: Maybe make individual functions or move this out of the
            // threading area?
            let mut r2o = OnDisk::new(self.prefix.clone(), self.buffer_size);
            let x = i.clone();
            let child = thread::spawn(move || {
                let backoff = Backoff::new();
                loop {
                    if let Ok(command) = channel.pop() {
                        if let ThreadCommand::Terminate = command {
                            // Flush out anything still in memory...
                            println!("Got terminate command! Cleaning... {}", x);
                            if r2o.number_of_value > 0 {
                                r2o.clean_buffer(&output_channel).expect("Unable to clean buffer");
                            }
                            println!("Returning {}", x);
                            return
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

                            // JGG: TODO: Just get length from the file directly, instead of using hashmap at all to keep track of it...
                            // Probably use sled as high performance backend...
                
                            let ovl_a = (util::str2u32(&result[2]).unwrap(), util::str2u32(&result[3]).unwrap());
                            let ovl_b = (util::str2u32(&result[7]).unwrap(), util::str2u32(&result[8]).unwrap());
                
                            r2o.add_overlap(id_a, ovl_a).unwrap();
                            r2o.add_overlap(id_b, ovl_b).unwrap();

                            if r2o.number_of_value >= r2o.buffer_size {
                                r2o.clean_buffer(&output_channel).expect("Unable to clean buffer");
                            }
                    
                        }
                    } else {
                        // Nothing to do, go ahead and clean buffer...
                        backoff.snooze();
                        if r2o.number_of_value > 5000 {
                            r2o.clean_buffer(&output_channel).expect("Unable to clean buffer");
                        }
                    }
                }
            });

            children.push(child);
        }

        let chunk_size = 4096;
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
                    println!("Chunks Buffer full, waiting...");
                    backoff.snooze();
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

        println!("Snoozing until no jobs left... {} currently left", channel.len());
        while channel.len() > 0 {
            backoff.snooze();
        }

        println!("issuing terminate to children");

        for _ in 0..children.len() {
            let mut result = channel.push(ThreadCommand::Terminate);
            backoff.spin();
            backoff.spin();
            backoff.spin();
            backoff.spin();
            backoff.spin();
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

        println!("Returning...");
        println!("2: {}", now.elapsed().as_secs());

        Ok(())
    }
}
