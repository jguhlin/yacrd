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

// JGG: TODO: Benchmark Mimalloc vs. Jemalloc
// Note: Mimalloc is usually faster for me, but has issues compiling on Windows
use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/* Jemalloc is usually faster than rust default */
// extern crate jemallocator;
// #[global_allocator]
// static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

/* crate use */
use anyhow::{anyhow, Context, Result};
use clap::Clap;
extern crate snap;
extern crate fs3;
extern crate crossbeam;

use std::time::{Duration, Instant};

/* mod declaration*/
mod cli;
mod editor;
mod error;
mod reads2ovl;
mod stack;
mod util;

use cli::SubCommand as cmd;

fn main() -> Result<()> {
    env_logger::init();

    let now = Instant::now();

    let params = cli::Command::parse();

    /* Get bad region of reads */
    let mut reads2badregion: Box<dyn stack::BadPart> =
        if Some(util::FileType::Yacrd) == util::get_file_type(&params.input) {
            /* Read bad part from yacrd report */
            Box::new(stack::FromReport::new(&params.input)?)
        } else {
            /* Get bad part from overlap */
            let mut reads2ovl: Box<dyn reads2ovl::Reads2Ovl> = match params.ondisk.clone() {
                Some(prefix) => Box::new(reads2ovl::OnDisk::new(
                    prefix,
                    util::str2u64(&params.ondisk_buffer_size)?,
                )),
                None => Box::new(reads2ovl::FullMemory::new()),
            };

            reads2ovl.init(&params.input)?;

            println!("Finished reading in PAF file, converting to stack::FromOverlap");
            println!("{}", now.elapsed().as_secs());

            // JGG: TODO: Make multi-threaded
            Box::new(stack::FromOverlap::new(reads2ovl, params.coverage))
        };

    println!("{}", now.elapsed().as_secs());
    println!("Finished reading in PAF file...");
    println!("Writing out report...");

    /* Write report */
    let raw_out = Box::new(std::io::BufWriter::new(
        std::fs::File::create(&params.output).with_context(|| error::Error::CantWriteFile {
            filename: params.output.clone(),
        })?,
    ));

    let mut out = niffler::get_writer(
        raw_out,
        niffler::compression::Format::No,
        niffler::compression::Level::One,
    )?;

    // JGG: TODO: Make multi-threaded...
    for read in reads2badregion.get_reads() {
        let (bads, len) = reads2badregion.get_bad_part(&read)?;
        editor::report(&read, *len, bads, params.not_coverage, &mut out)
            .with_context(|| anyhow!("Filename: {}", &params.output))?;
    }

    println!("Report written...");
    println!("{}", now.elapsed().as_secs());
    println!("Running scrubb...");

    /* Run post operation on read or overlap */
    match params.subcmd {
        Some(cmd::Scrubb(s)) => editor::scrubbing(
            &s.input,
            &s.output,
            &mut *reads2badregion,
            params.not_coverage,
        )?,
        Some(cmd::Filter(f)) => editor::filter(
            &f.input,
            &f.output,
            &mut *reads2badregion,
            params.not_coverage,
        )?,
        Some(cmd::Extract(e)) => editor::extract(
            &e.input,
            &e.output,
            &mut *reads2badregion,
            params.not_coverage,
        )?,
        Some(cmd::Split(s)) => editor::split(
            &s.input,
            &s.output,
            &mut *reads2badregion,
            params.not_coverage,
        )?,
        None => (),
    };

    println!("{}", now.elapsed().as_secs());

    if let Some(prefix) = params.ondisk {
        for read in reads2badregion.get_reads() {
            let path = reads2ovl::ondisk::prefix_id2pathbuf(&prefix, &read);
            if path.is_file() {
                std::fs::remove_file(&path).with_context(|| anyhow!("We failed to remove file {:?}, yacrd finish analysis but temporary file isn't removed", path.clone()))?;
            }

            if let Some(parent_path) = path.parent() {
                if path.is_dir() {
                    std::fs::remove_dir_all(parent_path).with_context(|| {
                        error::Error::PathDestructionError {
                            path: parent_path.to_path_buf(),
                        }
                    })?;
                }
            }
        }
    }

    Ok(())
}
