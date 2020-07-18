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

/* local use */
use crate::reads2ovl;
use crate::util;
use crate::error;

pub struct FullMemory {
    reads2ovl: std::collections::HashMap<String, (Vec<(u32, u32)>, usize)>,
    no_overlap: Vec<(u32, u32)>,
}

impl FullMemory {
    pub fn new() -> Self {
        FullMemory {
            reads2ovl: std::collections::HashMap::new(),
            no_overlap: Vec::new(),
        }
    }
}

impl reads2ovl::Reads2Ovl for FullMemory {
    fn overlap(&self, id: &str) -> Result<Vec<(u32, u32)>> {
        if let Some((vec, _)) = self.reads2ovl.get(&id.to_string()) {
            Ok(vec.to_vec())
        } else {
            Ok(self.no_overlap.to_vec())
        }
    }

    fn length(&self, id: &str) -> usize {
        if let Some((_, len)) = self.reads2ovl.get(&id.to_string()) {
            *len
        } else {
            0
        }
    }

    fn add_overlap(&mut self, id: String, ovl: (u32, u32)) -> Result<()> {
        self.reads2ovl
            .entry(id)
            .or_insert((Vec::new(), 0))
            .0
            .push(ovl);

        Ok(())
    }

    fn init_paf(&mut self, input: Box<dyn std::io::Read>) -> Result<()> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(false)
            .from_reader(input);

        for record in reader.records() {
            let result = record.with_context(|| error::Error::ReadingErrorNoFilename {
                format: util::FileType::Paf,
            })?;

            if result.len() < 9 {
                bail!(error::Error::ReadingErrorNoFilename {
                    format: util::FileType::Paf,
                });
            }

            let id_a = result[0].to_string();
            let id_b = result[5].to_string();

            let len_a = util::str2usize(&result[1])?;
            let len_b = util::str2usize(&result[6])?;

            let ovl_a = (util::str2u32(&result[2])?, util::str2u32(&result[3])?);
            let ovl_b = (util::str2u32(&result[7])?, util::str2u32(&result[8])?);

            self.add_length(id_a.clone(), len_a);
            self.add_length(id_b.clone(), len_b);

            self.add_overlap(id_a, ovl_a)?;
            self.add_overlap(id_b, ovl_b)?;
        }

        Ok(())
    }

    fn add_length(&mut self, id: String, length: usize) {
        self.reads2ovl.entry(id).or_insert((Vec::new(), 0)).1 = length;
    }

    fn get_reads(&self) -> std::collections::HashSet<String> {
        self.reads2ovl.keys().map(|x| x.to_string()).collect()
    }
}
