# JGG
This is a faster version of yacrd with modifications by Joseph Guhlin
Experimental - Do not use - blah blah blah blah

I've got a ~10Gb compressed PAF file to process for a very large genome

On 100k PAF matches -> 90 minutes (vanilla yacrd)
TwoxHash -> 85 minutes
Multi-threaded PAF reader -> 70 minutes (slower hardware too, had to change servers)


# Yet Another Chimeric Read Detector for long reads

Using all-against-all read mapping, yacrd performs:

1. computation of pile-up coverage for each read
2. detection of chimeras

Chimera detection is done as follows:

1. for each region where coverage is smaller or equal than `min_coverage` (default 0), yacrd creates a _bad region_.
2. if there is a _bad region_ that starts at a position strictly after the beginning of the read and ends strictly before the end of the read, the read is marked as `Chimeric`
3. if total _bad region_ length > 0.8 * read length, the read is marked as `NotCovered`
4. if a read isn't `Chimeric` or `NotCovered` is `NotBad`

## Rationale

Long read error-correction tools usually detect and also remove chimeras. But it is difficult to isolate or retrieve information from just this step.

DAStrim (from the [DASCRUBBER suite](https://github.com/thegenemyers/DASCRUBBER) does a similar job to yacrd but relies on a different mapping step, and uses different (likely more advanced) heuristics. Yacrd is simpler and easier to use.

This [repository](https://github.com/natir/yacrd-and-fpa-upstream-tools-for-lr-genome-assembly) contains a set of scripts to evaluate yacrd against other similar tools such as [DASCRUBBER](https://github.com/thegenemyers/DASCRUBBER/) and [miniscrub](https://bitbucket.org/berkeleylab/jgi-miniscrub) on real data sets.

## Input

Any set of long reads (PacBio, Nanopore, anything that can be given to [minimap2](https://github.com/lh3/minimap2)).
yacrd takes the resulting PAF (Pairwise Alignement Format) from minimap2 or BLASR m4 file from some other long reads overlapper as input.

## Requirements

- [Rust](https://www.rust-lang.org/) in stable channel
- libgz
- libbzip2
- liblzma

## Instalation

### With conda

yacrd is avaible in [bioconda channel](https://bioconda.github.io/)

if bioconda channel is setup you can run :

```
conda install yacrd
```

### From source

```
git clone https://github.com/natir/yacrd.git
cd yacrd
git checkout v0.6.1

cargo build
cargo test
cargo install --path .
```

## How to use Yacrd

### Find chimera

```
minimap2 reads.fq reads.fq > overlap.paf
yacrd -i overlap.paf -o reads.yacrd
```

### Post-detection operation

yacrd can perform some post-detection operation:

- filter: for sequence or overlap file, record with reads marked as Chimeric or NotCovered isn't write in output
- extract: for sequence or overlap file, record contains reads marked as Chimeric or NotCovered is write in output
- split: for sequence file bad region in middle of reads are removed, NotCovered read is removed
- scrubb: for sequence file all bad region are removed, NotCovered read is removed

```
minimap2 reads.fq reads.fq > mapping.paf
yacrd -i mapping.paf -o reads.yacrd filter -i reads.fasta -o reads.filter.fasta
yacrd -i mapping.paf -o reads.yacrd extract -i reads.fasta -o reads.extract.fasta
yacrd -i mapping.paf -o reads.yacrd split -i reads.fasta -o reads.split.fasta
yacrd -i mapping.paf -o reads.yacrd scrubb -i reads.fasta -o reads.scrubb.fasta
```

### Read scrubbing overlapping recommended parameter

For nanopore data, we recommend using minimap2 with all-vs-all nanopore preset with a maximal distance between seeds fixe to 500 (option `-g 500`) to generate overlap. We recommend to run yacrd with minimal coverage fixed to 4 (option `-c`) and minimal coverage of read fixed to 0.4 (option `-n`).

This is an exemple of how run a yacrd scrubbing:
```
minimap2 -x ava-ont -g 500 reads.fasta reads.fasta > overlap.paf
yacrd -i overlap.paf -o report.yacrd -c 4 -n 0.4 scrubb -i reads.fasta -o reads.scrubb.fasta
```

For pacbio P6-C4 data, we recommend to use minimap2 with all-vs-all pacbio preset with a maximal distance between seeds fixe to 800 (option `-g 800`) to generate overlap. We recommend to run yacrd with minimal coverage fixed to 4 (option `-c 4`) and minimal coverage of read fixed to 0.4 (option `-n 0.4`).

```
minimap2 -x ava-pb -g 800 reads.fasta reads.fasta > overlap.paf
yacrd -i overlap.paf -o report.yacrd -c 4 -n 0.4 scrubb -i reads.fasta -o reads.scrubb.fasta
```

For pacbio Sequel data, we recommend to use minimap2 with all-vs-all pacbio preset with a maximal distance between seeds fixe to 5000 (option `-g 5000`) to generate overlap. We recommand to run yacrd with minimal coverage fixed to 3 (option `-c 3`) and minimal coverage of read fixed to 0.4 (option `-n 0.4`).

```
minimap2 -x ava-pb -g 5000 reads.fasta reads.fasta > overlap.paf
yacrd -i overlap.paf -o report.yacrd -c 3 -n 0.4 scrubb -i reads.fasta -o reads.scrubb.fasta
```

### Important note

#### Extension

yacrd use extension to detect format file if your filename contains (anywhere):
- `.paf`: file is consider has minimap file
- `.m4`, `.mhap`: file is consider has blasr m4 file (mhap output)
- `.fa`, `.fasta`: file is consider has fasta file
- `.fq`, `.fastq`: file is consider has fastq file
- `.yacrd`: file is consider has yacrd output file

#### Compression

yacrd automatically detect file if is compress or not (gzip, bzip2 and lzma compression is available). For post-detection operation, if input is compressed output have the same compression format.

#### Use yacrd report as input

You can use yacrd report as input in place of overlap file, `ondisk` option are ignored if you use yarcd report has input.

## Output

```
type_of_read    id_in_mapping_file  length_of_read  length_of_gap,begin_pos_of_gap,end_pos_of_gap;length_of_gap,be…
```

### Example

```
NotCovered readA 4599    3782,0,3782
```

Here, readA doesn't have sufficient coverage, there is a zero-coverage region of length 3782bp between positions 0 and 3782.

```
Chimeric    readB   10452   862,1260,2122;3209,4319,7528
```

Here, readB is chimeric with 2 zero-coverage regions: one between bases 1260 and 2122, another between 4319 and 7528.

## Citation

If you use yacrd in your research, please cite the following publication:

```
Pierre Marijon, Rayan Chikhi, Jean-Stéphane Varré, yacrd and fpa: upstream tools for long-read genome assembly, Bioinformatics, btaa262, https://doi.org/10.1093/bioinformatics/btaa262
```

bibtex format:
```
@article {@article{Marijon_2020,
	doi = {10.1093/bioinformatics/btaa262},
	url = {https://doi.org/10.1093%2Fbioinformatics%2Fbtaa262},
	year = 2020,
	month = {apr},
	publisher = {Oxford University Press ({OUP})},
	author = {Pierre Marijon and Rayan Chikhi and Jean-St{\'{e}}phane Varr{\'{e}}},
	editor = {Inanc Birol},
	title = {yacrd and fpa: upstream tools for long-read genome assembly},
	journal = {Bioinformatics}
}
```
