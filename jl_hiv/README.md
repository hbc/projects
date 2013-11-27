## Variant calling in deeply sequenced viral populations

This set of scripts provides an automated pipeline for identifying mutations in
viral populations using Illumina deep sequencing. Coupled with
[our deep sequencing variant calling tool][1], it provides a full workflow to
move from Illumina fastq sequences to a VCF file of population mutations and
amino acid changes.

## Usage

Configuration files in [YAML][3] input format define all inputs to the process.
An [example configuration file][2] is a useful starting point. With this file,
the entire run process consists of a single commandline:

     python scripts/variant_identify.py <your_config.yaml>
     
This creates a `variation` directory containing files named
`raw_your-run-name-sort-realign.tsv` which has detailed statistics about each
position with aligned reads. These values feed directly into the
[variant calling framework][1].

## What does it do?

The build script performs the following steps to prepare for variant calling:

- Collapses the input fastq reads into unique reads. At high sequencing depth,
  we expect extensive read duplication, and this step avoids uncessary overhead of
  aligning identical reads multiple times.

- Aligns collapsed fastq reads to reference genome. This handles ambiguous reference
  genomes with IUPAC characters, which is useful for error matching in viral
  populations with known variant regions.
  
- Re-aligns reads, avoiding inconsistent and incorrect alignments due to indels.

- Summarizes unique reads at each position with read quality score, alignment
  quality score and percent representation of the k-mer surrounding region.
  These metrics feed directly into variant calling.

## Installation

The pipeline leverages these freely available tools:

- [novoalign][5] -- alignment to the reference genome
- [Picard][8] -- Manipulation of BAM alignment files
- [GATK][6] -- re-alignment of reads around indels
- [khmer][7] -- count k-mer regions surrounding each variant

The [CloudBioLinux][4] project provides automated installation scripts with all
of these dependencies.

Following installation of these, run:

    $ python setup.py build
    $ sudo python setup.py install
    
to install required Python libraries.

[1]: https://github.com/hbc/projects/tree/master/snp-assess
[2]: https://github.com/hbc/projects/blob/master/jl_hiv/config/20120111.yaml
[3]: https://en.wikipedia.org/wiki/YAML
[4]: http://cloudbiolinux.org
[5]: http://www.novocraft.com/main/index.php
[6]: http://www.broadinstitute.org/gatk/
[7]: https://github.com/ged-lab/khmer
[8]: http://picard.sourceforge.net/

