# Variant calling in deeply sequenced viral populations

This library provides variant calling for viral populations provided Illumina
inputs. Given read statistics calculated by
[our alignment and quality assessment pipeline][0-0], it identifies low
frequency mutations (> 0.1%) in the population by separating real variations
from sequencing errors. It uses an input control library with known low
frequency variations to train a random forest classifier that filters sequencing
errors.

The detailed steps in variant calling are:

- Use a control population with known low frequency variations to identify
  positions with real variations and positions with false positives caused by
  sequencing errors.
  
- Extract classification metrics from the raw read statistics:

   - Read quality. Phred scaled quality scores given to the individual called
     base by the sequencer.
   - Mapping score. The Novoalign score assigned to a sequencing read positioned 
     in the genome.
   - K-mer abundance. A metric describing the uniqueness of the 13bp region
     surrounding the read. Rare k-mers with low frequencies are more likely to be
     sequencing errors.

- Train a [random forest ensemble classifier][0-4] using the classification
  metrics and true and false positive examples identified.
  
- Use this classifier to identify variants in each of the sample populations.
  Following filtering with the random forest classifier, remaining low frequency
  changes are real variants.
  
- Annotate variations with amino acid changes and map these changes to known
  drug resistance mutations.

- Report resulting variations, annotations and metrics in a standard
  [VCF format][0-3] file.

## Usage

To get started, you need four things:

- A processing directory created by the [alignment pipeline][0-0], which will
  have a `variation` subdirectory containing detailed read information.
- A run configuration file with detailed information about the experiments to
  process: [example configuration file][0-1].
- A parameter configuration file specifying the filtering and calling approaches
  to use: [example configuration file][0-2].
- The [Leiningen][2] build tool, which will install required dependencies.

Then run the entire variant calling process with a single command:

    $ lein snp-call <run config YAML> <parameter config YAML> <work directory>
    
Each `variation/raw_*tsv` file of read statistics will have an associated
`variation/raw_*calls-protein-annotate.vcf` file containing variants and amino
acid changes in [VCF format][0-3].

[0-0]: https://github.com/hbc/projects/tree/master/jl_hiv
[0-1]: https://github.com/hbc/projects/blob/master/snp-assess/config/20120111-variant.yaml
[0-2]: https://github.com/hbc/projects/blob/master/snp-assess/config/20120111-classification.yaml
[0-3]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[0-4]: https://en.wikipedia.org/wiki/Random_forest

# Parallelized Hadoop utilities

In addition to variant calling, this library provides parallel utilities for
mining and evaluating next-generation sequencing data. It uses [Cascalog][1]
queries for assessing deep sequencing variation call statistics with
[Hadoop][3]. The goal is to examine problematic call positions and assess how
tweaking filtering parameters can help to distinguish sequencing noise from real
low frequency variations in mixed populations.

The test directory contains example variation data and positions. The
variation data file is a combination of positions, individual read
bases and quality, kmer frequency and alignment scores associated with
those read bases. Each of these can be tweaked to improve calling
accuracy.

Requires:

* [Leiningen][2]
* [Hadoop][3], only if you'd like to distribute the queries across a
  Hadoop cluster. This can also be run without Hadoop installed.

To run standalone:

        % lein snp-data /directory/of/varation/data /directory/of/positions

To run on Hadoop:

        % lein uberjar
        % hadoop fs -mkdir /tmp/snp-assess/data
        % hadoop fs -mkdir /tmp/snp-assess/positions
        % hadoop fs -put your_variation_data.tsv /tmp/snp-assess/data
        % hadoop fs -put positions_of_interest.tsv /tmp/snp-assess/positions
        % hadoop jar snp-assess-0.0.1-SNAPSHOT-standalone.jar
                     snp_assess.core /tmp/snp-assess/data /tmp/snp-assess/positions

Outputs summarize counts, mean kmer frequency, mean quality and mean
alignment score at each position and base:

        HXB2_IUPAC_93-5 951 A     1 8.8e-05 32.0 66.0
        HXB2_IUPAC_93-5 951 C     1 2.2e-06 8.0 50.0
        HXB2_IUPAC_93-5 951 G    83 5.0e-02 28.4 137.8
        HXB2_IUPAC_93-5 951 T     3 2.0e-04 24.7 55.7
        HXB2_IUPAC_93-5 953 A    10 1.6e-02 23.1 175.5
        HXB2_IUPAC_93-5 953 C     1 1.4e-04 28.0 53.0
        HXB2_IUPAC_93-5 953 G   126 9.9e-02 19.6 59.1
        HXB2_IUPAC_93-5 953 T    14 7.5e-04 10.1 61.4

[1]: http://github.com/nathanmarz/cascalog
[2]: https://github.com/technomancy/leiningen#readme
[3]: http://www.cloudera.com/hadoop/

## Off-target variation

The off-target limit of detection assesses what range of variation can
be accurately detected by looking at the distribution of minority
variants in positions without any expected variation. The approach to
plot this is:

* Start with list of positions without expected variation and
  expected base at that position.
* Calculate frequency of non-expected variations at each position.
* Calculate frequency of non-expected variations after filtration.
* Collect counts of off-target frequencies within bins.
* Plot as histogram

## Required coverage

by sample multiplexing, we would like to be able to reduce coverage to
maximize samples processed but still be able to detect variants of
interest. To do this, plot percentage of variants detected at
different amounts of coverage:

* Start with list of minority variant positions, expected base
  and frequency of variation.
* Define function to determine if variant is detected.
* Reduce coverage by randomly removing reads. Repeat until we lose
  ability to detect. Can repeat multiple times to sample.
* Report minimum detection coverage.
* Plot distribution of minimum detection coverage at each frequency
  of interest.

Another useful graph is a histogram of coverage at each position,
to assess the coverage distribution for a given number of total
reads.

## Linear classifier for evaluating multi-variant positions

Expected variable positions and raw data feed a linear classifier, implemented
in [Weka][4] and wrapped for Clojure by [clj-ml][5]:

    lein snp-classify raw-data-file position-file work-directory

`raw-data-file` is a tab-delimited file of calls plus associated quality
metrics; `test/data/raw/raw_variations.tsv`. `position-file` contains known
expected variations and their frequencies: `test/data/coverage_pos/pos.tsv`.
`work-directory` is the directory to write output information to. The
`classifier` subdirectory contains the build classifier and raw details about
detected variations as a YAML dump. Configuring the classifier by adjusting
`config.clj`; this could take a YAML input at the commandline as well.

To summarize read calls in a region:

    lein snp-classify-eval raw-data-file work-directory

[4]: http://www.cs.waikato.ac.nz/~ml/weka/
[5]: https://github.com/leadtune/clj-ml
