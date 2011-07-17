[Cascalog][1] query for assessing deep sequencing variation call
statistics with [Hadoop][3]. The goal is to examine problematic
call positions and assess how tweaking filtering parameters can help
to distinguish sequencing noise from real low frequency variations
in mixed populations.

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

        % lein deps
        % lein run :snp-data /directory/of/varation/data /directory/of/positions

To run on Hadoop:

        % lein deps
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
