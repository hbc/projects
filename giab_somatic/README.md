## Somatic truth sets from Genome in a Bottle samples

A mixture of two [Genome in a Bottle](https://github.com/genome-in-a-bottle)
samples -- NA12878 and NA24385 -- to emulate a somatic-like tumor-normal set.
Known calls from these two samples can be used to estimate true and false
positives from somatic variant callers. Input BAMs and a complete description
of the experimental setup are available from the Genome in a Bottle FTP site:

ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/

`giab_truthset.snakefile` has the full commands for creating the truth sets.

NA12878 and NA24385 are from Genome in a Bottle v3.3.2 calls.

The somatic truth sets contain two annotations:

- FREQ: Specifying the expected frequency of the mutation: 30% and 15% for
  NA12878 homozygotes and heterozygotes from the 30% NA12878 70% NA24385
  tumor-like mixture.
- SOMTYPE: A classification for the somatic type (high_freq_somatic and
  mod_freq_somatic) useful in training.

## Build 37

* Somatic VCF truth set
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-GRCh37-truth.vcf.gz
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-GRCh37-truth.vcf.gz.tbi
* Callable BED regions
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-GRCh37-truth-regions.bed
* Germline VCF truth set
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-germline-GRCh37-truth.vcf.gz
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-germline-GRCh37-truth.vcf.gz.tbi

# Build 38

* Somatic VCF truth set
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-hg38-truth.vcf.gz
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-hg38-truth.vcf.gz.tbi
* Callable BED regions
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-somatic-hg38-truth-regions.bed
* Germline VCF truth set
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-germline-hg38-truth.vcf.gz
  - https://s3.amazonaws.com/bcbio/giab/NA12878-NA24385/2018-07-05/na12878-na24385-germline-hg38-truth.vcf.gz.tbi
