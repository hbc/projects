# Germline validation using rtg vcfeval and post-comparison
# extraction of false positives and negatives to characterize
# concordance and discordance.
#
# bwa version differences (0.7.12 versus 0.7.5a) explain the majority of 
# discordance. We see two classes of issues:
# 
# - Filter differences: these are variants identified in both callsets that fall
#   in lower mapping quality regions. Depending on the bwa version, they get
#   assigned different mapping qualities, resulting in triggering or not
#   triggering the GATK MappingQuality and LowQualityDepth hard filters. The
#   difference in filter categorizes the majority of the false positives and
#   some of the false negatives.
# 
# - Low complexity regions: these are variants that fall into low complexity
#   regions identified by Heng Li
#   (http://bioinformatics.oxfordjournals.org/content/30/20/2843). These are 2%
#   of the genome but contribute to a large percentage of concordance
#   differences (section 3.2 in the paper). It looks like bwa became more strict
#   aligning in those regions, as a large percentage of the remaining
#   false negatives fall in these regions.
# 
# Eliminating discordance from filter differences and low complexity regions this
# brings the comparison in line with expectations:
# 
# False positives: 69752 - 41847 (filter differences)
# False negatives: 126599 - 42607 (filter differences) - 52709 (low complexity regions)
# 
# Threshold    True-pos  False-pos  False-neg  Precision  Sensitivity  F-measure
# ----------------------------------------------------------------------------
# bwa adjusted  5345611      27905     31283     0.9948      0.9942     0.9945
# original      5345611      69752     126599    0.9871      0.9769     0.9820


configfile: 'germline_concordance-config.yaml'

rule all:
    input:
      expand("{analysis}/counts-fn_missing.txt", analysis=config["analysis"]),
      expand("{analysis}/counts-discordant.txt", analysis=config["analysis"])

rule count_fn_missing:
    """Count missing cases in LCRs due to bwa version differences."""
    input: "{analysis}/compare/fn-missing.vcf.gz"
    output: "{analysis}/counts-fn_missing.txt"
    shell:
      "bedtools intersect -wa -a {input} -b {config[ref][lcr]} | wc -l > {output}"

rule count_discordant_overlap:
    """Count details on missing cases with overlap between callers."""
    input: "{analysis}/compare/fn-overlap.vcf.gz", "{analysis}/compare/fp-overlap.vcf.gz"
    output: "{analysis}/counts-discordant.txt"
    shell:
      "zgrep -c -v ^# {input} > {output}"

rule prep_calls:
    output: "{analysis}/calls.vcf.gz"
    shell: "wget -O - {config[url]} | bgzip -c > {output} && tabix -p vcf {output}"

rule rtg_vcfeval:
    input: "{analysis}/calls.vcf.gz", "{analysis}/compare/baseline-sample.vcf.gz",
           config["regions"], config["ref"]["rtg"]
    output: "{analysis}/germline-rtg", "{analysis}/germline-rtg/fn.vcf.gz", "{analysis}/germline-rtg/fp.vcf.gz"
    shell:
      "rm -rf {output[0]} && "
      "rtg vcfeval -c {input[0]} -b {input[1]} --bed-regions {input[2]} "
      "-t {input[3]} -o {output[0]}"

rule extra_fns_missing:
    input: "{analysis}/calls.vcf.gz", "{analysis}/germline-rtg/fn.vcf.gz"
    output: "{analysis}/compare/fn-missing.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=1 -w 2 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fns:
    input: "{analysis}/calls.vcf.gz", "{analysis}/germline-rtg/fn.vcf.gz"
    output: "{analysis}/compare/fn-overlap.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=2 -w 1 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fps_missing:
    input:  "{analysis}/compare/baseline-sample.vcf.gz", "{analysis}/germline-rtg/fp.vcf.gz"
    output: "{analysis}/compare/fp-missing.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=1 -w 2 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fps:
    input:  "{analysis}/compare/baseline-sample.vcf.gz", "{analysis}/germline-rtg/fp.vcf.gz"
    output: "{analysis}/compare/fp-overlap.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=2 -w 1 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule baseline_sample:
    input: config["baseline"]
    output: "{analysis}/compare/baseline-sample.vcf.gz"
    shell:
      "bcftools view {input} -s {config[sample]} -o {output} -O z "
      "&& tabix -f -p vcf {output}"

