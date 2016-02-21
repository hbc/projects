# Germline validation using rtg vcfeval and post-comparison
# extraction of false positives and negatives to characterize
# concordance and discordance.

configfile: 'germline_concordance-config.yaml'

rule all:
    input: "germline-rtg", "compare/fn-overlap.vcf.gz", "compare/fp-overlap.vcf.gz",
           "compare/fn-missing.vcf.gz", "compare/fp-missing.vcf.gz"

rule rtg_vcfeval:
    input: config["calls"], config["baseline"], config["regions"], config["ref"]["rtg"]
    output: "germline-rtg", config["eval"]["fn"], config["eval"]["fp"]
    shell:
      "rmdir germline-rtg && "
      "rtg vcfeval -c {input[0]} -b {input[1]} --bed-regions {input[2]} "
      "-t {input[3]} --sample={config[sample]} -o {output[0]}"

rule extra_fns_missing:
    input: config["calls"], config["eval"]["fn"]
    output: "compare/fn-missing.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=1 -w 2 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fns:
    input: config["calls"], config["eval"]["fn"]
    output: "compare/fn-overlap.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=2 -w 1 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fps_missing:
    input:  "compare/baseline-sample.vcf.gz", config["eval"]["fp"]
    output: "compare/fp-missing.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=1 -w 2 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule extract_fps:
    input:  "compare/baseline-sample.vcf.gz", config["eval"]["fp"]
    output: "compare/fp-overlap.vcf.gz"
    shell:
      "bcftools isec {input[0]} {input[1]} -n=2 -w 1 -o {output} -O z "
      "&& tabix -f -p vcf {output}"

rule baseline_sample:
    input: config["baseline"]
    output: "compare/baseline-sample.vcf.gz"
    shell:
      "bcftools view {input} -s {config[sample]} -o {output} -O z "
      "&& tabix -f -p vcf {output}"

