# Somatic validation using rtg vcfeval
#
# Uses scripts from:
# https://github.com/hbc/projects/tree/master/giab_somatic
# to pre-process calls from multiple callers into a single output
# for comparison


configfile: 'somatic_concordance-config.yaml'

rule all:
    input: expand("{analysis}/somatic-rtg", analysis=config["analysis"])

rule prep_calls:
    input: expand("{analysis}/{calls}.vcf", analysis=config["analysis"], calls=config["calls"])
    output: "{analysis}/calls.vcf.gz"
    shell: "bgzip {input} -c > {output} && tabix -p vcf {output}"

rule combine_calls:
    input:  "{analysis}/calls.vcf.gz"
    output: "{analysis}/calls-cleaned.vcf"
    shell: "python {config[code][giab]}/scripts/combine_calls.py {input}"

rule prep_combined:
    input:  "{analysis}/calls-cleaned.vcf"
    output: "{analysis}/calls-cleaned.vcf.gz"
    shell: "bgzip {input} && tabix -p vcf {output}"

rule rtg_vcfeval:
    input: "{analysis}/calls-cleaned.vcf.gz", config["baseline"], config["ref"]["rtg"]
    output: "{analysis}/somatic-rtg"
    shell:
      "rm -rf {output[0]} && "
      "rtg vcfeval -c {input[0]} -b {input[1]} -t {input[2]} -o {output[0]}"
