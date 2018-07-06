# Prepare a truth set for NA24385/NA12878 mixture sequencing
# Downloads Genome in a Bottle calls for both genomes, cleans them
# merges the VCFs and extracts unique calls in NA12878 as target somatics.

configfile: 'giab_truthset-hg38-config.yaml'

# Download and prepare VCF files for two genomes using GiaB data

rule all:
    input:
        expand("work-{build}/calls/{genome}.vcf.gz.tbi", genome=config["genome"], build=config["build"]),
        expand("na12878-na24385-somatic-{build}-truth.vcf.gz", build=config["build"]),
        expand("na12878-na24385-somatic-{build}-truth-regions.bed", build=config["build"]),
        expand("na12878-na24385-germline-{build}-truth.vcf.gz", build=config["build"])

rule create_workdir:
    output: expand("work-{build}", build=config["build"])
    shell:
        "mkdir -p {output}"

rule get_na12878_vcf:
    output: "work-{build}/calls/na12878-giab.vcf.gz"
    shell:
        r"""wget -O {output} {config[na12878_vcf]}"""

rule get_na12878_bed:
    output: "work-{build}/calls/na12878-giab.bed"
    shell:
        r"""wget -O {output} {config[na12878_bed]} """

rule get_na24385_vcf:
    output: "work-{build}/calls/na24385-giab.vcf.gz"
    shell:
        r"""wget -O {output} {config[na24385_vcf]}"""

rule get_na24385_bed:
    output: "work-{build}/calls/na24385-giab.bed"
    shell:
        r"""wget -O - {config[na24385_bed]} > {output}"""

rule make_indices:
    input: "work-{build}/calls/{genome}.vcf.gz"
    output: "work-{build}/calls/{genome}.vcf.gz.tbi"
    shell: "tabix -f -p vcf {input}"

rule make_indices2:
    input: "work-{build}/{genome}-prep.vcf.gz"
    output: "work-{build}/{genome}-prep.vcf.gz.tbi"
    shell: "tabix -f -p vcf {input}"

# Merge files and extract somatic truth set

rule prep_na12878:
    input: "work-{build}/calls/na12878-giab.vcf.gz", "work-{build}/calls/na12878-giab.bed", "work-{build}/calls/na12878-giab.vcf.gz.tbi"
    output: "work-{build}/na12878-giab-prep.vcf.gz"
    shell:
        "bcftools view -R {input[1]} -o - {input[0]} | "
        "sed 's/INTEGRATION/NA12878/' | sed 's/HG001/NA12878/' | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "vt sort -m local -w 10000 -o {output} -"

rule prep_na24385:
    input: "work-{build}/calls/na24385-giab.vcf.gz", "work-{build}/calls/na24385-giab.bed", "work-{build}/calls/na24385-giab.vcf.gz.tbi"
    output: "work-{build}/na24385-giab-prep.vcf.gz"
    shell:
        "bcftools view -R {input[1]} -o - {input[0]} | "
        "sed 's/INTEGRATION/NA24385/' | sed 's/HG002/NA24385/' | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "vt sort -m local -w 10000 -o {output} -"

rule truth_regions:
    input: "work-{build}/calls/na12878-giab.bed", "work-{build}/calls/na24385-giab.bed"
    output: "na12878-na24385-somatic-{build}-truth-regions.bed"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} > {output}"

rule merge_both:
    input: "work-{build}/na12878-giab-prep.vcf.gz", "work-{build}/na24385-giab-prep.vcf.gz",
           "na12878-na24385-somatic-{build}-truth-regions.bed",
           "work-{build}/na12878-giab-prep.vcf.gz.tbi", "work-{build}/na24385-giab-prep.vcf.gz.tbi",
    output: "work-{build}/na12878-na24385.vcf.gz"
    shell: "bcftools merge {input[0]} {input[1]} -R {input[2]} -O z -o {output}"

rule extract_somatic:
    input: "work-{build}/na12878-na24385.vcf.gz"
    output: "work-{build}/na12878-na24385-somatic.vcf.gz"
    shell:
        "zcat {input} | "
        """vawk --header '(S$NA12878$GT!="0/0" && S$NA12878$GT!="./." """
        """ && S$NA12878$GT!=S$NA24385$GT) && (S$NA24385$GT=="0/0" || S$NA24385$GT=="./.")' | """
        "bgzip -c > {output}"

rule extract_germline:
    input: "work-{build}/na12878-na24385.vcf.gz"
    output: "na12878-na24385-germline-{build}-truth.vcf.gz"
    shell:
        "zcat {input} | "
        """vawk --header '(S$NA12878$GT!="0/0" && S$NA12878$GT!="./." """
        """ && S$NA24385$GT!="0/0" && S$NA24385$GT!="./.")' | """
	"bcftools view -s NA12878 -o - | "
        "vt sort -m local -w 10000 - | "
        "sed 's#1/1#0/1#' | bgzip -c > {output} && "
	"tabix -f -p vcf {output}"

rule create_annotation:
    input: "work-{build}/na12878-na24385-somatic.vcf.gz"
    output: "work-{build}/na12878-na24385-somatic-annotations.bed.gz"
    shell:
      "zcat {input} | "
      """vawk '{{if (S$NA12878$GT=="0/1") {{print $1,$2-1,$2,"mod_freq_somatic",15}} """
      """else {{print $1,$2-1,$2,"high_freq_somatic",30}} }}' """
      "| bgzip -c > {output} && "
      "tabix -f -p bed {output}"

rule truth_set:
    """Prepare truth set from extracted somatic.

    - Converts homozygous ref calls to het, since they'll be hets when mixed with NA24385
      in tumor-like validation set.
    """
    input: "work-{build}/na12878-na24385-somatic.vcf.gz", "work-{build}/na12878-na24385-somatic-annotations.bed.gz", "work-{build}"
    output: "na12878-na24385-somatic-{build}-truth.vcf.gz"
    shell:
        "bcftools view {input[0]} -s NA12878 -o - | "
        "vt sort -m local -w 10000 - | "
	"vcfanno -base-path {input[2]} giab_truthset-annotate.conf /dev/stdin | "
        "sed 's#1/1#0/1#' | bgzip -c > {output} && "
        "tabix -f -p vcf {output}"
