# Prepare a truth set for NA24385/NA12878 mixture sequencing
# Downloads Genome in a Bottle calls for both genomes, cleans them
# merges the VCFs and extracts unique calls in NA12878 as target somatics.

configfile: 'giab_truthset-config.yaml'

# Download and prepare VCF files for two genomes using GiaB data

rule all:
    input:
        expand("calls/{genome}.vcf.gz.tbi", genome=config["genome"]),
        "na12878-na24385-somatic-truth.vcf.gz",
        "na12878-na24385-somatic-truth-regions.bed"

rule get_na12878_vcf:
    output: "calls/na12878-giab.vcf.gz"
    shell:
        r"""wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.2.2/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz"""

rule get_na12878_bed:
    output: "calls/na12878-giab.bed"
    shell:
        r"""wget -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.2.2/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed > {output}"""

rule get_na24385_vcf:
    output: "calls/na24385-giab.vcf.gz"
    shell:
        r"""wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.2.2/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf.vcf.gz"""

rule get_na24385_bed:
    output: "calls/na24385-giab.bed"
    shell:
        r"""wget -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.2.2/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf.bed > {output}"""

rule make_indices:
    input: "calls/{genome}.vcf.gz"
    output: "calls/{genome}.vcf.gz.tbi"
    shell: "tabix -f -p vcf {input}"

rule make_indices2:
    input: "work/{genome}-prep.vcf.gz"
    output: "work/{genome}-prep.vcf.gz.tbi"
    shell: "tabix -f -p vcf {input}"

# Merge files and extract somatic truth set

rule prep_na12878:
    input: "calls/na12878-giab.vcf.gz", "calls/na12878-giab.bed", "calls/na12878-giab.vcf.gz.tbi"
    output: "work/na12878-giab-prep.vcf.gz"
    shell:
        "bcftools view -R {input[1]} -o - {input[0]} | "
        "sed 's/INTEGRATION/NA12878/' | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "vt sort -m local -w 10000 -o {output} -"

rule prep_na24385:
    input: "calls/na24385-giab.vcf.gz", "calls/na24385-giab.bed", "calls/na24385-giab.vcf.gz.tbi"
    output: "work/na24385-giab-prep.vcf.gz"
    shell:
        "bcftools view -R {input[1]} -o - {input[0]} | "
        "sed 's/INTEGRATION/NA24385/' | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "vt sort -m local -w 10000 -o {output} -"

rule truth_regions:
    input: "calls/na12878-giab.bed", "calls/na24385-giab.bed"
    output: "na12878-na24385-somatic-truth-regions.bed"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} > {output}"

rule merge_both:
    input: "work/na12878-giab-prep.vcf.gz", "work/na24385-giab-prep.vcf.gz",
           "na12878-na24385-somatic-truth-regions.bed",
           "work/na12878-giab-prep.vcf.gz.tbi", "work/na24385-giab-prep.vcf.gz.tbi",
    output: "work/na12878-na24385.vcf.gz"
    shell: "bcftools merge {input[0]} {input[1]} -R {input[2]} -O z -o {output}"

rule extract_somatic:
    input: "work/na12878-na24385.vcf.gz"
    output: "work/na12878-na24385-somatic.vcf.gz"
    shell:
        "zcat {input} | "
        """vawk --header '(S$NA12878$GT!="0/0" && S$NA12878$GT!="./." """
        """ && S$NA12878$GT!=S$NA24385$GT) && (S$NA24385$GT=="0/0" || S$NA24385$GT=="./.")' | """
        "bgzip -c > {output}"

rule truth_set:
    """Prepare truth set from extracted somatic.

    - Converts homozygous ref calls to het, since they'll be hets when mixed with NA24385
      in tumor-like validation set.
    """
    input: "work/na12878-na24385-somatic.vcf.gz"
    output: "na12878-na24385-somatic-truth.vcf.gz"
    shell:
        "bcftools view {input} -s NA12878 -o - | "
        "vt sort -m local -w 10000 - | "
        "sed 's#1/1#0/1#' | bgzip -c > {output} && "
        "tabix -f -p vcf {output}"

