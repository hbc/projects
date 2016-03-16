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
        r"""wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz"""

rule get_na12878_bed:
    output: "calls/na12878-giab.bed"
    shell:
        r"""wget -O - ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz | gunzip -c > {output}"""

rule get_na24385_vcf_rtg:
    output: "calls/na24385-rtg.vcf.gz"
    shell:
        "wget -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/Rutgers_IlluminaHiSeq300X_rtg_11052015/rtg_allCallsV2.vcf.gz | "
        "gunzip | "
        "sed s/HG002/NA24385/ | "
        "bcftools view --exclude-uncalled -f 'PASS,.' -s NA24385 | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "bgzip -c > {output} "
        "&& tabix -p vcf -f {output}"

rule get_na24385_vcf_cgi:
    """Retrieve CGI calls for NA24385, removing partially called, uncalled and filtered variants.

    Also converts genotypes into non-phased normalized format.
    """
    output: "calls/na24385-cgi.vcf.gz"
    shell:
        "wget -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/CompleteGenomics_newLFR_CGAtools_06122015/vcfBeta-GS000039526-ASM-NA24385_WellCount.vcf | "
        "sed s/GS000039526-ASM/NA24385/ |"
        r"grep -v '0/\.:' | grep -v '0|\.:' | grep -v '1/\.:' | grep -v '1|\.:' | "
        r"grep -v '\.|0:' | grep -v '\.|1:' | "
        r"sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        r"sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        r"sed 's#2|1#1/2#' | sed 's#1|2#1/2#' | sed 's#2/1#1/2#' | "
        """bcftools view --exclude-uncalled -i 'FORMAT/FT="PASS" || FORMAT/FT="."' -O z -o {output} """
        "&& tabix -f -p vcf {output}"

rule make_indices:
    input: "calls/{genome}.vcf.gz"
    output: "calls/{genome}.vcf.gz.tbi"
    shell: "tabix -f -p vcf {input}"

# Merge files and extract somatic truth set

rule prep_na12878:
    input: "calls/na12878-giab.vcf.gz", "calls/na12878-giab.bed", "inputs/fix-na12878-header.txt", "calls/na12878-giab.vcf.gz.tbi"
    output: "calls/na12878-giab-prep.vcf.gz"
    shell:
        "bcftools annotate -h {input[2]} -R {input[1]} -o - {input[0]} | "
        "sed 's#0|0#0/0#' | sed 's#0|1#0/1#' | sed 's#1|1#1/1#' | "
        "sed 's#1|0#0/1#' | sed 's#1/0#0/1#' | "
        "vt sort -m local -w 10000 -o {output} -"

rule union_na24385:
    input: "calls/na24385-rtg.vcf.gz", "calls/na24385-cgi.vcf.gz", config["ref"]["seq"]
    output: "work/na24385-rtg-cgi.vcf.gz"
    shell:
        "bcbio-variation-recall ensemble -n 1 --names rtg,cgi "
        "{output} {input[2]} {input[0]} {input[1]}"

rule merge_both:
    input: "calls/na12878-giab-prep.vcf.gz", "work/na24385-rtg-cgi.vcf.gz", "calls/na12878-giab.bed",
           "calls/na12878-giab-prep.vcf.gz.tbi"
    output: "work/na12878-na24385.vcf.gz"
    shell: "bcftools merge {input[0]} {input[1]} -R {input[2]} -O z -o {output}"

rule extract_somatic:
    input: "work/na12878-na24385.vcf.gz"
    output: "na12878-na24385-somatic.vcf.gz"
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
    input: "na12878-na24385-somatic.vcf.gz"
    output: "na12878-na24385-somatic-truth.vcf.gz"
    shell:
        "bcftools view {input} -s NA12878 -o - | "
        "vt sort -m local -w 10000 - | "
        "sed 's#1/1#0/1#' | bgzip -c > {output} && "
        "tabix -f -p vcf {output}"

rule callable_regions_na24385:
    output: "calls/na24385-callable.bed"
    shell:
        r"""wget -O - https://s3.amazonaws.com/bcbio/giab/AJTrio/NA24385-callable.bed.gz | gunzip -c > {output}"""

rule truth_regions_na24385:
    input: "calls/na24385-callable.bed"
    output: "work/na24385-truth-regions.bed"
    shell:
        "cat {input} | "
        r"""grep -v ^M | grep -v _random | grep -v ^Un_ """
        "> {output}"

rule truth_regions:
    input: "calls/na12878-giab.bed", "work/na24385-truth-regions.bed"
    output: "na12878-na24385-somatic-truth-regions.bed"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} > {output}"
