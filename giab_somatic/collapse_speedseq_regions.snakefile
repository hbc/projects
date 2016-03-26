# Combine parallelization regions from speedseq to be in smaller total regions

rule all:
    input: "regions/ceph18-b37-include-2014-01-15-stats.txt"

rule region_stats:
    input: "regions/ceph18-b37-include-2014-01-15-collapse.bed"
    output: "regions/ceph18-b37-include-2014-01-15-stats.txt"
    shell:
        "wc -l {input} > {output} && "
        "awk '{{print $3 - $2}}' {input} | datamash --header-out min 1 q1 1 median 1 q3 1 max 1 >> {output}"

rule collapse_regions:
    input: "regions/ceph18.b37.include.2014-01-15.bed"
    output: "regions/ceph18-b37-include-2014-01-15-collapse.bed"
    #shell:
    #    "bedtools merge -d 1000 -i {input} > {output}"
    run:
        target_size = 2e7
        merge_size = 100000
        with open(input[0]) as in_handle:
            with open(output[0], "w") as out_handle:
                cur_chrom = None
                for line in in_handle:
                    chrom, start, end = line.strip().split()[:3]
                    start = int(start)
                    end = int(end)
                    if chrom != cur_chrom:
                        if cur_chrom:
                            out_handle.write("%s\t%s\t%s\n" % (cur_chrom, cur_start, cur_end))
                        cur_chrom = chrom
                        cur_start = start
                        cur_end = end
                    elif start - cur_end > merge_size or end - cur_start > target_size:
                        if cur_chrom:
                            out_handle.write("%s\t%s\t%s\n" % (cur_chrom, cur_start, cur_end))
                            cur_chrom = chrom
                            cur_start = start
                            cur_end = end
                    else:
                        cur_end = end
                if cur_chrom:
                    out_handle.write("%s\t%s\t%s\n" % (cur_chrom, cur_start, cur_end))


rule get_original:
    output: "regions/ceph18.b37.include.2014-01-15.bed"
    shell:
        "wget -O - https://github.com/hall-lab/speedseq/raw/master/annotations/ceph18.b37.include.2014-01-15.bed "
        "| sort -V -k1,1 -k2,2n | grep -v ^GL > {output}"
