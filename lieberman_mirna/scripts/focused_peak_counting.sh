for file in *filtered.bam
    do
        overlap_file=`basename $file .bam`_overlap.bam
        bedtools intersect -abam $file -b filtered_exon.gtf -u > $overlap_file
        tmp_file=`basename $overlap_file .bam`.tmp
        bedtools intersect -abam $overlap_file -b filtered_exon.gtf -wb -bed | cut -f1,2,3,4,21 | cut -f1 -d';' | sort | uniq> $tmp_file
        bed_file=`basename $overlap_file .bam`.bed
        python ~/cache/projects/lieberman_mirna/scripts/get_sequence.py $tmp_file $overlap_file > $bed_file
    done

#bedtools intersect -abam control_v2_v3.sorted.mapped.filtered.bam -b filtered.gtf -u > control_overlap.bam
#1249bedtools intersect -abam control_overlap.bam -b filtered.gtf -wb -bed | cut -f1,2,3,4,21 | cut -f1 -d';' | sort | uniq> control_overlap.tmp
