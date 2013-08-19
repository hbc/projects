for file in `find ~/cache/projects/wen_hiv/results/cutadapt -name miseq*_trimmed.fastq`
do
bsub -q hsph bowtie2 -S ~/cache/projects/wen_hiv/results/align/`basename $file .fastq`.sam ~/hsph/biodata/genomes/Hsapiens/GRCh37/bowtie2/GRCh37 $file
done
