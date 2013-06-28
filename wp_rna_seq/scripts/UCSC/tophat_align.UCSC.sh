## 2012-05-08 - first draft
## 2012-05-14 - got rid of repetitive parts of script
##              added RNA-seqQC prep
###jnhutchinson
###bash script to run tophat alignments for WP-rna-seq

###hard coded

#basedataDir    ---untrimmed
#               ---trimmed
#resultsDir ---untrimmed
#               ---sample(s)
#               ---RNAseqQC
#           ---trimmed
#               ---sample(s)
#               ---RNAseqQC         
                
#################################################################################
###VARIABLES
#################################################################################

basedataDir="/n/scratch00/hsph/projects/wp_rna_seq/data/fastq"
alignDir="/n/scratch00/hsph/projects/wp_rna_seq/results/tophataligns"
tmpDir="/n/scratch00/hsph/tmp"
queue="hsph"
cores=7
gap=150
GTFfile="/n/scratch00/hsph/biodata/genomes/Mmusculus/mm9/iGenomes/UCSC/mm9/Annotation/Genes/genes.gtf"
bowtiegenome="/n/scratch00/hsph/biodata/genomes/Mmusculus/mm9/iGenomes/UCSC/mm9/Sequence/BowtieIndex/genome"
rRNAfile="/n/home08/jhutchin/resources/mouse_rRNA/mouse_all_rRNA.fasta"
samplesuffix="fastq"
trim_size=$1 #first base to keep -1

##samples
samples="LIB003615_TRA00004588_CGATGT_L005
LIB003615_TRA00004674_TGACCA_L005
LIB003615_TRA00004675_ACAGTG_L005
LIB003615_TRA00004676_GCCAAT_L005
LIB003615_TRA00004677_CAGATC_L005
LIB003615_TRA00004678_TAGCTT_L005
LIB003615_TRA00004679_CTTGTA_L005
LIB003616_TRA00004589_CGATGT_L006
LIB003616_TRA00004680_TGACCA_L006
LIB003616_TRA00004681_ACAGTG_L006
LIB003616_TRA00004682_GCCAAT_L006
LIB003616_TRA00004683_CAGATC_L006
LIB003616_TRA00004684_TAGCTT_L006
LIB003616_TRA00004685_CTTGTA_L006"
/genomes/Mmusculus/mm9/iGenomes

##to trim or not to trim-t
if [ $trim_size -eq 0 ]
then
    trimvarDir="untrimmed"
else 
    trimvarDir="trimmed"
fi

####################################################################################################################################
##FUNCTIONS
####################################################################################################################################

wait.on.lsf() { ## wait on jobs{
    n=`bjobs -P $projectID | awk 'NR>1' | wc -l`
    while [ $n -ne 0 ]
    do
        n=`bjobs -P $projectID | awk 'NR>1' | wc -l`
        # number of running
        echo "Running: "`bjobs -P $projectID | grep RUN | wc -l`
        # number of pending
        echo "Pending: "`bjobs -P $projectID | grep PEND | wc -l`
        sleep 60 ##raising this to 120 will cause LSF to drop the script
    done
}

 
###############################################################################################################################
## 5'TRIM
###############################################################################################################################
  
echo "trim fastqs"
projectID="trims${trim_size}"
if [ $trim_size -gt 0 ]
then
  for sample in $samples
  do
	bsub -P $projectID -e ${basedataDir}/${trimvarDir}/err_trims -o ${basedataDir}/${trimvarDir}/out_trims -q ${queue} fastx_trimmer -f $trim_size -Q33 -i ${basedataDir}/untrimmed/${sample}_R1.${samplesuffix} -o ${basedataDir}/${trimvarDir}/${sample}_R1.${samplesuffix} 
    bsub -P $projectID -e ${basedataDir}/${trimvarDir}/err_trims -o ${basedataDir}/${trimvarDir}/out_trims -q ${queue} fastx_trimmer -f $trim_size -Q33 -i ${basedataDir}/untrimmed/${sample}_R2.${samplesuffix} -o ${basedataDir}/${trimvarDir}/${sample}_R2.${samplesuffix} 
  done
fi
# wait on jobs
wait.on.lsf

  
#################################################################################################################################
# ALIGN
#################################################################################################################################
  
echo "do alignments"
projectID="aligns${trim_size}"
for sample in $samples
do 
   if [ -e ${basedataDir}/${trimvarDir}/${sample}_R1.${samplesuffix} ] && [ -e ${basedataDir}/${trimvarDir}/${sample}_R2.${samplesuffix} ] # check files
   then
	  bsub -P $projectID -q ${queue} -e ${alignDir}/${trimvarDir}/${sample}/err -o ${alignDir}/${trimvarDir}/${sample}/out \
      tophat -p $cores -r $gap \
      -G $GTFfile \
      -o ${alignDir}/${trimvarDir}/${sample} \
      $bowtiegenome \
      ${basedataDir}/${trimvarDir}/${sample}_R1.${samplesuffix} ${basedataDir}/${trimvarDir}/${sample}_R2.${samplesuffix}
    fi  
done
#wait on jobs
wait.on.lsf
  
##################################################################################################################################
# PREP BAM FILES
##################################################################################################################################

echo "redo headers on Tophat bam output"
projectID="headers${trim_size}"
for sample in $samples
do 
  bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
  java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/AddOrReplaceReadGroups.jar \
  INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.bam \
  OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.bam \
  RGID=$sample \
  RGCN=CCCB \
  RGLB=HeartPU \
  RGPL=ILLUMINA \
  RGPU=NA \
  RGSM=accepted_hits.bam
done
##wait on jobs
wait.on.lsf

####################################################################################################################################
echo "co-ordinate sort"
projectID="coordinate_sort${trim_size}"
for sample in $samples
do
  bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
  java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/SortSam.jar \
  INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.bam \
  OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.bam \
  SO=coordinate
done
##wait on jobs
wait.on.lsf

for sample in $samples
do
  rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.bam  
done

####################################################################################################################################
echo "reorder reads to match contigs in the reference"
projectID="reorder_contigs${trim_size}"
for sample in $samples
do    
  bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
  java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/ReorderSam.jar \
  INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.bam \
  OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.bam \
  REFERENCE=${bowtiegenome}.fa
done
##wait on jobs
wait.on.lsf

for sample in $samples
do
  rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.bam  
#done

####################################################################################################################################
echo "mark duplicates"
projectID="markdups${trim_size}"
for sample in $samples
do    
  bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
  java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/MarkDuplicates.jar \
  INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.bam \
  OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam \
  METRICS_FILE=${alignDir}/${trimvarDir}/${sample}/dup_metrics.txt
done
##wait on jobs
wait.on.lsf

for sample in $samples
do
  rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.bam
done

####################################################################################################################################
 echo "index"
 projectID="index${trim_size}"
 for sample in $samples
 do    
   bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
   samtools index ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam
 done
 ##wait on jobs
 wait.on.lsf

####################################################################################################################################
##RNA-seqQC
####################################################################################################################################
#Jun4.2012 will fail if directory not present
 projectID="RNAseqQC${trim_size}"
 echo "Sample_ID	Bam_File	Notes" >${alignDir}/${trimvarDir}/${projectID}/sample.list
 echo "make sample list"
 for sample in $samples
 do
   echo -e "$sample\t${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam\tNA" >>${alignDir}/${trimvarDir}/${projectID}/sample.list
 done
 wait 
 
 echo "run RNA-seqQC"
 bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${projectID}/err.${projectID} -o ${alignDir}/${trimvarDir}/${projectID}/out.${projectID} \
 java -Xmx4g -Djava.io.tmpdir=${tmpDir} -jar ~/.local/bin/RNA-SeQC_v1.1.7.jar \
 -o ${alignDir}/${trimvarDir}/${projectID}/ \
 -r ${bowtiegenome}.fa \
 -s ${alignDir}/${trimvarDir}/${projectID}/sample.list \
 -t $GTFfile \
 -BWArRNA $rRNAfile \
 -rRNAdSampleTarget 100000

####################################################################################################################################
##HTseq-prep
####################################################################################################################################
 projectID="remove_chrM_hits${trim_size}"
 echo "removing chrM hits"
 for sample in $samples
 do
     bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
     samtools view -b \
     ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam \
     chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \
     -o ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.bam
 done
 #wait on jobs
 wait.on.lsf

####################################################################################################################################
 echo "index2"
 projectID="index2${trim_size}"
 for sample in $samples
 do    
   bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
   samtools index ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.bam
 done
 ##wait on jobs
 wait.on.lsf

####################################################################################################################################
 projectID="sort_bam_by_read_name${trim_size}"
 echo "sorting bam file by read name"
 for sample in $samples
 do
     bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
     java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/SortSam.jar \
     INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.bam \
     OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.bam \
     SORT_ORDER=queryname
 done
 #wait on jobs
 wait.on.lsf
 
 #cleanup
 for sample in $samples
 do
 	rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.bam
   rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.bam.bai
 done

####################################################################################################################################
 projectID="convert_bam_to_sam${trim_size}"
 echo "converting bam file to sam file"
 for sample in $samples
 do
     bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
     java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar /n/HSPH/local/share/java/picard/SamFormatConverter.jar \
     INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.bam \
 	OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.sam
 done
 #wait on jobs
 wait.on.lsf
 
 #cleanup
 for sample in $samples
 do
 	rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.bam
 done

####################################################################################################################################
##HTseq-count
####################################################################################################################################
#need to activate Python virtual environment 'env2' - use source instead of bash when running file
 projectID="htseq-count${trim_size}"
 echo "htsseq-count"
 source /n/home08/jhutchin/.virtualenvs/env2/bin/activate
 for sample in $samples
 do
     bsub -P $projectID -q $queue -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/HTseq-counts.tab \
 	htseq-count -q --stranded=no \
 	${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.sam \
 	$GTFfile
 done 
 #wait on jobs
 wait.on.lsf

 #cleanup
 for sample in $samples
 do
 	rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.nochrM.readnamesorted.sam
 done
 
for sample in $samples
do
	sed -n '/Rik/,$p' ${alignDir}/${trimvarDir}/${sample}/HTseq-counts.tab |sed '/^$/q' >${alignDir}/${trimvarDir}/${sample}/HTseq-counts2.tab
	mv ${alignDir}/${trimvarDir}/${sample}/HTseq-counts2.tab ${alignDir}/${trimvarDir}/${sample}/HTseq-counts.tab
done

