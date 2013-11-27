## 2012-05-14 - first draft
###jnhutchinson
###bash script to run RNA-seqQC on a group of alignments

###hard coded

#resultsDir	---untrimmed
#				---sample
#				---RNAseqQC
#			---trimmed
#				---sample
#				---RNAseqQC			
				
#################################################################################
###VARIABLES
alignDir="/n/home08/jhutchin/scratch/projects/wp_rna_seq/results/tophataligns"
GTFfile="/n/home08/jhutchin/scratch/tmp/JH/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
bowtiegenome="/n/home08/jhutchin/scratch/tmp/JH/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome"
trim_size=11

if [ $trim_size -eq 0 ]
then
	trimvarDir="untrimmed"
else 
	trimvarDir="trimmed"
fi

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


wait.on.lsf() { ## wait on jobs{
	n=`bjobs -P $projectID | awk 'NR>1' | wc -l`
	while [ $n -ne 0 ]
	do
		n=`bjobs -P $projectID | awk 'NR>1' | wc -l`
	 	# number of running
		echo "Running: "`bjobs -P $projectID | grep RUN | wc -l`
	 	# number of pending
		echo "Pending: "`bjobs -P $projectID | grep PEND | wc -l`
		sleep 60
	done
}


#################################################################################
##PREP BAM FILES
#################################################################################
echo "redo headers on Tophat bam output"
 
projectID="headers"
for sample in $samples
do 
	bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
	java -Xmx2g -jar /n/HSPH/local/share/java/picard/AddOrReplaceReadGroups.jar \
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


#################################################################################
echo "co-ordinate sort"
projectID="coordinate_sort"
for sample in $samples
do
	bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
	java -Xmx2g -jar /n/HSPH/local/share/java/picard/SortSam.jar \
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

#################################################################################
echo "reorder reads to match contigs in the reference"
projectID="reorder_contigs"
for sample in $samples
do	
	bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
	java -Xmx2g -jar /n/HSPH/local/share/java/picard/ReorderSam.jar \
	INPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.bam \
	OUTPUT=${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.bam \
	REFERENCE=${bowtiegenome}.fa
done
##wait on jobs
wait.on.lsf

for sample in $samples
do
	rm ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.bam  
done

#################################################################################
echo "mark duplicates"
projectID="markdups"
for sample in $samples
do	
	bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
	java -Xmx2g -jar /n/HSPH/local/share/java/picard/MarkDuplicates.jar \
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

#################################################################################
echo "index"
projectID="index"
for sample in $samples
do	
	bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${sample}/err.${projectID} -o ${alignDir}/${trimvarDir}/${sample}/out.${projectID} \
	samtools index ${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam
done
##wait on jobs
wait.on.lsf

#################################################################################
##RNA-seqQC
projectID="RNAseqQC"
echo "make sample list"
for sample in $samples
do
	echo -e "$sample\t${alignDir}/${trimvarDir}/${sample}/accepted_hits.readgroups.sorted.reordered.dedup.bam\tNA" >>${alignDir}/${trimvarDir}/${projectID}/sample.list
done
wait

echo "run RNA-seqQC"
bsub -P $projectID -q normal_serial -e ${alignDir}/${trimvarDir}/${projectID}/err.${projectID} -o ${alignDir}/${trimvarDir}/${projectID}/out.${projectID} \
java -Xmx2g -jar ~/.local/bin/RNA-SeQC_v1.1.6.jar \
-o ${alignDir}/${trimvarDir}/${projectID}/ \
-r ${bowtiegenome}.fa \
-s ${alignDir}/${trimvarDir}/${projectID}/sample.list \
-t /n/home08/jhutchin/scratch/tmp/JH/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
-ttype 3 \
-gld





