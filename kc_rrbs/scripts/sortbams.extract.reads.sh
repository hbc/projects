####################################################################################################################################
##VARIABLES
####################################################################################################################################

BAMDIR="/n/scratch00/hsph/projects/kc_rrbs/data/bams"
FASTQDIR="/n/scratch00/hsph/projects/kc_rrbs/data/fastq"
SAMPLES="S1_Mouse_ES_V16_p8_RRBS
S2_Mouse_iPS_TTF_derived_T20_p8_RRBS
S3_Mouse_iPS_TTF_derived_T21_p8_RRBS
S4_Mouse_iPS_TTF_derived_T23_p8_RRBS
S5_Mouse_iPS_VM_derived_C61_p8_RRBS
S6_Mouse_iPS_VM_derived_C63_p8_RRBS
S7_Mouse_iPS_VM_derived_C64_p8_RRBS
S8_Mouse_tailtip_fibroblast_RRBS
S9_Mouse_ventricular_myocyte_RRBS"

queue="hsph"

tmpDir="/n/scratch00/hsph/tmp"

####################################################################################################################################
##FUNCTIONS
####################################################################################################################################

wait.on.lsf() { 
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

####################################################################################################################################
##ANALYSES
####################################################################################################################################

##################################################################

##revert sam file
echo "reverting Sam"
projectID="revert_sam"	

for SAMPLE in $SAMPLES
do 
	bsub -P $projectID -q ${queue} -e ${BAMDIR}/err.revert.sam -o ${BAMDIR}/out.revert.sam  \
	java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar ~/HSPH/local/share/java/picard/RevertSam.jar \
	INPUT=${BAMDIR}/${SAMPLE}.bam \
	OUTPUT= ${BAMDIR}/${SAMPLE}.reverted.bam 
done
wait.on.lsf

##################################################################
		
##EXTRACT reads to fastq
echo "extracting fastqs"
projectID="extract_fastq"	

for SAMPLE in $SAMPLES
do 
	bsub -P $projectID -q ${queue} -e ${FASTQDIR}/err.extract.fastq -o ${FASTQDIR}/out.extract.fastq  \
	java -Xmx2g -Djava.io.tmpdir=${tmpDir} -jar ~/HSPH/local/share/java/picard/SamToFastq.jar \
	INPUT=${BAMDIR}/${SAMPLE}.reverted.bam \
	FASTQ= ${FASTQDIR}/${SAMPLE}.fastq \
	INCLUDE_NON_PF_READS=TRUE \
	INCLUDE_NON_PRIMARY_ALIGNMENTS=TRUE
done
wait.on.lsf

