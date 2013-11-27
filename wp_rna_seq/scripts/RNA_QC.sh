#################################################################################
###VARIABLES
alignDir="/n/home08/jhutchin/scratch/projects/wp_rna_seq/results/tophataligns"
queue="normal_serial"
GTFfile="/n/home08/jhutchin/scratch/tmp/JH/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
bowtiegenome="/n/home08/jhutchin/scratch/tmp/JH/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome"
trim_size=11 #first base to keep

##sample
samples="LIB003615_TRA00004674_TGACCA_L005
LIB003615_TRA00004588_CGATGT_L005
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

##to trim or not to trim
if [ $trim_size -eq 0 ]
then
	trimvarDir="untrimmed"
else 
	trimvarDir="trimmed"
fi


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
bsub -P $projectID -q ${queue} -e ${alignDir}/${trimvarDir}/${projectID}/err.${projectID} -o ${alignDir}/${trimvarDir}/${projectID}/out.${projectID} \
java -Xmx2g -jar ~/.local/bin/RNA-SeQC_v1.1.6.jar \
-o ${alignDir}/${trimvarDir}/${projectID}/ \
-r ${bowtiegenome}.fa \
-s ${alignDir}/${trimvarDir}/${projectID}/sample.list \
-t ${GTFfile} \
-ttype 3 \
-gld
