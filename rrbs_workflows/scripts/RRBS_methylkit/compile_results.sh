filename=$1
sampleID=${filename%"quant_meth_methlykit.md"}
newfile=${sampleID}temp.txt


sed 's/$/  /g' ${sampleID}.fastq_trimming_report.txt >${sampleID}.trimming.report.tmp
sed 's/$/  /g' ${sampleID}.trimmed.fq_Bismark_mapping_report.txt >${sampleID}.mapping.report.tmp



echo -e "\n#RRBS QC Report for sample" "${sampleID}  \n" >>${sampleID}_QC_report.md
echo "---" >>${sampleID}_QC_report.md


echo -e "\n##Quality and Adapter Trimming Report  \n" >>${sampleID}_QC_report.md
cat ${sampleID}.trimming.report.tmp  >>${sampleID}_QC_report.md
echo -e "\n[FASTQC before trimming](./fastqc/${sampleID}.fastq/pretrim/${sampleID}_fastqc/fastqc_report.html)\n" >>${sampleID}_QC_report.md
echo -e "\n[FASTQC after  trimming](./fastqc/${sampleID}.fastq/posttrim/${sampleID}.trimmed.fq_fastqc/fastqc_report.html)\n" >>${sampleID}_QC_report.md

echo -e "\n---\n" >>${sampleID}_QC_report.md

echo -e "\n##Bismark Alignment Report  \n" >>${sampleID}_QC_report.md
cat ${sampleID}.mapping.report.tmp >>${sampleID}_QC_report.md
echo -e "\n---\n" >>${sampleID}_QC_report.md 

cat $filename >>${sampleID}_QC_report.md

pandoc -f markdown -t html -c http://dl.dropbox.com/u/4253254/CSS/GitHub.css  ${sampleID}_QC_report.md -o ${sampleID}_QC_report.html
