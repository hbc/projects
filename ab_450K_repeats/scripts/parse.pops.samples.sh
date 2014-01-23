POPS="AFR
AMR
ASN
EUR"

SUBPOPS="ASW
CEU
CHB
CHS
CLM
FIN
GBR
IBS
JPT
LWK
MXL
PUR
TSI
YRI"

SNPs=$(cut -f6 Illumina.450K.intersect.1000G_omni2.5.hg19.JH2.xls)

for POP in $POPS
	do
		samples=$(grep $POP phase1_integrated_calls.20101123.ALL.panel.panel | cut -f1 | tr '\n' ',' |sed 's/,$//g')
		vcftools --vcf 1000G_omni2.5.hg19.vcf --keep ${samples} --snps ${SNPs} --frq --out Illumina.450K.intersect.1000G_omni2.5.hg19.${POP}
		echo $POP
done


#for SUBPOP in $SUBPOPS
#	do
#		samples=$(grep $SUBPOP phase1_integrated_calls.20101123.ALL.panel.panel | cut -f1 | tr '\n' ',' |gsed 's/,$//g')
#		vcf-subset -c ${samples} 1000G_omni2.5.hg19.vcf > ${SUBPOP}.1000G_omni2.5.hg19.vcf 
#		echo $SUBPOP
#done
#			
			
			
			
#			vcf-subset -c NA0001,NA0002 LIG1.vcf > pop.LIG1.vcf




#vcftools --vcf 1000G_omni2.5.hg19.vcf --keep ${samples} --snps ${SNPs} --frq --out Illumina.450K.intersect.1000G_omni2.5.hg19.${POP}