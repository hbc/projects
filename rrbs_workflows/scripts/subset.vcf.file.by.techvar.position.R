## load in Alan's technical variation results
dataDir="~/projects/rrbs_workflows/sample_data/AB_data/BisSNP/concordance_checks/1_no_change"
vcffiles=file.path(dataDir, c("B1850_ATCACG_L002_R1.trimmed.fq_bismark.RG.RO.recal1.bam.rawcpg.vcf","C1850_CAGATC_L003_R1.trimmed.fq_bismark.RG.RO.recal1.bam.rawcpg.vcf"))

techvar=read.table(file.path(dataDir, "concordance_B1850_C1850_pipeline_original.tab"), header=T)

techvar$end=techvar$start+1

for (n in 0:2){
  subset.pos.file=tempfile()
  subset=techvar[which(techvar$techrep==n),c("chr", "start", "end")]
  names(subset)=c("chr", "start", "end")
  write.table(subset, file=subset.pos.file, col.names=T, row.names=F, quote=F, sep="\t")
  
  for (vcffile in vcffiles){
    vcf.subset.output = paste(sub(".vcf$", "", vcffile), n, "techvar", "vcf", sep=".")
    system(paste("vcftools --vcf", vcffile,"--bed", subset.pos.file, "--recode-to-stream --recode-INFO-all >", vcf.subset.output  ))
  }
  unlink(subset.pos.file)          
}