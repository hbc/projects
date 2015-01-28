baseDir <- "/Volumes/ody/projects/ab_rrbs_comparison/results"
chrMDir <- file.path(baseDir, "measuring_bisulfite_conversion", "chrM")
CpGIDir <-   file.path(baseDir, "measuring_bisulfite_conversion", "CpG_islands")

QCresults <- read.csv(file.path(baseDir, "QCsummary.csv"))
names(QCresults) <- sub("X._", "", names(QCresults))
row.names(QCresults) <- QCresults$sampleID

CpGInonCpG <- read.table(file.path(CpGIDir, "nonCpGmeth.txt" )) 
CpGInonCpG$V2 <- as.numeric(sub("%", "", CpGInonCpG$V2))
row.names(CpGInonCpG) <- CpGInonCpG$V1

CpGICpG <- read.table(file.path(CpGIDir, "CpGmeth.txt" )) 
CpGICpG$V2 <- as.numeric(sub("%", "", CpGICpG$V2))
row.names(CpGICpG) <- CpGICpG$V1

chrMnonCpG <- read.table(file.path(chrMDir, "nonCpGmeth.txt" ))
chrMnonCpG$V2 <- as.numeric(sub("%", "", chrMnonCpG$V2))
row.names(chrMnonCpG) <- chrMnonCpG$V1

chrMCpG <- read.table(file.path(chrMDir, "CpGmeth.txt" )) 
chrMCpG$V2 <- as.numeric(sub("%", "", chrMCpG$V2))
row.names(chrMCpG) <- chrMCpG$V1


methmeasures <- do.call(cbind, list(CpGICpG, CpGInonCpG, chrMCpG, chrMnonCpG))
methmeasures <- methmeasures[,grep("V2", names(methmeasures))]
names(methmeasures) <- c("CpGI_CpG_meth", "CpGI_nonCpG_meth", "chrM_CpG_meth", "chrM_nonCpG_meth")

results <- cbind(QCresults, methmeasures)
results$shortids <- unlist(lapply(strsplit(row.names(results),"_"), function(x) paste(x[1:3], collapse="_")))


write.csv(results, "QC_regionalmeth_results.csv" )

pdf(file.path(baseDir, "nonCpG_methylation_and_mapping_efficiency.pdf"), width=11, height=8.5)

p <- ggplot(results, aes(x=mapping_efficiency, y=CpGI_nonCpG_meth, color=methylated_CHG))+  geom_point(size=5)
p <- p+geom_text(data=subset(results, CpGI_nonCpG_meth > 2.7), 
            aes(mapping_efficiency, CpGI_nonCpG_meth,label=shortids),
            color="black", size=2, hjust=-0.2, angle=15,alpha=0.8, positions=position_jitter())
p <- p+ggtitle("CpG Island non-CpG methylation (unconverted cytosines) and mapping efficiency")+xlab("% reads mapped")+ylab("% non-CpG methylation in CpG Islands")
p+ scale_color_gradient(name="% non-CpG \nmethylation\n(global)")

p <- ggplot(results, aes(x=mapping_efficiency, y=chrM_nonCpG_meth, color=methylated_CHG))+  geom_point(size=5)
p <- p+geom_text(data=subset(results, chrM_nonCpG_meth > 5), aes(mapping_efficiency, chrM_nonCpG_meth,label=shortids),color="black", size=2,  hjust=-0.2, alpha=0.8, positions=position_jitter())
p <- p+ggtitle("Mitochondrial non-CpG methylation (unconverted cytosines) and mapping efficiency")+xlab("% reads mapped")+ylab("% non-CpG methylation on mitochondrial DNA")
p+ scale_color_gradient(name="% non-CpG \nmethylation\n(global)")
dev.off()
