## Load Libraries for analyses

## Setup directories and variables
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/ds_PTGIR_promoter")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/ds_PTGIR_promoter"
} else {
  baseDir <- "/Volumes/ody/consults/ds_PTGIR_promoter"
}
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results")


array.data = read.csv(file.path(dataDir,"57480002-6_DiffExp_Results_All_6_Samples_Simultaneously.csv" ),  header=T)
write(as.character(array.data$SYMBOL), file=file.path(dataDir, "background.gene.set.txt")) 

array.data.sig <- array.data[which(abs(array.data$logFC)>=1 & array.data$adj.P.Val<=0.1),]

enriched.TFs=c("Fev","Hoxa5","Mef2A","Runx1","Tal1","Gata1","Hnf1A","Evi1","Pbx1","Plag1","Gabpafe","T","Esrrb","Stat1","Nkx3-1","Egr1","Nfatc2","E2F1","Elk1","Tal1","Tcf3","Myf","Stat3","Plag1","Esr2","Lhx3","Foxa1","Pbx1","Spib","Pparg","Rxra","Hnf1B","Hoxa5","Nfatc2","Sp1","Nfya","Tal1","Gata1","Stat1","Nkx3-1","Foxd3","Srf","Pdx1","Rreb1")

write.table(array.data[(array.data$SYMBOL %in% enriched.TFs),], quote=F, sep="\t", col.names=T, row.names=F, file=file.path(resultsDir, "enriched.TFs.expression.levels.tab"))



