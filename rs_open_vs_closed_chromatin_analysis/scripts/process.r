##LIBRARIES
library(pheatmap)
library(RColorBrewer)

##FUNCTIONS
intersectSeveral <- function(...) { Reduce(intersect, list(...)) } 
setdiffSeveral <- function(...) { Reduce(setdiff)}
get.objectname <- function(x) deparse(substitute(x))


##VARIABLE
setwd("/n/home08/jhutchin/HSPH/projects/JH/rs_open_closed_chromatin/")
baseDir <- getwd()
dataDir <- paste(baseDir, "/data/wg/", sep="") 
resultsDir <- paste(baseDir, "/results/wg/", sep="")
suffix <- ".hg18.out"
samples <- c("adenoma1","adenoma2", "adenoma3", "adenoma4","adenoma5","carcinoma1","carcinoma2", "carcinoma3","carcinoma4", "carcinoma5")

##DATALOAD
files <- list.files(pattern=".out", path=dataDir)
for (file in files) {
  temp.data <- read.delim(paste(dataDir, file, sep=""))
  names(temp.data) <- c("q.int", "r.int", "q.len", "dist", "num", "denom", "p.value")
  assign(sub(suffix, "", file), temp.data)
  }
rm(file, temp.data)

##pvalue distribution, compared to uniform distribution?
pdf(paste(resultsDir, "pvalues.dist.pdf", sep=""))
  plot(0,0 , ylim=c(0,2), xlim=c(0,1), col="white", main="pvalue.distributions")
 for (object in ls()[grep("adenoma..adenoma", ls())]) {
   lines(density(get(object)$p.value), col="blue")
   }
for (object in ls()[grep("carcinoma..carcinoma", ls())]) {
   lines(density(get(object)$p.value), col="black")
   }
for (object in ls()[grep("adenoma..carcinoma", ls())]) {
   lines(density(get(object)$p.value), col="red")
   }
for (object in ls()[grep("carcinoma..adenoma", ls())]) {
   lines(density(get(object)$p.value), col="pink")
   }
lines(density(runif(1000, 0,1)), col="grey")
legend("topright", col=c("blue", "black", "red", "pink", "grey"), pch=15, legend=c("intra-adenoma", "intra-carcinoma", "adenoma-carcinoma", "carcinoma-adenoma", "uniform.dist"))
dev.off()


#qqplots against normal distribution
objects <- sort(unique(c(paste(combn(samples, 2)[1,], combn(samples, 2)[2,], sep="."), paste(combn(samples, 2)[2,], combn(samples, 2)[1,], sep="."))))
plotme <- data.frame() 
for(objectname in objects){
   data.temp <- get(objectname)
   data.temp$id <- objectname
  plotme <- rbind(plotme, data.temp)
   }

pdf(paste(resultsDir, "qqplots.pdf", sep=""))
for (x in samples){
  p <- qplot( sample=-log(p.value), data=plotme[grep(paste("^", x,sep=""), plotme$id),], distribution = qexp) + 
  stat_abline(intercept=0,slope=1, col="red") + 
  facet_wrap(~id, ncol=3) + opts(aspect.ratio = 1) + 
  scale_y_continuous("Observed -log10(pvalue)") +
  scale_x_continuous("Theoretical -log10(pvalue)")
  print(p)
}
dev.off()


###similarity score heatmap
p.value.cutoff <- 0.01
sim.matrix <- matrix(nrow=length(samples),ncol=length(samples))
dimnames(sim.matrix)[[1]]=samples
dimnames(sim.matrix)[[2]]=samples
for (sample1 in samples){
  for (sample2 in samples) {
    temp.data <- get(paste(sample1, sample2,sep="."))
    sim <- length(which(temp.data$p.value<p.value.cutoff))
    sim.matrix[sample1, sample2] <- sim
  }
}
diag(sim.matrix) <- NA
pheatmap(sim.matrix, col=brewer.pal(9,"Greys"), cluster_cols=F, cluster_rows=F, filename=paste(resultsDir, "heatmap.", p.value.cutoff, ".pdf", sep=""))
rm(temp.data)

##bring it all comparison data together into single file
objects <- sort(unique(c(paste(combn(samples, 2)[1,], combn(samples, 2)[2,], sep="."), paste(combn(samples, 2)[2,], combn(samples, 2)[1,], sep="."))))
for (sample in samples) {
   sample.objects <- objects[grep(paste("^", sample, sep=""),  objects)]
   cols=c("q.int", "r.int", "dist", "p.value")
   output <- get(sample.objects[1])[,cols]
   for (i in 2:length(sample.objects)) {
        output <- merge(output, get(sample.objects[i])[,cols], by="q.int", all=T)
   }
  names(output) <- c(paste(sample, "peak", sep="."), as.vector(t(outer(sub(paste(sample, ".", sep=""), "", sample.objects), c("peak","dist", "p.value"), paste, sep="."))))
  assign(paste(sample, "peak.dist.pvalue", sep="."), output)
}
for (objectname in ls()[grep("dist.pvalue", ls())]) {
  data.temp <- get(objectname)
  write.table(data.temp, file=paste(resultsDir, objectname, ".xls", sep=""), sep="\t", col.names=T, row.names=F)
}




##########################################################################################################################################################
##Potential Analyses
##########################################################################################################################################################

objects.pvalues <- ls()[grep("peak.dist.pvalue", ls())]
for (object.pvalue in objects.pvalues){
        temp.data <- get(object.pvalue)
         names.temp <- names(temp.data)
         adenoma.pvalues <- temp.data[,c(1, grep("adenoma.+p.value", names.temp))]
         carcinoma.pvalues <-temp.data[,c(1,grep("carcinoma.+p.value", names.temp))]

         ad.temp <- adenoma.pvalues[,2:ncol(adenoma.pvalues)]
         carc.temp <- carcinoma.pvalues[,2:ncol(carcinoma.pvalues)]
         data.temp <- cbind(ad.temp, carc.temp)
         wilcox.pvalues <- apply(data.temp, 1, function(n) wilcox.test(as.numeric(n[1:ncol(ad.temp)]), as.numeric(n[(ncol(ad.temp)+1):ncol(data.temp)]))$p.value)
         ks.pvalues <- apply(data.temp, 1, function(n) ks.test(as.numeric(n[1:ncol(ad.temp)]), as.numeric(n[(ncol(ad.temp)+1):ncol(data.temp)]))$p.value)
         
         adenoma.pvalues.melt <- melt(adenoma.pvalues)
         adenoma.pvalues.melt$sampletype <- "adenoma"
         carcinoma.pvalues.melt <- melt(carcinoma.pvalues)
         carcinoma.pvalues.melt$sampletype <- "carcinoma"
         plotme.data <- rbind(adenoma.pvalues.melt, carcinoma.pvalues.melt)
         plotme.data.subset <- plotme.data[which(plotme.data[,1] %in% temp.data[which(ks.pvalues<0.05),1]),]
        plotme.data.subset$id <- plotme.data.subset[,1]
        pdf(paste(object.pvalue), "pdf", sep=".") 
        qplot(log(value), data=plotme.data.subset,  geom="density", colour=sampletype) + facet_wrap(~id, scales="free_y")  
        dev.off()
}
        
                                                       
         
         

           
            