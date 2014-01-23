##RUN script from base directory, 
# with subdirectories named "results" and "data" (contents are descrived by title)


################################
##LIBRARIES
library(affy)
library(simpleaffy)
library(arrayQualityMetrics)
library(limma)
library(mouse430a2.db)
library(pheatmap)
library(RColorBrewer)
library(sva)



####################################################################################################################
##LOAD DATA
  #load in phenotype data from covars.desc in "data" subdirectory
  #this file contains the names and descriptions of CEL files contained in same directory 
  #this info will then be used to load in the CEL files
  #stores all this data in a single microarray "Affybatch" object
  
mic.data <- read.affy('covars.desc', path='./data', verbose=T);
  #dump unwanted first column that read.affy() is labelling "sample"
pData(mic.data)=pData(mic.data)[,-1]


####################################################################################################################
##QA/QC

# use raw data to create first QA/QC report
arrayQualityMetrics(expressionset=mic.data, outdir='./results/report_raw', force=TRUE, do.logtransform=TRUE)

#---------------------------------------------------------------------------------------------------------------------------------------
# Normalize 
mic.edata <- call.exprs(mic.data, "rma"); #mic.edata now contains the normalized microarray data
# Create QA/QC report for normalized data
arrayQualityMetrics(expressionset=mic.edata, outdir='./results/report_rma', force=TRUE, do.logtransform=FALSE)

#---------------------------------------------------------------------------------------------------------------------------------------
## Primary Clustering # are the samples clustering by sample type?
allData <- exprs(mic.edata) #extract expression values from normalized "Affybatch" object into a separate matrix 
colnames(allData) <- pData(mic.edata)$Sample
# Primary PCA
myPca <- prcomp(t(allData))
  # SD of components
	pdf("./results/PCA.nobatchcorrection.pdf")
	plot(myPca$sdev, xlab="PCA", ylab="sddev")
  # Plot samples in all combinations in the 5 first primary components
	tmpPCAData <- as.data.frame(myPca$x[,1:5])
	colors <- c('darkgrey','darkgrey', 'blue', 'cyan', 'red','pink')
	plot(tmpPCAData, col=colors, pch=row.names(myPca$x))
	dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------
# Batch correction
  # create model with Batch as factor variable
mod <- model.matrix(~as.factor(Batch) , data=pData(mic.edata))
mod <-  model.matrix(~as.factor(Condition) , data=pData(mic.edata))
batch <- pData(mic.edata)$Batch
  # Modify expression matrix accounting for Batch effect using ComBat method
edata <- exprs(mic.edata)

combat_edata <- ComBat(dat=edata,
                       batch=batch,
                       mod=mod,
                       numCovs=NULL,
                       par.prior=TRUE, 
                       prior.plots=TRUE)
mic.batch.edata <- mic.edata
exprs(mic.batch.edata) <- combat_edata
	# Post batch correction QA/QC report-----------------------------------------------------------------------------------------------
	arrayQualityMetrics(expressionset=mic.batch.edata, outdir='./results/batchcorrected_report_rma', force=TRUE, do.logtransform=FALSE)
  	# Post batch correction clustering-------------------------------------------------------------------------------------------------
	allData <- exprs(mic.batch.edata)
	colnames(allData) <- pData(mic.batch.edata)$Sample
  	# Post batch correction PCA--------------------------------------------------------------------------------------------------------
	myPca <- prcomp(t(allData))
    	# SD of components
		pdf("./results/PCA.batchcorrection.pdf")
		plot(myPca$sdev, xlab="PCA", ylab="sddev")
    	# Plot samples in all combinations in the 5 first components
		tmpPCAData <- as.data.frame(myPca$x[,1:5])
		colors <- c('darkgrey','darkgrey', 'blue', 'cyan', 'red','pink')
		plot(tmpPCAData, col=colors, pch=row.names(myPca$x))
		dev.off()


##############################################################################################################################
## ANALYSIS
# limma

# before doing analysis, output all the normalized and batch corrected expression values to table
# add gene name, symbol and EntrezID information before outputting data
eset.all.annotated <- as.data.frame(exprs(mic.batch.edata))
eset.all.annotated$Symbol <- unlist(mget(row.names(eset.all.annotated), mouse430a2SYMBOL, ifnotfound=NA))
eset.all.annotated$GeneName <- unlist(mget(row.names(eset.all.annotated), mouse430a2GENENAME, ifnotfound=NA))
eset.all.annotated$EntrezID <- unlist(mget(row.names(eset.all.annotated), mouse430a2ENTREZID, ifnotfound=NA))
write.table(eset.all.annotated, file='./results/all.exprs.xls', row.names=T, col.names=NA, sep='\t')

# Create appropriate design matrix (makes a matrix with arrays as rows, sample groups as columns)
# a one or a zero indicate respectively, that a sample either belongs or does not belong to the sample group, 
# in this case it looks like thi, where 1-6 are the 6 microarrays named in the covars.desc file in the same order as in the file
#   KO OV WT
# 1  1  0  0
# 2  1  0  0
# 3  0  0  1
# 4  0  0  1
# 5  0  1  0
# 6  0  1  0

design <- model.matrix(~ -1+factor(pData(mic.batch.edata)$Condition))
# make sure the headings match
colnames(design) <- c("KO", "OV", "WT")

# Fit a linear model for each gene based on the given series of arrays
eset <- exprs(mic.batch.edata)
fit <- lmFit(eset, design) 

# Create contrast matrix to perform specified pairwise comparisons
# the following command will form this table, where columns are contrasts/comparisons and rows are sample groups. 
# a zero denotes that the sample group is not involved in the contrast, a 1 denotes that it has higher expression in the contrast and a -1 denotes lower expression in the contrast
#  Contrasts
#Levels OV-WT WT-KO OV-KO
#    KO     0    -1    -1
#    OV     1     0     1
#    WT    -1     1     0
contrast.matrix <- makeContrasts("OV-WT", "WT-KO","OV-KO", levels=c("KO", "OV", "WT"))

# Computes estimated coefficients and standard errors for contrasts
fit2 <- contrasts.fit(fit, contrast.matrix) 

# Computes moderated t-statistics and log-odds of differential expression 
# by empirical Bayes shrinkage of the standard errors towards a common value.
fit2 <- eBayes(fit2) 

# to make sure you are getting a consistent increase or decrease in expression from KO->WT->OV, use decideTests function
# this function will assign a value of -1,0 or 1 depending on how well the data at a probeset fits with the contrast, a 0 indicates no significant difference between the two groups, a 1 indicates a significant difference *in the direction specified in the contrast matrix* and a -1 indicates a significant difference in opposite direction as specified by the contrast matrix
# For example: given our contrast matrix, in the OV-WT contrast, a 0 result from decideTests indicates no difference between OV and WT samples, a 1 indicates that the OV sample has significantly higher expression than WT and a -1 indicates that WT samples have significantly higher expression than OV samples   
# we want probesets that have consistentely increasing or decreasing significant differences between KO->WT samples and WT->OV samples (which then transitively implies consistent significant differences from KO->OV samples) 
# So, given that a 1 indicates a significant difference and the +/- value indicates the direciton, we want probesets that give results of 1,1,1 or -1,-1,-1 for the OV-WT, WT-KO and OV-KO contrasts
         
# Obtain results with significant differences in all contrasts 
# we say here that a significant result has a multiple test adjusted pvalue (using false discovery rate or "fdr") of less than 0.05
results=decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)
results.up <- row.names(results)[which(results[,1]==1 & results[,2]==1 & results[,3]==1)]
results.down <- row.names(results)[which(results[,1]==-1 & results[,2]==-1 & results[,3]==-1)]

# get the actual significance values of the contrasts for all probesets
x1 <- topTable(fit2, coef=c(1,2,3), adjust="fdr", sort.by="B", p.value=1,   number=nrow(fit2$genes), genelist=fit2$genes)
# add gene name and symbol annotations (if found)
x1$Symbol=unlist(mget(x1$ID, mouse430a2SYMBOL, ifnotfound=NA))
x1$GeneName=unlist(mget(x1$ID, mouse430a2GENENAME, ifnotfound=NA))
x1$EntrezID=unlist(mget(x1$ID, mouse430a2ENTREZID, ifnotfound=NA))
# get sigfnificance values for best results, ie. those with  significant consistent differences in all contrasts 
x1.sig <- x1[c(results.up,results.down),]
x1.sig <- x1.sig[order(x1.sig$adj.P.Val),]

# output results to files
# all F-statistic values
write.table(x1, file="./results/all.Fstats.xls", sep="\t", row.names=F, col.names=T)
# Fstatistic and expression values for probesets with consistent direction significant results across all contrasts 
write.table(x1.sig, file='./results/KO-WT-OV-3way.sig.Fstats.xls', row.names=F, sep='\t')
write.table(eset.all.annotated[x1.sig$ID,], file='./results/KO-WT-OV-3way.sig.exprs.xls', row.names=T, col.names=NA, sep="\t")


##############################################################################################################################
##VISUALIZATION
# Heatmaps

# generate matrices of normalized expression values for the top differentially expressed  hits
m1 = exprs(mic.edata[x1.sig[, "ID"], ]) 
colnames(m1) = pData(mic.data)$Sample
## add in gene symbols for each Affy probe, append to rowname with "paste", so you will be able to see it in the plot
row.names(m1)=paste(row.names(m1)," (", x1.sig$Symbol, ")", sep="")
# write out heatmaps with Affy probeset IDs
library(gplots)
pheatmap(m1, main=do.call(paste,as.list(colnames(fit2$coefficients))), color=rev(redgreen(40)), cluster_cols=F, fontsize_row=8, cellwidth=24,cellheight=12, lwd=0.1,filename="./results/KO-WT-OV-3way.sig.heatmap.pdf")
