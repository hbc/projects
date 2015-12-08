

## ----libraries-----------------------------------------------------------
library(knitr)
library(ggplot2) # pretty plots
library(minfi)
library(xlsx) # load Excel files
library(wateRmelon) #pfilter and BMIQ funcitons
library(IlluminaHumanMethylation450k.db)
library(devtools)
install_git("git://github.com/hbc/CHBUtils.git") # misc personal utilities
library(CHBUtils)
library(methylumi)
library(beanplot)


## ----variables-----------------------------------------------------------
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_450k/")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_450k/"
  } else {
  baseDir <- "/Volumes/ody/consults/lk_450k/"
  }
dataDir <- file.path(baseDir, "data")
metaDir <- file.path(baseDir, "meta")
resultsDir<-file.path(baseDir, "results")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## ----functions-----------------------------------------------------------
densityBeanPlot.JH <- function (dat, sampGroups = NULL, sampNames = NULL, main = NULL, pal = branrewer.pal(12, "Set3"), numPositions = 10000, label.size=0.8) 
{
    if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
        b <- getBeta(dat)
    }
    else if (is(dat, "matrix")) {
        b <- dat
    }
    else {
        stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet'  or matrix.")
    }
    n <- ncol(b)
    if (!is.null(sampNames)) 
        colnames(b) <- sampNames
    if (is.null(main)) 
        main <- "Beta"
    if (is.null(sampGroups)) 
        sampGroups <- rep(1, n)
    sampGroups <- as.factor(sampGroups)
    col <- lapply(sampGroups, function(x) rep(pal[x], 4))
    if (is.null(numPositions)) 
        idx <- 1:dim(dat)[1]
    else idx <- sample(nrow(b), numPositions)
    x <- melt(b[idx, ], varnames = c("cpg", "sample"))
    o <- order(colnames(b))
    beanplot(value ~ sample, horizontal = TRUE, what = c(0, 1, 
        1, 0), log = "", las = 1, ylim = c(0, 1), xlab = "Beta", 
        main = main, col = col[o], data = x, cex.axis=label.size, cex.lab = 0.9, beanlinewd = 1, 
        border = NA)
    abline(h = 1:(n + 1) - 0.5, lty = 3, col = "grey70")
}


## ----metadataimport, results='hide'--------------------------------------
samplesheet <- read.delim(file.path(metaDir, "Samples_Table.txt" ), sep="\t")
metadata <- as.matrix(read.xlsx(file.path(metaDir, "human-methylation450-96samples-sample-sheet_Kobzik.xlsx"), sheetIndex = 2, header=F))

# fix the name of the  replicate a5 sample in the metadata so you can match to samplesheet
metadata[which(metadata[,"X1"]=="b12" & metadata[,"X2"]=="a5"), 2] <- "a5_2"
metadata <- as.data.frame(metadata)
# capitalize the array positionson the metadata so you can match by it as well
metadata$X1 <- toupper(metadata$X1)

# sanity check to make sure you have metadata for all samples
all(samplesheet$Sample.ID %in% metadata$X2)
all(metadata$X2 %in% samplesheet$Sample.ID)

# merge all sample metadata
merged_metadata <- merge(metadata, samplesheet, by.x=c("X1", "X2"), by.y=c("Sample_Well", "Sample.ID"))

# setup for minfi & methylumi import functions
merged_metadata$barcode<-paste(merged_metadata$Sentrix.Barcode, merged_metadata$Sample.Section, sep="_")
merged_metadata$Basename <- file.path(dataDir, merged_metadata$barcode)

names(merged_metadata) <- sub("Sentrix.Barcode", "Slide", names(merged_metadata))
names(merged_metadata) <- sub("Sample.Section", "Array", names(merged_metadata))
names(merged_metadata) <- sub("X1", "Position", names(merged_metadata))
names(merged_metadata) <- sub("X2", "SampleID", names(merged_metadata))
names(merged_metadata) <- sub("X3", "group", names(merged_metadata))
# dump extra columns
merged_metadata <- merged_metadata[,!grepl("Detected|Signal|Pool_ID|Sample_Plate|Index|Sample.Group", names(merged_metadata))]
# add this so that methylumi can get the sampleNames on iDAT import
row.names(merged_metadata) <- merged_metadata$barcode


## ----printmetadata, results='asis'---------------------------------------
kable(metadata)


## ----dataimport, cache=TRUE----------------------------------------------
#minfi 
## method, use etended option to get an RGChannelSEtextended, necessary for pfiltering in next stepsz
if(file.exists(file.path(resultsDir, "RDATA.RGset"))){
  load(file.path(resultsDir, "RDATA.RGset"))
  } else {
    RGset <- read.450k.exp(base=dataDir, targets=merged_metadata, extended = TRUE)
    save(RGset, file=file.path(resultsDir, "RDATA.RGset"))
    }

#methylumi
barcodes <- unique(sub("_Grn.idat|_Red.idat", "", list.files(path = dataDir, pattern = "idat$")))

if(file.exists(file.path(resultsDir, "RDATA.mldat"))){
  load(file.path(resultsDir, "RDATA.mldat"))
  } else {  
    mldat <- methylumIDAT(barcodes=barcodes, idatPath=dataDir, parallel=TRUE, oob=TRUE)
    save(mldat, file=file.path(resultsDir, "RDATA.mldat"))
    }
# merge in sample metadata
pData(mldat) <- merge(pData(mldat), merged_metadata)
# add in Ilumina annotations, here we grab them from the array facility's GenomeStudio output
annots <- read.delim(file.path(dataDir, "Group Methylation Profile.txt"), sep="\t")
annots <- annots[,!grepl("Group", names(annots))]
fData(mldat) <- merge(fData(mldat), annots, by.x="Probe_ID", "TargetID")
q()
