# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Feb 2 15:04:20 EST 2018
# Derived from GEO2R

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(rio)

# load series and platform data from GEO

gset <- getGEO("GSE101584", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6885", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "XXXXXXX111XXXX00000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)

# > design
# G0 G1
# GSM2706387  0  1
# GSM2706388  0  1
# GSM2706389  0  1
# GSM2706394  1  0
# GSM2706395  1  0
# GSM2706396  1  0
# GSM2706397  1  0
# GSM2706398  1  0

## Here the G1 sampels are the N1 neutrophils () and the GO are the N2s.
## I set up the contrast so that positive fold changes are when a gene is higher in the N1s (i.e G1/G0 or N1/N2)



colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G0-G1, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex))

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
export(tT, file=file.path("data", "GSE101584_N2_over_N1.stats.csv"))


