## ----setup---------------------------------------------------------------
if(file.exists("/Users/johnhutchinson/projects")){
  baseDir  <- "/Users/johnhutchinson/projects/lk_FOY/avon450K"
} else if (file.exists("/n/home08/jhutchin/regal")){
baseDir <- "/n/regal/hsph_bioinfo/jhutchin/lk_FOY/avon450k"
}
dataDir <- file.path(baseDir, "data//ARIES_data_B1361")
resultsDir <- file.path(baseDir, "results")
metaDir <- file.path(baseDir, "meta")

## ----libraries-----------------------------------------------------------
library(ff)

## ----dataload, cache=TRUE------------------------------------------------
load(file.path(dataDir, "B1361.rawBetas.F7releasev2_apr2014.Rdata"))
betas.7 <- betas
colnames(betas.7) <-  paste(colnames(betas.7), "7", sep="_")
rm(betas)


betas.ff=ff(betas.7)

ffsave.image("test")
rm(betas.ff)
gc()
q()
