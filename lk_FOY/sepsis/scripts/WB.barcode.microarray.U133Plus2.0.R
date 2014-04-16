## ----barcode_general_libraries-------------------------------------------
source("http://bioconductor.org/biocLite.R") # BioConductor script necessary for installing new BioC libraries with biocLite()
library(xtable) # table generation for reports
library(plyr) # library for iteratively working with data structures
library(ggplot2) # library for plotting 
library(RColorBrewer) # library for generating color palettes
library(googleVis) # library for presenting tables


## ----barcode_general_directories-----------------------------------------
if (file.exists("/n/hsphS10/hsphfs1/chb/projects/lk_FOY/sepsis")) {
  baseDir <- "/n/hsphS10/hsphfs1/chb/projects/lk_FOY/sepsis"
  }  else if (file.exists("/Volumes/home08/jhutchin/consults/lk_FOY/sepsis")) {
    baseDir <- "/Volumes/home08/jhutchin/consults/lk_FOY/sepsis"
    } else {
      baseDir <- "/Volumes/ody/consults/lk_FOY/sepsis"
      }
dataDir <- file.path(baseDir, "data", "WB")
resultsDir <- file.path(baseDir, "results", "WB", "U133Plus2")
metaDir <- file.path(baseDir, "meta", "WB")
# colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


## ----barcode_load_affybatch, echo=FALSE, eval=FALSE----------------------
load(file.path(resultsDir, "RDATA.raw_and_normalized_microarray.data_U133Plus2.0"))


## ----barcode_libraries---------------------------------------------------
# frozen RMA and barcoding library
library(frma) 
# previously analyzed dataset from same array
library(hgu133plus2frmavecs)
# for pretty dendrograms
library(ggdendro)
# for contrast matrix
library(limma)
# pretty heatmaps
library(pheatmap) 
# annotations for the hgU1332.0Plus array
library(hgu133plus2.db) 


## ----barcode_frma_run, cache=TRUE, eval=FALSE----------------------------
 mic.frma <- frma(mic.raw, summarize="random_effect")
 save.image(file.path(resultsDir, "RDATA.frma.normalized.U133Plus2.0"))
 bc <- barcode(mic.frma)
 q()
