# Load libraries
library(biomaRt)

# Set directories
baseDir <- "~/projects/ab_nucmito/"
dataDir <- file.path(baseDir, "data")
resultsDir <- file.path(baseDir, "results")
metaDir <- file.path(baseDir, "meta")


# Load data
annots.450k <- read.csv(file.path(dataDir, "14_April_2014_Updated_humanmethylation450_15017482_v1-2.csv"), skip=7)
annots.mito <- read.delim(file.path(dataDir, "mitoproteome.tsv"))

# Find genomic locations of mitochondrial genes using biomaRt

## setup biomart
ensemblmart <-  useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensemblmart)
filters <- listFilters(ensemblmart)

## pull down chromosomal positions using uniprot ids from mito annotations
mito.pos <- getBM(annots.mito$ACCESSION,filters="uniprot_swissprot_accession", attributes=c("uniprot_swissprot_accession","chromosome_name","start_position", "end_position", "strand"), mart=ensemblmart)

# subset to human chromosomes, dicard all mitochondrial, patches and polymorphic regions
chrs <- c(seq(1,22), "X", "Y")
mito.pos <- mito.pos[mito.pos$chromosome_name %in% chrs,]
