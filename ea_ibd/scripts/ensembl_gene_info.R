#!/usr/bin/env Rscript
#
# Retrieve gene identifiers and descriptions from Ensembl biomaRt
# given mouse transcript IDs
#
# Usage:
#  ensembl_gene_info.R <out file> <comma separated list of transcript IDs>

library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
out.file <- args[1]
tx.ids <- strsplit(args[2], ",")[[1]]

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
attrs <- c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_id", "description")
filters <- c("ensembl_transcript_id")
result <- getBM(attributes=attrs, filters=filters, values=tx.ids, mart=mart)

write.table(result, file=out.file, sep=",", col.names=FALSE, row.names=FALSE)
