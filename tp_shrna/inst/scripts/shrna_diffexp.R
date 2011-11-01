#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <in_file> <out_base>

library(HBCshRNA)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
out.base <- args[2]

config <- c(min_count = 500)
in_data <- read.csv(infile, header=TRUE)
in.data <- prepareInputs(in.data, config$min_count)
cds <- estimateVariance(in.data, out.base)
res <- callDifferentialExpression(cds, in.data, out.base)
