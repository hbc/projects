#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <in_file> <out_base>

library(HBCshRNA)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
out.base <- args[2]

in.data <- prepareInputs(infile)
cds <- estimateVariance(in.data, out.base)
res <- callDifferentialExpression(cds, in.data, out.base)
