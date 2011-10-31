#!/usr/bin/env RScript

# Plot spread of target RNA counts for single accession to examine combining
# Usage:
#   plot_target_spread.r <in_file> <out_base>

library(HBCshRNA)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
out_base <- args[2]

in_data <- read.csv(infile, header=TRUE)
reorg_data <- loadByTarget(in_data)
out_file <- paste(out_base, "spread.pdf", sep="-")
plotSpreadDistribution(reorg_data, out_file)
