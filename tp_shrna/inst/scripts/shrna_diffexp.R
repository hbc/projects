#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <yaml config file>

library(yaml)
library(hbc.shrna)

args <- commandArgs(trailingOnly=TRUE)
config_file <- args[1]
config <- yaml.load_file(config_file)

in_data <- read.csv(config$infile, header=TRUE)
if (!is.null(in_data$Replicate)) {
  in_data <- subset(in_data, (!in_data$Replicate %in% config$drop_replicates))
}
if (config$id_name == "accession") {
  work_info <- prepareByAccession(in_data, config)
} else {
  work_info <- prepareByTarget(in_data, config)
}

writeTsvOutput <- function(res_sig, in_data, out_base, config) {
  res_sig_merge <- mergeGenes(in_data, res_sig, config)

  out_file <- paste(out_base, "diffexp.tsv", sep="-")
  write.table(res_sig_merge$merged, file=out_file, row.names=FALSE, sep="\t")
  back_file <- paste(out_base, "diffexp-background.tsv", sep="-")
  write.table(res_sig_merge$background, file=back_file, row.names=FALSE, sep="\t")
}
  
ready <- estimateVariance(work_info, config)
called <- callDifferentialExpression(ready$cds, ready$cds_multi, work_info, config)
writeTsvOutput(called$res, in_data, config$out_base, config)
if (!is.null(called$res_multi)) {
  writeTsvOutput(called$res_multi, in_data, paste(config$out_base, "multi", sep="-"), config)
}
