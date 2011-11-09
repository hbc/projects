#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <yaml config file>

library(yaml)
library(HBCshRNA)

args <- commandArgs(trailingOnly=TRUE)
config_file <- args[1]
config <- yaml.load_file(config_file)

in_data <- read.csv(config$infile, header=TRUE)
#in_data <- head(in_data, 90)
if (config$id_name == "accession") {
  work_info <- prepareByAccession(in_data, config)
} else {
  work_info <- prepareByTarget(in_data, config)
}
  
cds <- estimateVariance(work_info, config$out_base)
res_sig <- callDifferentialExpression(cds, work_info, config$out_base, config)
res_sig_merge <- mergeGenes(in_data, res_sig, config)

out_file <- paste(config$out_base, "diffexp.tsv", sep="-")
write.table(res_sig_merge$merged, file=out_file, row.names=FALSE, sep="\t")
back_file <- paste(config$out_base, "diffexp-background.tsv", sep="-")
write.table(res_sig_merge$background, file=back_file, row.names=FALSE, sep="\t")
