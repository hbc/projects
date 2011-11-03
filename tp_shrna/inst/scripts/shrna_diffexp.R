#!/usr/bin/env RScript

# Calculate differential expression of shortRNA experimental data using DEseq.
# Usage:
#   shrna_diffexp.R <in_file> <out_base>

library(HBCshRNA)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
out_base <- args[2]

#config <- list(min_count = 500, id_name="shrna.id",
#               fdr_thresh = 0.1)
config <- list(min_count = 500, id_name="accession",
               fdr_thresh = 0.1)

in_data <- read.csv(infile, header=TRUE)
#in_data <- head(in_data, 90)
if (config$id_name == "accession") {
  work_info <- prepareByAccession(in_data, config)
} else {
  work_info <- prepareByTarget(in_data, config)
}
  
cds <- estimateVariance(work_info, out_base)
res_sig <- callDifferentialExpression(cds, work_info, out_base, config)
res_sig_genes <- mergeGenes(in_data, res_sig, config)

out_file <- paste(out_base, "diffexp.tsv", sep="-")
write.table(res_sig_genes, file=out_file, row.names=FALSE, sep="\t")
