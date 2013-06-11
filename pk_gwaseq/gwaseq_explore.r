#!/usr/bin/env RScript

# Exploratory plots of GWASeq combined results

library(ggplot2)
library(knitr)

getProgram <- function(args) {
  sub("--file=", "", args[grep("--file=", args)])
}
args <- commandArgs(trailingOnly=TRUE)
program <- getProgram(commandArgs(trailingOnly=FALSE))
in_file <- args[1]
low_variant_cutoff <- 10
qual_cutoff <- 10000
base_name <- "gwaseq_explore"

d <- read.csv(in_file, header=TRUE)
d$total_variants <- d$HET + d$HOM_VAR
d$is_low <- ifelse(d$total_variants < low_variant_cutoff, "low", "high")
d$is_bait <- ifelse(d$bait == 1, "bait", "nonbait")
d$qual_trunc <- ifelse(d$qual < qual_cutoff, d$qual, qual_cutoff)
d$pct_hets <- d$HET / d$total_variants
d.low <- subset(d, d$total_variants < low_variant_cutoff)

knit(paste(file.path(dirname(program), base_name), ".Rmd", sep=""))
header <- paste(file.path(dirname(program), base_name), "-header.tex", sep="")
system(paste("markdown2pdf -H ", header, " ", base_name, ".md", sep=""))
