#!/usr/bin/env RScript

# Plot called variant frequencies along with filtered background frequency
# Usage:
#   plot_afback.r <in VCF file>

library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
#library(VariantAnnotation)

args <- commandArgs(trailingOnly=TRUE)
in_vcf_file <- args[1]
out_file <- paste(str_split(in_vcf_file, ".vcf")[[1]], "pdf", sep=".")
max_freq <- 2.5

#in_vcf <- readVcf(in_vcf_file, "HXB2", ScanVcfParam(fixed=NA, geno=NA, info=c("AF", "AFBACK")))
#print(info(in_vcf))

extract_info <- function(regexp) {
  llply(str_split(str_extract(in_tbl$V8, regexp), "="),
        function(x) as.numeric(x[2]))
}

in_tbl <- read.table(in_vcf_file, header=FALSE, sep="\t")
in_tbl$AF <- extract_info("AF=[0-9.]+")
in_tbl$AFBACK <- extract_info("AFBACK=[0-9.]+")

sub_tbl <- subset(in_tbl, select=c("V2", "AF", "AFBACK"))
names(sub_tbl) <- c("pos", "called", "background")
melt_tbl <- melt(sub_tbl, id=c("pos"), measured=c("called", "background"))
ready_tbl <- subset(melt_tbl, !is.na(melt_tbl$value))
ready_tbl <- subset(ready_tbl, ready_tbl$value < max_freq)
names(ready_tbl) <- c("position", "variant", "frequency")

print(head(ready_tbl))

p <- ggplot(ready_tbl, aes(x=position, y=frequency)) + geom_point(aes(color=variant)) +
     opts(legend.position="bottom", title="Low frequency variations")
ggsave(out_file, p, width=6, height=4)
