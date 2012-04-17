#!/usr/bin/env RScript

# Plot called variant frequencies along with filtered background frequency
# Usage:
#   plot_afback.r <in VCF file> [<detailed YAML call file>]

library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(yaml)
#library(VariantAnnotation)

args <- commandArgs(trailingOnly=TRUE)
in_vcf_file <- args[1]
in_detail_file <- args[2]
out_file <- paste(str_split(in_vcf_file, ".vcf")[[1]], "png", sep=".")
max_freq <- 2.5

#in_vcf <- readVcf(in_vcf_file, "HXB2", ScanVcfParam(fixed=NA, geno=NA, info=c("AF", "AFBACK")))
#print(info(in_vcf))

extract_info <- function(regexp) {
  llply(str_split(str_extract(in_tbl$V8, regexp), "="),
        function(x) as.numeric(x[2]))
}

read_detailed_fps <- function(in_file) {
  if (!is.na(in_file)) {
    details <- yaml.load_file(in_file)
    fps <- unlist(llply(details,
                        function(x) {
                          if (x$class == "false-positive")
                            x$pos[[2]] + 1
                        }))
  }
}

fps <- read_detailed_fps(in_detail_file)

in_tbl <- read.table(in_vcf_file, header=FALSE, sep="\t")
in_tbl$AF <- extract_info("AF=[0-9.]+")
in_tbl$AFBACK <- extract_info("AFBACK=[0-9.]+")

sub_tbl <- subset(in_tbl, select=c("V2", "AF", "AFBACK"))
names(sub_tbl) <- c("pos", "called", "background")
melt_tbl <- melt(sub_tbl, id=c("pos"), measured=c("called", "background"))
ready_tbl <- subset(melt_tbl, !is.na(melt_tbl$value))
ready_tbl <- subset(ready_tbl, ready_tbl$value < max_freq)
names(ready_tbl) <- c("position", "orig", "frequency")


ready_tbl <- ddply(ready_tbl, .(position, orig), transform,
                   variant = if (orig == "called" && position %in% fps) {
                              "false.positive"
                             } else {
                              orig})

print(head(ready_tbl))

p <- ggplot(ready_tbl, aes(x=position, y=frequency)) + geom_point(aes(color=variant)) +
     opts(legend.position="bottom", title="Low frequency variations") +
     scale_colour_hue(breaks=c("background", "called", "false.positive"),
                      labels=c("background", "called", "false positive"))
ggsave(out_file, p, width=6, height=5)
