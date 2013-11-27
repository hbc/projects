#!/usr/bin/env RScript

# Plot true and false positive rates for randomly downsampled reads for
# minority variant assessment.

library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- str_c(str_split(in_file, "[.]")[[1]][1], ".pdf")
orig_read_count <- 120000

d <- read.csv(in_file, header=TRUE)
d$tpr <- d$true.positive / (d$true.positive + d$false.negative)
d$fpr <- d$false.positive /  (d$false.positive + d$true.negative)
d$reads <- d$downsample * orig_read_count

d2 <- ddply(d, .(reads),
            .fun = function(df) {
              data.frame(true.positive.rate = mean(df$tpr),
                         false.positive.rate = mean(df$fpr),
                         tpr.std = sd(df$tpr),
                         fpr.std = sd(df$fpr))
            })
print(d2)

d3 <- melt(d2, id=c("reads", "tpr.std", "fpr.std"),
           measured=c("true.positive.rate", "false.positive.rate"))
d3$std <- ifelse(d3$variable=="true.positive.rate", d3$tpr.std, d3$fpr.std)
d3$ymin <- d3$value - d3$std
d3$ymax <- d3$value + d3$std
print(d3)

facet_labeller <- function(var, val) {
  ifelse(val=="true.positive.rate", "True positive rate",
         "False positive rate")}
p <- ggplot(d3, aes(x=reads, y=value, ymin=ymin, ymax=ymax)) + geom_line() +
  geom_errorbar() + facet_grid(variable ~ ., scales="free", labeller=facet_labeller) +
  opts(axis.title.y = theme_blank()) + xlab("Downsampled reads")
ggsave(out_file, p, width=6, height=6)

