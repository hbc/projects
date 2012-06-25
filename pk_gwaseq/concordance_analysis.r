#!/usr/bin/env RScript

# Concordance analysis: subset inputs + cumulative frequency plot

library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- args[2]

con_df_orig <- read.csv(in_file, header=TRUE)
con_df <- subset(con_df_orig, select=c("sample", "callable_concordance",
                                "concordant", "nonref_concordant"))
print(summary(con_df))
print(summary(subset(con_df, callable_concordance < 90.0)))
write.table(con_df, file=paste(out_file, ".csv", sep=""), row.names=FALSE, sep=",")

# Cumulative frequency
con_df_all <- summarize(con_df, callable_concordance = unique(callable_concordance),
                        cdf = ecdf(callable_concordance)(unique(callable_concordance))
                              * length(callable_concordance))
print(head(con_df_all))
p <- ggplot(con_df_all, aes(cdf, callable_concordance)) + geom_point() +
            opts(title="Callable Concordance")
ggsave(paste(out_file, "-cdf.png", sep=""), p, width=5, height=5)
