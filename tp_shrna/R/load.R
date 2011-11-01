#' Load and reorganize input data for analysis
#' This has lots of variable references for the current experiment
#' and should be generalized for other inputs.

library(plyr)
library(reshape2)


#' Filter input data requiring minimum read counts
#' @imports plyr
filterDfByCounts <- function(in_df, min_count){
  with_flag <- ddply(in_df, .(shrna.id),
                     .fun = function(df) {
                       transform(df,
                                 keep = max(df$d.3) >= min_count)
                     })
  filter_df <- subset(with_flag, keep == TRUE)
  filter_df$keep <- NULL
  filter_df
}

# Custom cleaning of input headers
cleanData <- function(in_data, min_count) {
  names(in_data)[2:4] = c("d.3", "w.3", "accession")
  names(in_data) <- tolower(names(in_data))
  if (!is.null(min_count)) {
    in_data <- filterDfByCounts(in_data, min_count)
  }
  in_data
}

#' Prepare input file in data.frame suitable for DEseq
#' Does 3 transformations:
#'  1. Collapse multiple probes targetting the same acession by summing
#'  2. Reshape the data.frame so columns correspond to replicates and conditions
#'  3. Removes any samples with missing values.
#' This is specific for input data file and would need to be generalized for
#' other inputs.
#' @export
#' @imports plyr, reshape2
prepareInputs <- function(in.data, min_count = NULL) {
  in.data <- cleanData(in.data, min_count)
  #in.data <- head(in.data, 90)
  #print(head(in.data))
  collapsed.data <- ddply(in.data, .(accession, replicate), function (df)
                          data.frame(sum.d.3 = sum(df$d.3),
                                     sum.w.3 = sum(df$w.3),
                                     gene.symbol=df$gene.symbol[1]))
  #print(head(collapsed.data))
  collapsed.data$gene.symbol <- NULL
  melt.data <- melt(collapsed.data, id=c("accession", "replicate"),
                    measured=c("sum.d.3", "sum.w.3"))
  reshape.data <- dcast(melt.data, accession ~ variable + replicate)
  row.names(reshape.data) <- reshape.data$accession
  reshape.data$accession <- NULL
  reshape.data <- na.exclude(reshape.data)
  #print(head(reshape.data))
  list(counts = reshape.data,
       conditions = c(rep("d.3", ncol(reshape.data) / 2),
                      rep("w.3", ncol(reshape.data ) / 2)))
}

#' Organize input data table to examine individual shRNAs by accession
#' The goal is to assess target variability to determine a best pooling
#' strategy.
#' Different stratifications around accession:
#'  targets (3) -> keep as different points
#'  replicates (3) -> collapse with mean,
#'  positions(2) -> separate as two different comparisons
#' After organizing around these, then collapse based on percentage variation.
#' @export
#' @imports plyr, reshape2
loadByTarget <- function(in.data, min_count = NULL) {
  in.data <- cleanData(in.data, min_count)
  # first collapse the replicates
  norep.data <- ddply(in.data, .(shrna.id), function(df)
                      data.frame(mean.A = mean(df$d.3),
                                 mean.B = mean(df$w.3),
                                 accession = df$accession[1]))
  melt.data <- melt(norep.data, id = c("shrna.id", "accession"),
                    measured = c("mean.A", "mean.B"))
  percent.data <- ddply(melt.data, .(accession), function(df)
                        data.frame(spread = (max(df$value) - min(df$value)) /
                                             max(df$value)))
}
