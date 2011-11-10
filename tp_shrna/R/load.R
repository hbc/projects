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
prepareByAccession <- function(in_data, config) {
  in_data <- cleanData(in_data, config$min_count)
  config$model <- data.frame(config$model)
  config$model$condition <- as.character(config$model$condition)
  in_data$gene.symbol <- NULL
  in_data$log2..3w.3d. <- NULL
  #print(head(in_data))
  collapsed.data <- ddply(in_data, .(accession, replicate), function (df)
                          data.frame(sum.1 = sum(df[[config$model$condition[1]]]),
                                     sum.2 = sum(df[[tail(config$model$condition, 1)[1]]])))
  #print(head(collapsed.data))
  melt.data <- melt(collapsed.data, id=c("accession", "replicate"),
                    measured=c("sum.1", "sum.2"))
  reshape.data <- dcast(melt.data, accession ~ variable + replicate)
  row.names(reshape.data) <- reshape.data$accession
  reshape.data$accession <- NULL
  row.names(config$model) <- names(reshape.data)
  reshape.data <- na.exclude(reshape.data)
  #print(head(reshape.data))
  list(counts = reshape.data,
       model = config$model)
}

#' Prepare input data organized by target instead of accession
#' @export
#' @imports plyr, reshape2
prepareByTarget <- function(in_data, config) {
  in_data <- cleanData(in_data, config$min_count)
  config$model <- data.frame(config$model)
  config$model$condition <- as.character(config$model$condition)
  in_data$gene.symbol <- NULL
  in_data$log2..3w.3d. <- NULL
  in_data$accession <- NULL
  melt.data <- melt(in_data, id=c("shrna.id", "replicate"),
                    measured=c(config$model$condition[1],
                               tail(config$model$condition, 1)[1]))
  reshape.data <- dcast(melt.data, shrna.id ~ variable + replicate)
  row.names(reshape.data) <- reshape.data$shrna.id
  reshape.data$shrna.id <- NULL
  row.names(config$model) <- names(reshape.data)
  reshape.data <- na.exclude(reshape.data)
  list(counts = reshape.data,
       model = config$model)
}

#' Prepare data table mapping shRNA ids to accession and gene symbols
#' @export
mergeGenes <- function(orig_data, remap_data, config = NULL) {
  orig_data <- cleanData(orig_data, config$min_count)
  scols <- c("accession", "gene.symbol")
  if (config$id_name == "shrna.id") {
    scols <- c("shrna.id", scols)
  }
  genemap <- unique(subset(orig_data, select = scols))
  merged_data <- merge(remap_data, genemap)
  sorted_data <- merged_data[order(merged_data$pval), ]
  list(merged=sorted_data,
       background=genemap)
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
