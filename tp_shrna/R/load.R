library(plyr)
library(reshape2)

# Prepare input file in data.frame suitable for DEseq
# Does 3 transformations:
#  1. Collapse multiple probes targetting the same acession by summing
#  2. Reshape the data.frame so columns correspond to replicates and conditions
#  3. Removes any samples with missing values.
# This is specific for input data file and would need to be generalized for
# other inputs.
prepareInputs <- function(infile) {
  in.data <- read.csv(infile, header=TRUE)
  #in.data <- head(in.data, 90)
  names(in.data)[2:4] = c("counts.d.3", "counts.w.3", "accession")
  names(in.data) <- tolower(names(in.data))
  #print(head(in.data))
  collapsed.data <- ddply(in.data, .(accession, replicate), function (df)
                        data.frame(sum.d.3 = sum(df$counts.d.3),
                                   sum.w.3 = sum(df$counts.w.3),
                                   gene.symbol=df$gene.symbol[1]))
  #print(head(collapsed.data))
  collapsed.data$gene.symbol <- NULL
  melt.data <- melt(collapsed.data, id=c("accession", "replicate"),
                    measured=c("sum.d.3", "sum.w.3"))
  reshape.data <- dcast(melt.data, accession ~ variable + replicate)
  row.names(reshape.data) <- reshape.data$accession
  reshape.data$accession <- NULL
  reshape.data <- na.exclude(reshape.data)
  print(head(reshape.data))
  list(counts = reshape.data,
       conditions = c(rep("d.3", ncol(reshape.data) / 2), rep("w.3", ncol(reshape.data ) / 2)))
}
