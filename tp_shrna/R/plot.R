library(ggplot2)

#' Plot histogram distribution of spread values
#' @export
#' @imports ggplot2
plotSpreadDistribution <- function(in_data, out_file) {
  p <- ggplot(in_data, aes(x=spread)) + geom_histogram(binwidth=0.05) +
       xlab("Percentage count spread for multiple targets at an accession")
  if (!is.null(out_file)) {
    ggsave(out_file, p, width=6, height=6)
  }
  p
}
