library(gProfileR)
library(treemap)
library(RCurl)
library(dplyr)
library(stringr)

myresContrastName <- function(res) {
  mcols(res)[2L, 2L] %>%
    str_replace("^.*:\\s", "") %>%
    str_replace_all("_", " ") %>%
    # Improve appearance for difference of differences
    str_replace_all("\\+", " \\+\n    ")
}

run_gprofiler_revigo <- function(my_results, 
                                 num_sig=200, 
                                 gp_org="mmusculus", 
                                 gp_ordered=FALSE, 
                                 gp_correction="fdr", 
                                 rev_cutoff=0.7, 
                                 rev_organism="Mus musculus")  {
  if( !file.exists("revigo.pl")){
    message("revigo.pl needs to be in the working directory")
  }
  
  contrast <- myresContrastName(my_results)
  fileStem <- snake(contrast)
  
  # get significant genes
  sig_genes <-     my_results %>% tbl_df() %>% 
    rownames_to_column(var="ensemblId") %>% 
    arrange(padj) %>% 
    dplyr::select(ensemblId) %>% 
    slice(., 1:num_sig) %>% 
    unlist() %>% 
    as.vector()
  
  #run gprofileR
  gprofiler_results <- gprofiler(query = sig_genes, 
                                 organism = gp_org, 
                                 ordered_query = gp_ordered, 
                                 exclude_iea = F, 
                                 max_set_size = 0, 
                                 correction_method = "fdr", 
                                 hier_filtering = "none", 
                                 domain_size = "annotated", 
                                 custom_bg = "")
  
  allterms <- gprofiler_results$term.id
  GOs <- allterms[grep("GO:", allterms)]
  pvals <- gprofiler_results$p.value[grep("GO:", allterms)]
  
  gprofiler_results.simplified <- gprofiler_results[, c("term.id", "term.name", "domain", "p.value", "term.size", "query.size", "overlap.size", "intersection")]
  names(gprofiler_results.simplified) <- c("term.id", "term.name",  "domain", "p.value","term.size", "query.size", "overlap.size", "assoc.gene.ids")
  
  # output to file
  write.csv(gprofiler_results.simplified, file.path(deDir, paste("gprofiler_results_", fileStem, ".csv", sep="" )))
  
  # Revigo 
  ## Parameters to change
  cutoff <- rev_cutoff #Allowed values: "0.90" "0.70" "0.50" "0.40" 
  organism <- "Mus musculus" #Allowed values: See organism.list below
  isPValue <- "yes" #Allowed values: "yes"  "no"
  whatIsBetter <- "higher" #Allowed values: "higher" "lower" "absolute" "abs_log"
  measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"
  domain <- "BP"
  
  #Do not change below
  organism.list <- list(
    "whole UniProt"=0, 
    "Homo sapiens"=9606,
    "Mus musculus"=10090,
    "Rattus norvegicus"=10116,
    "Bos taurus"=9913,
    "Gallus gallus"=9031,
    "Danio rerio"=7955,
    "Takifugu rubripes"=31033,
    "Xenopus laevis"=8355,
    "Drosophila melanogaster"=7227,
    "Caenorhabditis elegans"=6239,
    "Arabidopsis thaliana"=3702,
    "Oryza sativa"=39947,
    "Zea mays"=4577,
    "Saccharomyces cerevisiae"=4932,
    "Schizosaccharomyces pombe"=4896,
    "Dictyostelium discoideum"=44689,
    "Plasmodium falciparum"=5833,
    "Chlamydomonas reinhardtii"=3055,
    "Escherichia coli"=83333,
    "Bacillus subtilis"=1423,
    "Pseudomonas aeruginosa"=287,
    "Mycobacterium tuberculosis"=1773,
    "Mycoplasma genitalium"=2097,
    "Synechocystis sp."=1148
  )
  organism.db <- as.character(organism.list[organism])
  
  cat("**Biological Processes Domain**")
  
  domain="BP"
  mycommand=paste('revigo.pl -goterms', paste(GOs,collapse=","), '-gopvals', paste(pvals,collapse=","), '-cutoff', cutoff,  '-organism', organism.db, '-ispvalue', isPValue, '-whatisbetter', whatIsBetter, '-measure', measure, '-domain', domain,sep=" ")
  mytempfile <- tempfile()
  system2(command='perl', args=mycommand, stdout=mytempfile)
  source(mytempfile)
  
  cat("**Molecular Functions Domain**")
  
  domain="MF"
  mycommand=paste('revigo.pl -goterms', paste(GOs,collapse=","), '-gopvals', paste(pvals,collapse=","), '-cutoff', cutoff,  '-organism', organism.db, '-ispvalue', isPValue, '-whatisbetter', whatIsBetter, '-measure', measure, '-domain', domain,sep=" ")
  mytempfile <- tempfile()
  system2(command='perl', args=mycommand, stdout=mytempfile)
  source(mytempfile)
  
  cat("**Cellular Component Domain**")
  
  domain="CC"
  mycommand=paste('revigo.pl -goterms', paste(GOs,collapse=","), '-gopvals', paste(pvals,collapse=","), '-cutoff', cutoff,  '-organism', organism.db, '-ispvalue', isPValue, '-whatisbetter', whatIsBetter, '-measure', measure, '-domain', domain,sep=" ")
  mytempfile <- tempfile()
  system2(command='perl', args=mycommand, stdout=mytempfile)
  source(mytempfile)
}  


run_gprofiler <- function(my_results, 
                          num_sig=200, 
                          gp_org="mmusculus", 
                          gp_ordered=FALSE, 
                          gp_correction="fdr", 
                          rev_cutoff=0.7, 
                          rev_organism="Mus musculus")  {
  
  contrast <- myresContrastName(my_results)
  fileStem <- snake(contrast)
  
  # get significant genes
  sig_genes <-     my_results %>% tbl_df() %>% 
    rownames_to_column(var="ensemblId") %>% 
    arrange(padj) %>% 
    dplyr::select(ensemblId) %>% 
    slice(., 1:num_sig) %>% 
    unlist() %>% 
    as.vector()
  
  #run gprofileR
  gprofiler_results <- gprofiler(query = sig_genes, 
                                 organism = gp_org, 
                                 ordered_query = gp_ordered, 
                                 exclude_iea = F, 
                                 max_set_size = 0, 
                                 correction_method = "fdr", 
                                 hier_filtering = "none", 
                                 domain_size = "annotated", 
                                 custom_bg = "")
  
  allterms <- gprofiler_results$term.id
  GOs <- allterms[grep("GO:", allterms)]
  pvals <- gprofiler_results$p.value[grep("GO:", allterms)]
  
  gprofiler_results.simplified <- gprofiler_results[, c("term.id", "term.name", "domain", "p.value", "term.size", "query.size", "overlap.size", "intersection")]
  names(gprofiler_results.simplified) <- c("term.id", "term.name",  "domain", "p.value","term.size", "query.size", "overlap.size", "assoc.gene.ids")
  
  # output to file
  write.csv(gprofiler_results.simplified, file.path(deDir, paste("gprofiler_results_", fileStem, ".csv", sep="" )))
  kable(gprofiler_results.simplified)
}



plotVolcanoJH <- function(object, alpha, lfc, padj=TRUE, color, alphashade=1, pointsize=2){
  results <- as.data.frame(object) %>% tbl_df()
  if (padj==TRUE){
    results <- mutate(results, highlight=ifelse(padj<alpha & abs(log2FoldChange)>lfc, "YES", "NO"))
    ggplot(results, aes(y=-log10(padj), x=log2FoldChange, color=highlight, fill=highlight))+
      geom_point( pch=21, size=pointsize)+
      scale_color_manual(values = c("grey", color))+
      scale_fill_manual(values = alpha(c("grey", color),  alphashade))+
      guides(color=FALSE, fill=FALSE)
   } else{
     results <- mutate(results, highlight=ifelse(padj<alpha & abs(log2FoldChange)>lfc, "YES", "NO"))
     ggplot(results, aes(y=-log10(pvalue), x=log2FoldChange, color=highlight, fill=highlight))+
       geom_point( pch=21, size=pointsize)+
       scale_color_manual(values = c("grey", color))+
       scale_fill_manual(values = alpha(c("grey", color), alphashade))+
       guides(color=FALSE, fill=FALSE)
  }
}
 

plotDEGHeatmapJH <- function(
  results,
  counts,
  alpha = 0.01,
  lfc = 0,
  title = TRUE,
  ...) {
  results <- results %>%
    as.data.frame() %>%
    camel(strict = FALSE) %>%
    # Keep genes that pass alpha cutoff
    .[!is.na(.[["padj"]]), , drop = FALSE] %>%
    .[.[["padj"]] < alpha, , drop = FALSE] %>%
    # Keep genes that pass log2 fold change cutoff
    .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
    .[.[["log2FoldChange"]] > lfc |
        .[["log2FoldChange"]] < -lfc, , drop = FALSE]
  if (nrow(results) == 0) {
    warning("No genes passed significance cutoffs", call. = FALSE)
    return(NULL)
  }
  genes <- rownames(results)
  if (isTRUE(title)) {
    title <- "deg"
  } else if (is.character(title)) {
    title <- title
  }
  if (length(genes) < 2) {
    message(paste(length(genes), "is too few to plot"))
  } else {
    plotGeneHeatmapJH(counts, genes = genes, title = title, ...)
  }
} 

library(viridis)
library(pheatmap)
plotGeneHeatmapJH <- function(
  counts,
  genes = NULL,
  annotationCol = NULL,
  title = NULL,
  color = inferno(256),
  legendColor = viridis,
  # Internal parameters
  scale = "row",
  cluster_cols=TRUE,
  gap = 0,
  showColnames=TRUE,
  cutreeRows=1) {
  counts <- as.matrix(assay(counts))
  
  # Check for identifier mismatch. Do this before zero count subsetting.
  if (!is.null(genes)) {
    if (!all(genes %in% rownames(counts))) {
      stop(paste(
        "Genes missing from counts matrix:",
        toString(setdiff(genes, rownames(counts)))),
        call. = FALSE)
    }
    counts <- counts %>%
      .[rownames(.) %in% genes, , drop = FALSE]
  }
  
  # Subset zero counts
  counts <- counts %>%
    .[rowSums(.) > 0, , drop = FALSE]
  if (!is.matrix(counts) |
      nrow(counts) < 2) {
    stop("Need at least 2 genes to cluster", call. = FALSE)
  }
  
  # Convert Ensembl gene identifiers to symbol names, if necessary
  if (nrow(counts) <= 100) {
    showRownames <- TRUE
  } else {
    showRownames <- FALSE
  }
  if (isTRUE(showRownames)) {
    counts <- gene2symbol(counts)
  }
  
  # Prepare the annotation columns
  if (!is.null(annotationCol)) {
    annotationCol <- annotationCol %>%
      as.data.frame() %>%
      # Coerce annotation columns to factors
      rownames_to_column() %>%
      mutate_all(factor) %>%
      column_to_rownames()
  }
  
  # Define colors for each annotation column, if desired
  if (!is.null(annotationCol) & !is.null(legendColor)) {
    annotationColors <- lapply(
      seq_along(colnames(annotationCol)), function(a) {
        col <- annotationCol[[a]] %>%
          levels()
        colors <- annotationCol[[a]] %>%
          levels() %>%
          length() %>%
          legendColor
        names(colors) <- col
        colors
      }) %>%
      setNames(colnames(annotationCol))
  } else {
    annotationColors <- NULL
  }
  
  # pheatmap will error if `NULL` title is passed as `main`
  if (is.null(title)) {
    title <- ""
  }
  
  # If `color = NULL`, use the pheatmap default
  if (is.null(color)) {
    color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  }
  
  pheatmap(
    counts,
    annotation_col = annotationCol,
    annotation_colors = annotationColors,
    border_color = NA,
    color = color,
    main = title,
    scale = scale,
    show_rownames = showRownames,
    show_colnames = showColnames,
    cluster_cols = cluster_cols,
    gaps_col = gap,
    cutree_rows = cutreeRows)
}
