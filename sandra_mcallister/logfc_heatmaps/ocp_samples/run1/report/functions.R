library(gProfileR)
library(treemap)
library(RCurl)
library(dplyr)
library(stringr)
library(knitr)
library(rio)
library(pheatmap)

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
    arrange(desc(abs(log2FoldChange))) %>% 
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
  


plotDEGHeatmapJH <- function(res_object, mylfc, myalpha, bcb_object, rld_object){
  plotgenes <- res_object %>% as.data.frame() %>% 
    tibble::rownames_to_column("ensgene")%>% 
    arrange(padj)  %>% 
    select(ensgene,log2FoldChange, padj) %>% 
    filter(abs(log2FoldChange)>mylfc & padj<myalpha) %>% 
    select(ensgene) %>%
    unlist() %>%  as.vector()
  myrld <- assay(rld_object)[plotgenes,, drop=FALSE]
    # get gene symbols
  mysymbols <- gene2symbol(bcb_object)
  row.names(myrld) <- mysymbols$symbol[match(row.names(myrld), mysymbols$ensgene)]
  pheatmap(myrld, scale="row", cluster_cols = TRUE, color=viridis(200))
}

output_results <- function(res_object,  bcb_object, keep_stat=TRUE, deDir="results/differential_expression", ngenes=20){
  myresults <- res_object %>% as.data.frame() %>% tibble::rownames_to_column("ensgene")
  if(keep_stat==TRUE) {
    myresults <- myresults %>% arrange(desc(padj))
  } else {
    myresults <- myresults %>% select(-pvalue, -padj, -stat) %>% arrange(abs(log2FoldChange))
  }
  # get annotations for genes
  myresults <- inner_join(myresults, metadata(bcb_object)$annotable, by="ensgene")
  
  print(kable(myresults[1:ngenes,]))
  comparison <- sub("Wald test p-value: ", "", res_object@elementMetadata[,2][5]) %>% gsub(" ", "_",.)
  myfilename <- paste(comparison, "all.xlsx", sep="")
  export(myresults, file=file.path(deDir, myfilename), format="xlsx")
}
