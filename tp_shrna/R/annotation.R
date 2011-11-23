# Prepare annotation data on genomic features using BioMart

library(IRanges)
library(ChIPpeakAnno)
#library(rtracklayer)
library(biomaRt)
library(RMySQL)

#' Retrieve annotations from a BioMart instance into GFF files.
#' Utilizes either ChIPpeakAnno default retrieval or a custom
#' function for non-coding RNAs.
#' @export
downloadMartData <- function(ftype, dataset, out_dir=NULL) {
  mart <- useMart("ensembl", dataset=dataset)
  if (!is.null(out_dir)) {
    save.file <- paste(out_dir,
                       paste(dataset, "-", ftype, ".gff3", sep=""),
                       sep="/")
  } else {
    save.file <- NULL
  }
  if (!is.null(save.file) && file.exists(save.file)) {
    mart.data <- import(save.file)
  } else {
    if (ftype %in% c("ncRNA", "lincRNA")) {
      mart.data <- getNcData(mart, ftype)
    } else {
      mart.data <- getAnnotation(mart, ftype)
    }
    export(mart.data, save.file)
  }
  list(data=mart.data, file=save.file)
}

#' Convert Ensembl GRCh37 chromosome names (1, 2) to
#' UCSC chromosome names (chr1, chr2) using the mapping
#' table provided by UCSC.
#' in_rd is a ranged data object.
convertChrNamesToUcsc <- function(in_rd) {
  ucsc_dbhost = "genome-mysql.cse.ucsc.edu"
  ucsc_dbuser = "genome"
  ucsc_dbname = "hg19"
  con = dbConnect(MySQL(), user=ucsc_dbuser, host=ucsc_dbhost,
                  dbname=ucsc_dbname)
  data = dbGetQuery(con, "select * from ucscToEnsembl")
  print(data)
}

#' Retrieve non-coding annotations not covered by ChIPpeakAnno
#' http://uswest.ensembl.org/info/docs/genebuild/ncrna.html
getNcData <- function(mart, ftype) {
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  other.nc <- c("IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
                "Mt_tRNA", "Mt_tRNA_pseudogene", "TR_V_pseudogene")
  if (ftype == "ncRNA") {
    nc.types <- c("lincRNA", "miRNA_pseudogene", "misc_RNA_pseudogene",
                  "polymorphic_pseudogene", "pseudogene",
                  "snoRNA", "snoRNA_pseudogene", "snRNA_pseudogene",
                  "processed_transcript")
  } else if (ftype == "lincRNA") {
    nc.types <- c("lincRNA")
  }
  attrs <- c("ensembl_gene_id", "chromosome_name", "start_position",
             "end_position", "strand", "description", "gene_biotype")
  result <- getBM(attributes=attrs, filters=c("biotype"), values=nc.types, mart=mart)
  # remove duplicated
  result <- unique(result)
  result <- result[order(result[,3]), ]              
  duplicated.id <- result[duplicated(result[,1]), 1]
  result <- result[!duplicated(result[,1]), ]
  result.rd <- RangedData(IRanges(start=as.numeric(result[,3]), end=as.numeric(result[,4]),
                                  names=as.character(result[,1])),
                          strand=result[,5], description=as.character(result[,6]),
                          space=as.character(result[,2]), biotype=as.character(result[,7]))
}
