library(edgeR)

extract_conditions = function(col_name) {
  # this will have to get rewritten as the column names change
  conditions = unlist(strsplit(col_name, "_", fixed=TRUE))[c(1,3)]
  return(conditions)
}

colnames_to_conditions = function(col_names) {
  df = do.call(rbind, Map(extract_conditions, col_names))
  rownames(df) = NULL
  return(data.frame(df))
 }

annotate_df = function(df, out_file) {
  require(biomaRt)
  join_column = "id"
  df[,join_column] = rownames(df)
  filter_type = "ensembl_gene_id"
  ensembl_gene = "hsapiens_gene_ensembl"
  gene_symbol = "hgnc_symbol"
  ensembl = useMart("ensembl", dataset=ensembl_gene)

  a = getBM(attributes=c(filter_type, gene_symbol, "description"),
    filters=c(filter_type), values=df[,join_column],
    mart=ensembl)
  m = merge(df, a, by.x=join_column, by.y=filter_type)
  write.table(m, out_file, quote=FALSE, row.names=FALSE, sep="\t")
  return(m)
}

extract_filename = function(in_file) {
  split = unlist(strsplit(in_file, "/", fixed=TRUE))
  filename = split[length(split)]
  return(unlist(strsplit(filename, ".", fixed=TRUE))[1])
}

args = commandArgs(TRUE)
count_file = args[1]
out_prefix = paste("results/edgeR/", extract_filename(count_file),
  sep="")

counts = read.table(count_file, header=TRUE, row.names=1)
col_names_conditions = c("mirna", "batch")
conditions = colnames_to_conditions(colnames(counts))
colnames(conditions) = col_names_conditions

# pairwise comparison
cds = DGEList(counts)
design = model.matrix(~mirna+batch, conditions)
cds = calcNormFactors(cds)
cds = estimateGLMCommonDisp(cds, design)
# variance between replicates
sqrt(cds$common.dispersion)
cds = estimateGLMTrendedDisp(cds, design)
cds = estimateGLMTagwiseDisp(cds, design)
fit = glmFit(cds, design)
lrt = glmLRT(fit)
#topTags(lrt)

mds_filename = paste(out_prefix, "_mds.pdf", sep="")
pdf(mds_filename)
plotMDS(cds)
dev.off(3)

#o <- order(lrt$table$PValue)
#cpm(counts)[o[1:10],]
#de = decideTestsDGE(lrt)
#detags = rownames(counts)[as.logical(de)]
#plotSmear(lrt)

lrt_annotated = lrt
lrt_annotated$table$DE = decideTestsDGE(lrt_annotated)
annotated_filename = paste(out_prefix, "_annotated.txt", sep="")
lrt_annotated$table =  annotate_df(lrt_annotated$table,
  annotated_filename)


# things to output
# table of counts of upregulated/downregulated genes
# variance between replicates. true abundance can vary as much as 20%
# top up and down regulated genes in a table
