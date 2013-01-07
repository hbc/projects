library(edgeR)
pb_file = "~/cache/projects/lu_rnaseq/results/deseq/Pb/control_vs_exposed/control_vs_exposed.counts.txt"

extract_conditions = function(col_name) {
  # this will have to get rewritten as the column names change
  field_string = unlist(strsplit(col_name, ".", fixed=TRUE))[5]
  conditions = unlist(strsplit(field_string, "_", fixed=TRUE))[1:3]
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

col_names_pb = c("control_2", "control_1", "control_3",
  "exposed_2", "exposed_1", "exposed_3")
col_names_conditions = c("sample", "treatment", "trial")

pb = read.table(pb_file, header=TRUE, row.names=1)
conditions = colnames_to_conditions(colnames(pb))
colnames(pb) = c("control_2", "control_1", "control_3",
          "exposed_2", "exposed_1", "exposed_3")
colnames(conditions) = c("sample", "treatment", "trial")

# pairwise comparison for pb
cds = DGEList(pb)
design = model.matrix(~trial+treatment, conditions)
cds = estimateGLMCommonDisp(cds, design)
cds = estimateGLMTrendedDisp(cds, design)
cds = estimateGLMTagwiseDisp(cds, design)
fit = glmFit(cds, design)
lrt = glmLRT(fit)
topTags(lrt)

plotMDS(cds)
pdf("Pb_mds.pdf")
plotMDS(cds)
dev.off(3)

o <- order(lrt$table$PValue)
cpm(pb)[o[1:10],]
de = decideTestsDGE(lrt)
detags = rownames(pb)[as.logical(de)]
plotSmear(lrt)


pbn_file = "~/cache/projects/lu_rnaseq/results/deseq/PbN/control_vs_exposed/control_vs_exposed.counts.txt"
col_names_pbn = c("control_2", "control_1", "control_3",
  "exposed_2", "exposed_1", "exposed_3")
col_names_conditions = c("sample", "treatment", "trial")

pbn = read.table(pbn_file, header=TRUE, row.names=1)
conditions = colnames_to_conditions(colnames(pbn))
colnames(pbn) = col_names_pbn
colnames(conditions) = c("sample", "treatment", "trial")

# pairwise comparison for pb
cds_pbn = DGEList(pbn)
design_pbn = model.matrix(~trial+treatment, conditions)
cds_pbn = estimateGLMCommonDisp(cds_pbn, design_pbn)
cds_pbn = estimateGLMTrendedDisp(cds_pbn, design_pbn)
cds_pbn = estimateGLMTagwiseDisp(cds_pbn, design_pbn)
fit_pbn= glmFit(cds_pbn, design_pbn)
lrt_pbn = glmLRT(fit_pbn)
topTags(lrt_pbn)

plotMDS(cds_pbn)
pdf("PbN_mds.pdf")
plotMDS(cds_pbn)
dev.off(3)

o_pbn <- order(lrt_pbn$table$PValue)
cpm(pbn)[o[1:10],]
de_pbn = decideTestsDGE(lrt_pbn)
detags = rownames(pbn)[as.logical(de_pbn)]
plotSmear(lrt)

lrt_annotated = lrt
lrt_annotated$table$DE = decideTestsDGE(lrt_annotated)
lrt_annotated$table =  annotate_df(lrt_annotated$table, "Pb_annotated.txt")
lrt_pbn_annotated = lrt_pbn
lrt_pbn_annotated$table$DE = decideTestsDGE(lrt_pbn_annotated)
lrt_pbn_annotated$table =  annotate_df(lrt_pbn_annotated$table,
  "PbN_annotated.txt")
