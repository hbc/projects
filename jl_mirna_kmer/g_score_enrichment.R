library(biomaRt)
library(ggplot2)

ensembl_gene = "hsapiens_gene_ensembl"
filter_type = "refseq_mrna"

work_dir = "/Users/rory/cache/g-quad"
g_scores_fn = "human-canonical-gscores.txt"
knowngene_refseq = "knowngene-to-refseq.txt"
oday_fn = "oday-lin28-binding.csv"
hefner_fn = "hefner-lin28-binding.csv"
cho_fn = "lin28-cho.csv"
tan_fn = "lin28-tan.csv"

setwd(work_dir)

ensembl_transcript_id_to_refseq = function(d) {
	require(biomaRt)
	ensembl = useMart('ensembl', dataset = ensembl_gene)
	a = getBM(attributes=c("ensembl_transcript_id", "refseq_mrna"),
      filters=c("ensembl_transcript_id"), values=d[,"id"],
      mart=ensembl)
	m = merge(d, a, by.x='id', by.y="ensembl_transcript_id")
	return(m)
}

entrez_id_to_refseq = function(d) {
	require(biomaRt)
	ensembl = useMart('ensembl', dataset = ensembl_gene)
	a = getBM(attributes=c("entrezgene", "refseq_mrna"),
      filters=c("entrezgene"), values=d[,"EntrezID"],
      mart=ensembl)
	m = merge(d, a, by.x='EntrezID', by.y="entrezgene")
	return(m)
}

add_biotype = function(d) {
    require(biomaRt)
    ensembl = useMart('ensembl', dataset=ensembl_gene)
    a = getBM(attributes=c("refseq_mrna", "gene_biotype"),
        filters = c("refseq_mrna"), values=d[,"hg19.kgXref.refseq"],
        mart=ensembl)
    m = merge(d, a, by.x="hg19.kgXref.refseq", by.y="refseq_mrna")
    return(m)
}

load_scores_dataset = function(filename) {
    x = read.table(filename, header=TRUE, sep="\t")
    x = x[complete.cases(x),]
    return(x)
}

load_hfner_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    x$id = x$TranscriptID
    x = ensembl_transcript_id_to_refseq(x)
    x = subset(x, "Sum.T.to.C.on.transcript" > 20)
    return(x[complete.cases(x),])
}

load_day_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    return(x[complete.cases(x),])
}

load_cho_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    x$id = x$Accession
    return(x[complete.cases(x),])
}

load_tan_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    x = entrez_id_to_refseq(x)
    x$id = x$refseq_mrna
    return(x[complete.cases(x),])
}

load_kg_to_refseq = function(filename) {
    x = read.table(filename, header=TRUE, sep="\t")
    x = x[complete.cases(x),]
    x = subset(x, hg19.kgXref.refseq != "")
    x = add_biotype(x)
    x = subset(x, gene_biotype == "protein_coding")
    return(x)
}

hfner = load_hfner_dataset(hefner_fn)
oday = load_oday_dataset(oday_fn)
kg_to_refseq = load_kg_to_refseq(knowngene_refseq)
cho = load_cho_dataset(cho_fn)
tan = load_tan_dataset(tan_fn)
scores = load_scores_dataset(g_scores_fn)

scores = merge(scores, kg_to_refseq, by.x="id", by.y="hg19.kgXref.kgID")
scores$id = scores$hg19.kgXref.refseq
scores = subset(scores, hg19.kgXref.refseq != "")


## oday analysis

oday = oday[, c("Gene.ID", "FoldChange", "WaldStat")]

in_oday = intersect(oday$Gene.ID, scores$hg19.kgXref.refseq)
not_in_oday = setdiff(scores$hg19.kgXref.refseq, oday$Gene.ID)

oday_scores = merge(oday, scores, by.x="Gene.ID", by.y="hg19.kgXref.refseq")
not_in_oday_scores = subset(scores, hg19.kgXref.refseq %in% not_in_oday)

oday_is_greater = function(oday_scores, not_in_oday_scores) {
    oday_mean = mean(oday_scores$length_normalized_score)
    s = sample(not_in_oday_scores$length_normalized_score,
        length(oday_scores$length_normalized_score))
    not_in_oday_mean = mean(s)
    return(oday_mean > not_in_oday_mean)
}
samples = replicate(10000, oday_is_greater(oday_scores, not_in_oday_scores))
p_value = 1 - sum(samples) / length(samples)
paste("Oday lin28a p-value that G-scores of lin28a binding transcripts are greater than the average G-score of protein coding transcripts in general:", p_value)

quartz()
nonzero_scores = subset(oday_scores, length_normalized_score > 0)
p = ggplot(nonzero_scores,
    aes(length_normalized_score, FoldChange)) + geom_point() +
    geom_smooth(method=lm)
p
ggsave("oday.pdf")

## hfner analysis
in_hfner = intersect(hfner$refseq_mrna, scores$hg19.kgXref.refseq)
not_in_hfner = setdiff(scores$hg19.kgXref.refseq, hfner$refseq_mrna)

hfner_scores = merge(hfner, scores, by.x="refseq_mrna", by.y="hg19.kgXref.refseq")
not_in_hfner_scores = subset(scores, hg19.kgXref.refseq %in% not_in_hfner)

hfner_is_greater = function(oday_scores, not_in_oday_scores) {
    oday_mean = mean(oday_scores$length_normalized_score)
    s = sample(not_in_oday_scores$length_normalized_score,
        length(oday_scores$length_normalized_score))
    not_in_oday_mean = mean(s)
    return(oday_mean > not_in_oday_mean)
}

samples = replicate(10000, hfner_is_greater(hfner_scores, not_in_hfner_scores))
p_value = 1 - sum(samples) / length(samples)
paste("Hefner lin28a p-value that G-scores of lin28a binding transcripts are greater than the average G-score of protein coding transcripts in general:", p_value)

nonzero_scores = subset(hfner_scores, length_normalized_score > 0)
p = ggplot(nonzero_scores,
    aes(length_normalized_score, ConversionEventCount)) + geom_point() +
    geom_smooth(method=lm)
p
ggsave("hfner.pdf")

# cho analysis
#in_cho = intersect(cho$id, scores$id)
#not_in_cho = setdiff(scores$id, cho$id)
#cho_scores = merge(cho, scores, by.x="id", by.y="id")

# tan analysis
in_tan = intersect(tan$id, scores$id)
not_in_tan = setdiff(scores$id, tan$id)
tan_scores = merge(tan, scores, by="id")
not_in_tan_scores = subset(scores, id %in% not_in_tan)
samples = replicate(10000, hfner_is_greater(tan_scores, not_in_tan_scores))
p_value = 1 - sum(samples) / length(samples)
