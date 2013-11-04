# script to do some bootstrap testing of the hypothesis that the lin28a
# targets are more G-quadraplex rich than a random sample of the
# transcripts in the genome

library(biomaRt)
library(ggplot2)

ensembl_gene = "hsapiens_gene_ensembl"
filter_type = "refseq_mrna"

work_dir = "/Users/rory/cache/g-quad"
g_scores_fn = "human-canonical-gscores.txt"
knowngene_refseq = "knowngene-to-refseq.txt"
oday_fn = "oday-lin28-binding.csv"
hefner_fn = "hefner-lin28-binding.csv"

setwd(work_dir)

load_scores_dataset = function(filename) {
    x = read.table(filename, header=TRUE, sep="\t")
    x = x[complete.cases(x),]
    x = subset(x, hg19.kgXref.refseq != "")
    return(x)
}

load_hfner_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    return(x[complete.cases(x),])
}

load_oday_dataset = function(filename) {
    x = read.csv(filename, header=TRUE, sep=",")
    return(x[complete.cases(x),])
}

load_kg_to_refseq = function(filename) {
    x = read.table(filename, header=TRUE, sep="\t")
    return(x[complete.cases(x),])
}

hfner = load_hfner_dataset(hefner_fn)
oday = load_oday_dataset(oday_fn)
kg_to_refseq = load_kg_to_refseq(knowngene_refseq)
scores = load_scores_dataset(g_scores_fn)

scores = merge(scores, kg_to_refseq, by.x="id", by.y="hg19.kgXref.kgID")

oday = oday[, c("Gene.ID", "FoldChange", "WaldStat")]

in_oday = intersect(oday$Gene.ID, scores$hg19.kgXref.refseq)
not_in_oday = setdiff(scores$hg19.kgXref.refseq, oday$Gene.ID)

oday_scores = merge(oday, scores, by.x="Gene.ID", by.y="hg19.kgXref.refseq")
not_in_oday_scores = subset(scores, hg19.kgXref.refseq %in% not_in_oday)

not_in_oday_scores$log_score = log(not_in_oday_scores$length_normalized_score)
ggplot(not_in_oday_scores, aes(log_score)) + geom_histogram()

oday_scores$log_score = log(oday_scores$length_normalized_score)
ggplot(oday_scores, aes(log_score)) + geom_histogram()


t.test(oday_scores$length_normalized_score, not_in_oday_scores$length_normalized_score,
       alternative="greater")

mean_oday = mean(oday_scores$length_normalized_score)
oday_is_greater = function(oday_scores, not_in_oday_scores) {
    oday_mean = mean(oday_scores$length_normalized_score)
    s = sample(not_in_oday_scores$length_normalized_score,
        length(oday_scores$length_normalized_score))
    not_in_oday_mean = mean(s)
    return(oday_mean > not_in_oday_mean)
}
samples = replicate(10000, oday_is_greater(oday_scores, not_in_oday_scores))
p_value = 1 - sum(samples) / length(samples)

## library(boot)
## oday_greater = function(data, indices) {
##     oday_mean = mean(oday_scores$length_normalized_score)
##     d = data[indices,]
##     return(mean(d$length_normalized_score) < oday_mean)
## }
## results = boot(not_in_oday_scores, statistic=oday_greater, R=10)

# "refseq_mrna"
# "ensembl_transcript_id"
# "ensembl_gene_id"
ensembl_transcript_id_to_refseq = function(d) {
	require(biomaRt)
	ensembl = useMart('ensembl', dataset = ensembl_gene)
	a = getBM(attributes=c("ensembl_transcript_id", "refseq_mrna"),
      filters=c("ensembl_transcript_id"), values=d[,"id"],
      mart=ensembl)
	m = merge(d, a, by.x='id', by.y="ensembl_transcript_id")
	return(m)
}

hfner$id = hfner$TranscriptID
