library(ggplot2)
library(reshape)
library(googlevis)
library(stringr)
library(edger)
library(tools)
wd = '/n/home05/kirchner/cache/projects/kobzik_influenza/scripts/'
setwd(wd)
in_dir = paste(wd, '../results/drugs/', sep='')
virus_vs_control_down_file = paste(in_dir, "virus_vs_control_down.tsv", sep='')
virus_vs_control_up_file = paste(in_dir, "virus_vs_control_up.tsv", sep='')
rofa_virus_vs_control_down_file = paste(in_dir, "rofa_and_virus_vs_virus_down.tsv", sep='')
rofa_virus_vs_control_up_file = paste(in_dir, "rofa_and_virus_vs_virus_up.tsv", sep='')

to_process = c(virus_vs_control_down_file, virus_vs_control_up_file,
  rofa_virus_vs_control_up_file, rofa_virus_vs_control_down_file)

read_file = function(in_file) {
  df = read.csv(in_file, sep="\t", header=FALSE)
  colnames(df) = c("short_name", "accession", "name", "length", "compound", "cas")
  return(df)
}

multi_drug_hits = function(df) {
  columns_to_compare = c("short_name", "name", "accession", "compound")
  unique_rows = unique(df[, columns_to_compare])
  drug_hits = table(unique_rows$compound)
  drugs_against_multiple = names(drug_hits[drug_hits > 1])
  return(subset(unique_rows, compound %in% drugs_against_multiple))
}

table_of_multi_hits = function(df) {
  tab = data.frame(table(factor(multi_drug_hits(df)$compound)))
  colnames(tab) = c("compound", "hits")
  return(tab)
}

write_tables = function(df, in_file) {
  out_multidrug = paste(file_path_sans_ext(in_file), "_multidrug_hits.tsv", sep='')
  out_multidrug_counts = paste(file_path_sans_ext(in_file), "_multidrug_counts.tsv", sep='')
  write.table(multi_drug_hits(df), out_multidrug, quote=FALSE, row.names=FALSE, sep="\t")
  write.table(table_of_multi_hits(df), out_multidrug_counts, quote=FALSE, row.names=FALSE,
              sep="\t")
}

process_file = function(in_file) {
  df = read_file(in_file)
  write_tables(df, in_file)
}
