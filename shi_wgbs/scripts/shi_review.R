library(readr)
library(ggvis)
library(dplyr)

setwd("/Volumes/orch/group_dir/shi_wgbs/ox_bs/work/cpg_split/s409/")
dat <- read_delim("s409_hmc_fixed.tsv", delim="\t")
temp <- dat[1:500000]


temp %>% 
  filter(., (C_counts+CT_counts)>=10, (ox_C_count+ox_CT_count)>=10) %>%  #coverage filter
  rowwise() %>%  # operate on rows
  mutate(., fisher=fisher.test(x=matrix(c(C_counts, CT_counts, ox_C_count, ox_CT_count),ncol=2), alternative="greater")$p.value)  %>% # fisher's eact test
  select(., chr, pos, strand, C_counts, CT_counts, ox_C_count, ox_CT_count, fisher) %>%
  mutate(., logxform=-log10(fisher)) %>%
  compute_density(~logxform) %>% ggvis(~pred_, ~resp_) %>% layer_lines() 
