library(readr)
library(ggvis)
library(dplyr)

setwd("~/Desktop/")
alldat <- read_delim("s409_hmc_fixed.tsv", delim="\t")

#dat <- filter(alldat, chr!="chrM")
dat <- alldat

dat <- filter(dat, (CT_counts)>=10, (ox_CT_count)>=10)   #coverage filter
dat <- rowwise(dat) %>%  mutate(., fisher=fisher.test(x=matrix(c(C_counts, (CT_counts-C_counts), ox_C_count, (ox_CT_count-ox_C_count)),ncol=2), alternative="greater")$p.value)
dat <-   dplyr::select(dat, chr, pos, strand, C_counts, CT_counts, ox_C_count, ox_CT_count, fisher) 
dat <- mutate(dat, ratiobs=(C_counts)/(CT_counts), ratiooxbs=(ox_C_count)/(ox_CT_count)  )
dat <- mutate(dat, diff=ratiobs-ratiooxbs)
dat$padjust <- p.adjust(dat$fisher, method="fdr")

#randomly sample table for plotting
ggplot(sample_n(dat,1e6), aes(x=-log10(padjust), y=CT_counts))+geom_point(alpha=0.05)

