library(illuminaHumanv4.db)
library(readr)
library(dplyr)
library(tidyr)
library(plyr)



lncRNAannots <- read_csv(file.path("~/Desktop","lncRNA Array Gene List.csv"))
names(lncRNAannots) <- gsub(" ", "", sub("\\#", "", names(lncRNAannots)))

# using BioC annotations

x <- illuminaHumanv4REFSEQ
refseqs <- mappedRkeys(x)

lncRNAannots <- select(lncRNAannots,GeneSymbol, Alias, Refseq, OfficialFullName, RT2CatalogNumber)
lncRNAannots_onarray <- filter(lncRNAannots, lncRNAannots$Refseq %in% refseqs)

probes <- mget(lncRNAannots_onarray$Refseq, revmap(x))

# check to see if these probes uniquely assay the lncRNA in question
temp <- lapply(probes, function(probelist) {
  unique(
    sapply(probelist, function(probe){
      get(probe, illuminaHumanv4REFSEQ)
    })
  )
})

probes <- ldply(probes, function(x) paste(x, collapse=","))
names(probes) <- c("Refseq", "IlluminID")

lncRNAannots <- full_join(lncRNAannots, probes, by="Refseq")

write.table(lncRNAannots, "lncRNA_Array_Gene_List_with_IlluminaIDs.xls", row.names=F, col.names=T, sep="\t", quote=F)


# using Illumina annotations
illmn_annots <- read_delim(file.path("~/Desktop", "HumanHT-12_V4_0_R2_15002873_B.txt"), delim="\t", skip=0)

lapply(lncRNAannots$Refseq, function(refseq) {
  grep("refseq", illmn_annots$RefSeq_ID)
})

