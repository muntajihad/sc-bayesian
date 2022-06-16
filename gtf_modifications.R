library(rtracklayer)
library(dplyr)
setwd("~/data")
#read gtf file and convert it to data frame
gtf <- as.data.frame(import("genes.gtf"))
# filter to choose only protein coding genes
gtf <- filter(gtf,grepl("NM_",transcript_id)) %>%
  group_by("gene"=gene_id, "chr"=seqnames) %>%
  summarise("start"=min(start),"end"=max(end)) %>%
  filter(!grepl("_",chr))

gtf <- data.frame(gtf) %>%
  mutate(length=end-start) %>%
  filter(length>2000)
  
promoter <- gtf %>%
  mutate(start = start -1000) %>%
  mutate(end=start + 2000) %>%
  mutate(length=end-start) %>%
  arrange(length)

gene.body <- gtf %>%
  mutate(start=start+1000) %>%
  mutate(length=end-start) %>%
  arrange(length)


