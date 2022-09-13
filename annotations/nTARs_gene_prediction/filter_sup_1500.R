library(tidyverse)

nTARs_no <- read.table("/shared/ifbstor1/projects/minomics/pgrognet/RNAseq_Podo/06-SCT_results/annotations/nTARs_no_merged.bed")

nTARs_no_above_1500 <- nTARs_no %>% 
  mutate(size = nTARs_no$V3 - nTARs_no$V2) %>%
  filter(size >= 1500) %>%
  select(V1,V2,V3)

write.table(x = nTARs_no_above_1500,
            file = "/shared/ifbstor1/projects/minomics/pgrognet/RNAseq_Podo/06-SCT_results/annotations/nTARs_gene_prediction/nTARs_1500.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)
