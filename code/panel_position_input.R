library(tidyverse)

## Prerequisite step
# LIMS: 1) Input list of genes from 42k results and literature studies
#       2) Export df containing HPRR and HPRA, save as lims_genes.xls

# Here we merge those gene lists, remove rows without HPRA and de-duplicate HPRR
lims_42k_hpra <- read_tsv("lims_42k_genes.xls") %>% 
  mutate(origin = "42k_genes")
lims_lit_genes <- read_tsv("lims_lit_genes.xls") %>% 
  mutate(origin = "lit_genes")

hprr2collect <- rbind(lims_42k_hpra, lims_lit_genes) %>%
  select(-X4) %>% 
  filter(`Ag name` != "") %>% 
  distinct(PrEST, .keep_all = T)


# This goes to position-finding script
write_tsv(hprr2collect, "AK_covid.tsv")
