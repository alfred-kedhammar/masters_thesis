library(tidyverse)

# Import hprr2gene
hprr2gene <- as_tibble(read.delim("../42k/data/42k_array.txt")) %>% 
  select(PrEST, Gene, Gene.desc)

# Make position df with gene info
plock <- read_tsv("AK_covid_pos.txt") %>%
  # Remove mostly empty columns occupying namespace
  select(-Gene, -PrEST.1) %>% 
  left_join(hprr2gene, by = "PrEST") %>% 
  select(PrEST, Ag.name, Gene, Gene.desc, everything())

# Make position df without gene info
hprr2pos <- plock %>% 
  select(-Gene, -Gene.desc)

pos_info <- names(hprr2pos)[2:length(hprr2pos)]


# 1) 42k hpra
d_42k <- read_tsv("42k_dedup.txt") %>% 
  # Rename "name" -> "PrEST
  rename(PrEST = name)
# This df has HPRA as unique key, meaning all other columns may be copied for some rows
d_42k_pos <- d_42k %>% 
  inner_join(hprr2pos, by = "PrEST") %>% 
  distinct(Ag.name, .keep_all = T) %>% 
  mutate(prio = "1_42k") %>% 
  select(prio, PrEST, db, sds, Gene, pos_info)

# 2) Literature genes
d_lit <- read_tsv("lims_lit_genes.xls") %>% 
  filter(`Ag name` != "") %>% 
  select(-X4, -`Ag name`)
d_lit_pos <- d_lit %>% 
  inner_join(hprr2pos, by = "PrEST") %>% 
  mutate(sds = NA, db = NA) %>%
  mutate(prio = "1_lit") %>% 
  select(prio, PrEST, db, sds, Gene, pos_info)

# 3) Fill upp remainder of panel on HPRR that:
# - were not reactive themselves
# - represents genes that showed other reactive HPRR
# - represents genes that were immunologically relevant and extracellularly accessible
# - are sorted in descending signal strength
bonus_genes <- unique(filter(d_42k, kw1 == T & kw2 == T)$Gene)
bonus_genes_sds <- d_42k %>% 
  filter(kw1 == T & kw2 == T) %>% 
  select(Gene, sds) %>% 
  distinct(Gene, .keep_all=T)
bl <- str_split(plock$Gene,";") %in% bonus_genes
d_bonus_pos <- plock[bl,] %>% 
  left_join(bonus_genes_sds, by = "Gene") %>%
  mutate(db = NA) %>% 
  mutate(prio = "2_bonus") %>% 
  select(prio, PrEST, db, sds, Gene, pos_info) %>% 
  arrange(desc(sds))

final <- rbind(d_42k_pos, d_lit_pos, d_bonus_pos) %>% 
  arrange(prio,
          ProjBox.4,
          ProjBox.3,
          ProjBox.2,
          ProjBox.1)


## WRITE RESULTS
# write_tsv(d_lit_pos, file = paste("plock_lit_",
# format(Sys.time(),"%y%m%d_%H%M%S"),".xls", sep = ""), na = "")
