library(tidyverse)
library(limma)
library(ggbeeswarm)
library(ggpubr)

file_names <- list.files('gpr/', pattern='.gpr')
file_paths <- paste('gpr/', file_names, sep="")


# All columns from STARTSCRIPT.R
gpr_columns <- c("Block", "Column", "Row", "Name", "ID", "X", "Y", "Dia.",
                 "F635 Median", "F635 Mean", "F635 SD", "F635 CV", "B635",
                 "B635 Median", "B635 Mean", "B635 SD", "B635 CV",
                 "% > B635+1SD", "% > B635+2SD", "F635 % Sat.", "F532 Median",
                 "F532 Mean", "F532 SD", "F532 CV", "B532", "B532 Median",
                 "B532 Mean", "B532 SD", "B532 CV", "% > B532+1SD",
                 "% > B532+2SD", "F532 % Sat.", "Ratio of Medians (635/532)",
                 "Ratio of Means (635/532)", "Median of Ratios (635/532)",
                 "Mean of Ratios (635/532)", "Ratios SD (635/532)",
                 "Rgn Ratio (635/532)", "Rgn R2 (635/532)", "F Pixels",
                 "B Pixels", "Circularity", "Sum of Medians (635/532)",
                 "Sum of Means (635/532)", "Log Ratio (635/532)",
                 "F635 Median - B635", "F532 Median - B532", "F635 Mean - B635",
                 "F532 Mean - B532", "F635 Total Intensity",
                 "F532 Total Intensity", "SNR 635", "SNR 532", "Flags",
                 "Normalize", "Autoflag")

# My subset of columns
gpr_subset <- c("Name", "ID",
                "F635 Median", "F635 Mean", "F635 SD",
                "B635 Median", "B635 Mean", "B635 SD",
                "F532 Median", "F532 Mean", "F532 SD",
                "B532 Median", "B532 Mean", "B532 SD",
                "F Pixels", "B Pixels", "Circularity",
                "F635 Median - B635", "F532 Median - B532",
                "Flags")

# Easy variable names for subset of columns (snake case)
gpr_colnames <- c("name", "id",
                  "f635_median", "f635_mean", "f635_sd",
                  "b635_median", "b635_mean", "b635_sd",
                  "f532_median", "f532_mean", "f532_sd",
                  "b532_median", "b532_mean", "b532_sd",
                  "f_pixels", "b_pixels", "circularity",
                  "fb635_median", "fb532_median",
                  "flags")

## DATA IMPORT
gpr = tibble()
pools = paste("pool", 1:4, sep="_")
sets = paste("set", 1:2, sep="_")

for(i in 1:8) {
  slide <- read.maimages(file_paths[i], source="genepix.median",
                         other.columns=gpr_subset)$other
  slide <- as.data.frame(slide)
  colnames(slide) = gpr_colnames
  
  slide <- as_tibble(slide) %>% 
    mutate(pool = pools[ceiling(i/2)],
           set = sets[2-i%%2])
  
  gpr <- rbind(gpr, slide)
}

## FILTERING AND DB%
db <- read_tsv("data/fragment_rsd.txt") %>% 
  select(PrEST, `%DB`) %>% 
  rename(name = PrEST, db = `%DB`)

d_clean <- gpr %>% 
  # Incorporate %DB
  left_join(db, by = "name") %>% 
  # Filter rows that are flagged, non-HPPR, too small or too weak
  filter(flags == 0,
         str_starts(name, "HPRR"),
         f_pixels >= 30,
         f532_median - b532_median > 10 * b532_sd) %>%
  # If there are replicates, only keep the strongest one
  arrange(name, pool, desc(f532_median)) %>%
  distinct(name, pool, .keep_all = T)


## TRANSFORM AND DEFINE SIGNAL
calc_stats <- function(d){
  mfi <- d$fb635_median
  mn <- mean(mfi)
  sd <- sd(mfi)
  d_new <- d %>% 
    mutate(mfi = mfi,
           sds = (mfi-mn) / sd)
  return(d_new)}

# Use per glass stats
d_rank = tibble()
for(p in pools){
  for(s in sets){
    glass <- calc_stats(filter(d_clean, pool == p, set == s))
    d_rank <- rbind(d_rank, glass)}}
d_rank <- d_rank %>%
  select(name, db, pool, set, sds, everything()) %>%
  arrange(desc(sds))

## Incorporate gene info
d_42k <- as_tibble(read.delim("data/42k_array.txt")) %>% 
  rename(name = PrEST)

d_id <- d_rank %>% 
  left_join(d_42k, by = "name")
d_id <- d_id %>%
  arrange(pool, desc(sds)) %>% 
  mutate(prank = c(1:nrow(filter(d_id, pool == "pool_1")),
                   1:nrow(filter(d_id, pool == "pool_2")),
                   1:nrow(filter(d_id, pool == "pool_3")),
                   1:nrow(filter(d_id, pool == "pool_4")))) %>%
  # Remove rows with no Uniprot ID
  filter(Uniprot != "") %>% 
  arrange(desc(sds)) %>%
  # Keep all except antigen names
  select(name, db, pool, set, sds, prank, Gene, Gene.desc, Uniprot,
         everything(), -Ag.name)

## Write Uniprot IDs
# to_write <- filter(d_id, sds > 4)
# write_csv(as.data.frame(to_write$Uniprot), col_names = F,
#           file = paste("data/uniprot_id_",
#                        format(Sys.time(),"%y%m%d_%H%M%S"),
#                        ".txt", sep = ""))
# -> then manually obtain GO-list from Uniprot and save as uniprot_GO.tsv

## Incorporate downloaded GO and flag keywords
go <- read_tsv("data/uniprot_GO.tsv") %>% 
  rename(Uniprot = Entry)

# Create regexp filters based on keywords
kw_bioprocess <- paste(
  read_delim("keywords/keywords_bioprocess.txt", "\t", col_names=F)$X1,
  collapse = "|")
kw_cellcomp <- paste(
  read_delim("keywords/keywords_cellcomp.txt", "\t", col_names=F)$X1,
  collapse = "|")

d_go <- d_id %>% 
  left_join(go, by = "Uniprot") %>%
  arrange(desc(sds)) %>%
  # kw1: immunological function
  # kw2: extracellular/membrane localization
  mutate(kw1 = str_detect(`Gene ontology (biological process)`,kw_bioprocess),
         kw2 = str_detect(`Gene ontology (cellular component)`,kw_cellcomp)) %>%
  select(name, db, sds, prank, pool, kw1, kw2, Gene,
         Uniprot, `Protein names`, names(go),
         everything(), -Gene.desc)

# Subset all >4 sd
d_sub <- filter(d_go, sds>4)

# Subset kw
d_kw_and <- filter(d_sub, kw1 == TRUE & kw2 == TRUE)
d_kw_or <- filter(d_sub,
                  !(kw1 == TRUE & kw2 == TRUE) & (kw1 == TRUE | kw2 == TRUE))

## Write results
# write_csv(d_sub, file = paste("results",format(Sys.time(),"%y%m%d_%H%M%S"),
#                               ".txt", sep = ""))

###############################################################################

## Write de-duplicated results
d_dedup <- d_sub %>%
  add_count(name) %>%
  distinct(name, .keep_all = TRUE)
# write_tsv(d_dedup, file = paste("dedup_results",
#                                 format(Sys.time(),"%y%m%d_%H%M%S"),".txt",
#                                 sep = ""))

## Write gene list for antigen collection
# gene_list <- str_split(d_sub$Gene, ";")
# gene_vector <- unlist(gene_list)
# write(gene_vector, "42k_genes.txt")