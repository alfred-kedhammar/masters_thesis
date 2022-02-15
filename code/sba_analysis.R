library(tidyverse)
library(ggbeeswarm)

########### Position information
## Convert 384 to 96 plate
# Greiner384 read up->down, left->right with plates 1-4 in quadrants 2,3,1,4
greiner2pcr_plate <- rep(c(rep(c(1,2),8), rep(c(3,4),8)),12)

greiner2pcr_well_seq <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8),2)
greiner2pcr_well <- c()
for(i in 0:11){
  greiner2pcr_well <- c(greiner2pcr_well, greiner2pcr_well_seq+i*8)
}

greiner2pcr_id <- (greiner2pcr_plate-1)*96+greiner2pcr_well

## Source positions

source2sba <- read_tsv("source2sba.txt")

########### Sample information
## eims

eims_14 <- read_tsv("p14_info.txt")
eims_16 <- read_tsv("p16_info.txt") %>% 
  mutate(source_plate = "EIMS_16")

eims <- eims_14 %>% 
  add_row(eims_16) %>% 
  add_column(phase = NA)

## covid

dspati_info <- read_tsv("dspati_info.txt") %>% 
  rename("sample_name" = sample_id)

dspati <- read_tsv("dspati.txt") %>% 
  left_join(dspati_info, by = "sample_name") %>% 
  mutate(sample_name = str_pad(sample_name, 3, side = "left", pad = "0"))
  

# Merge covid, eims and positions
all_samples <- full_join(dspati,eims)

sba_id_info <- tibble("sba_id" = greiner2pcr_id,
                      "sba_plate" = greiner2pcr_plate,
                      "sba_pos" = greiner2pcr_well,
                      "greiner_pos" = 1:384) %>%
  arrange(sba_id) %>% 
  left_join(source2sba, by = "sba_id") %>% 
  left_join(all_samples, by = c("source_plate", "source_pos")) %>% 
  mutate(class = if_else(source_plate == "control", "control",
                         if_else(source_plate == "empty", "empty", class)),
         sample_name = if_else(source_plate %in% c("control", "empty"),
                               source_plate, sample_name),
         simple_class = if_else(class == "covid",
                                "covid",
                                if_else(class == "MS" | class == "HC",
                                        "no_covid",
                                        class)),
         sample_id = paste(sample_name, phase, sep ="_"))


# Import Luminex data and transform layout

mfi <- read_tsv("mfi.txt", col_names = FALSE) %>%
  add_column(greiner_pos = 1:384) %>% 
  pivot_longer(X1:X384, names_to = "bead_id", values_to = "mfi") %>% 
  mutate(bead_id = rep(1:384, 384))

count <- read_tsv("count.txt", col_names = FALSE) %>% 
  add_column(greiner_pos = 1:384) %>% 
  pivot_longer(X1:X384, names_to = "bead_id", values_to = "count") %>% 
  mutate(bead_id = rep(1:384, 384))

# Merge mfi and count
luminex <- count %>% 
  add_column(mfi = mfi$mfi)


######### Get antigen info from different sources

# from PrEST to uniprot
d_42k <- as_tibble(read.delim("42k_array.txt")) %>% 
  rename("uniprot" = Uniprot) %>% 
  select(PrEST, uniprot)
# export and get uniprot info

# from uniprot to GO
go <- read_tsv("go.tsv") %>% 
  select(-1)
names(go) <- c("uniprot", "protein", "go_process", "go_cellcomp")

# keywords
kw_bioprocess <- paste(
  read_delim("keywords/keywords_bioprocess.txt", "\t", col_names=F)$X1,
  collapse = "|")
kw_cellcomp <- paste(
  read_delim("keywords/keywords_cellcomp.txt", "\t", col_names=F)$X1,
  collapse = "|")

antigen_info <- read_tsv("AK_SBA.txt") %>% 
  left_join(d_42k, by = "PrEST") %>% 
  left_join(go, by = "uniprot") %>% 
  mutate(bead_info = if_else(is.na(gene),
                             paste(bead_id, ag_name, sep="_"),
                             paste(bead_id, gene, sep="_")),
         bead_info = if_else(str_length(bead_info) > 20,
                                        paste(str_sub(bead_info, 1, 8), "..."),
                                        bead_info),
         kw1 = str_detect(go_process, kw_bioprocess),
         kw2 = str_detect(go_cellcomp, kw_cellcomp))



# merge information
d_full <- luminex %>% 
  left_join(sba_id_info, by = "greiner_pos") %>% 
  left_join(antigen_info, by = "bead_id")

######### Normalize by sample
# MAD (Median Absolute Deviation) = Median of absolute deviations from median

# MAD assumption does not apply for positive control beads
d_calc_stats <- d_full %>% 
  filter(bead_id %in% 1:382,
         count > 30)

sample_metrics <-  d_calc_stats%>% 
  group_by(sba_id) %>% 
  summarise(mfi_sample_median = median(mfi),
            mfi_sample_mad = mad(mfi, constant = 1))

ag_metrics <- d_calc_stats %>% 
  group_by(bead_id) %>% 
  summarise(mfi_ag_median = median(mfi),
            mfi_ag_mad = mad(mfi, constant = 1))

d_norm <- d_calc_stats %>% 
  left_join(sample_metrics, by = "sba_id") %>% 
  mutate(mfi_sub_sample = mfi - mfi_sample_median,
         mads_sample = mfi_sub_sample / mfi_sample_mad) %>% 
  left_join(ag_metrics, by = "bead_id") %>% 
  mutate(mfi_sub_ag = mfi - mfi_ag_median,
         mads_ag = mfi_sub_ag / mfi_ag_mad) %>% 
  mutate(mfi_sub_both = mfi - mfi_sample_median - mfi_ag_median,
         mads_both = mfi_sub_both / (mfi_sample_mad * mfi_ag_mad))

######### Trimming
## Exclusion criteria
# bead count < 30 (facility standard is remove <16, flag <30)
# set aside PCR-negative patients
# one EIMS sample has NA class and is therefore excluded via simple_class

d_trimmed <- d_norm %>% 
  filter(!(sample_name %in% c("009","059","063","085")),
         simple_class %in% c("covid", "no_covid"))

d_called <- d_trimmed %>% 
  mutate(reactive = mads_both > 1.5)

n_reac_per_sample <- d_called %>% group_by(sba_id) %>% count(reactive) %>% 
  pivot_wider(names_from = reactive, values_from = n) %>% 
  rename(n = `TRUE`, N = `FALSE`) %>% 
  select(sba_id, n, N)
sba_id_info_reac <- sba_id_info %>% 
  left_join(n_reac_per_sample, by = "sba_id")


######### Comparison full phase to control ###############

# The data of each phase are compared separately to the eims controls
for(phase_counter in 1:3){
  d_comp_current_phase <- d_called %>% 
    filter(phase == phase_counter | is.na(phase))
    
  cont <- count(d_comp_current_phase, bead_id, simple_class, reactive) %>% 
    mutate(reactive = if_else(as.character(reactive) == "TRUE", "n", "N")) %>% 
    pivot_wider(names_from = reactive, values_from = n, values_fill = 0) %>% 
    select(bead_id, simple_class, n, N) %>% 
    arrange(bead_id, simple_class)
  
  # Fisher
  f_val <- c()
  for(id in sort(unique(d_comp_current_phase$bead_id))){
    cont_table <- as.matrix(filter(cont, bead_id == id)[,3:4])
    dimnames(cont_table) = list(c("covid", "no_covid"), c("n", "N"))
    
    ftest <- fisher.test(x = cont_table, alternative = "greater")
    f_val <- c(f_val, ftest$p.value)
  }
  
  cont_wide <- cont %>% 
    pivot_wider(names_from = simple_class,
                values_from = c(n, N),
                names_sep = "_") %>%
    mutate(larger_nN_in =
             if_else( (n_covid / (n_covid + N_covid)) > 
                        (n_no_covid / (n_no_covid + N_no_covid)),
                      "covid", "no_covid"),
           f = f_val,
           #q = q_val,
           phase_comp = phase_counter) %>% 
    mutate(n_shorthand =
             paste(n_covid, "/", N_covid, " ",
                   n_no_covid, "/", N_no_covid, sep ="")) %>% 
    select(bead_id, phase_comp, f, #q,
           n_covid, N_covid, n_no_covid, N_no_covid, n_shorthand, larger_nN_in)
  
  if(phase_counter == 1){
    cont_all <- cont_wide
  }else{
    cont_all <- full_join(cont_all, cont_wide)
  }
}
  
antigen_info_f <- antigen_info %>% 
  full_join(cont_all, by = "bead_id") %>% 
  select(names(cont_all), db, bead_info, protein, kw1, kw2, go_process,
         everything()) %>% 
  arrange(f)

d_comp <- d_called %>%
  left_join(cont_all, by = "bead_id") %>%
  select(sba_id, sample_id, mfi, mads_sample, mads_ag, mads_both,
         names(cont_all), db, bead_info, protein, kw1, kw2,
         go_process, everything())

sig_phase <- filter(antigen_info_f, f < 0.05) %>% 
  select(bead_id, gene, kw1, kw2, protein, go_process,
         phase_comp, f, n_shorthand) %>% 
  arrange(phase_comp) %>% 
  pivot_wider(names_from = phase_comp,
              values_from = c(f,n_shorthand), names_prefix = "phase_") %>% 
  select(bead_id, gene,
         f_phase_1, n_shorthand_phase_1,
         f_phase_2, n_shorthand_phase_2,
         f_phase_3, n_shorthand_phase_3,
         everything()) %>% 
  arrange(is.na(f_phase_1), is.na(f_phase_2), is.na(f_phase_3),
          f_phase_1, f_phase_2, f_phase_3)

sig_beads <- unique(filter(antigen_info_f, f < 0.05)$bead_id)
sig_phase_full <- filter(antigen_info_f, bead_id %in% sig_beads) %>% 
  select(bead_id, gene, kw1, kw2, protein, uniprot, go_process, 
         phase_comp, f, n_shorthand) %>% 
  arrange(phase_comp) %>% 
  pivot_wider(names_from = phase_comp, 
              values_from = c(f,n_shorthand), names_prefix = "phase_") %>% 
  select(bead_id, gene,
         f_phase_1, n_shorthand_phase_1,
         f_phase_2, n_shorthand_phase_2,
         f_phase_3, n_shorthand_phase_3,
         everything()) %>% 
  arrange(is.na(f_phase_1), is.na(f_phase_2), is.na(f_phase_3),
          f_phase_1, f_phase_2, f_phase_3)

# write_tsv(sig_phase_full,
#         paste("results_",format(Sys.time(),"%y%m%d_%H%M%S"),".txt",sep = ""))


########## Numbering

x <- d_called %>% 
  filter(phase == 3,
         reactive == TRUE)

# 85080 obserations in covid patients
#   of which 2874 (3.4%) were reactive

# 43281 obs in phase 1
#   of which 1341 (3.1%) were reactive

# 22397 obs in phase 2
#   of which 828 (3.7%) were reactive

# 19402
#   of which 705 were reactive

######### Longitudinal macro

x <- d_called %>% 
  filter(!is.na(phase)) %>% 
  select(sample_name, phase, bead_id, reactive, mads_sample) %>% 
  pivot_wider(names_from = phase, values_from = c(reactive,mads_sample))

x1 <- x %>% 
  filter(!reactive_1, !is.na(reactive_2)) %>% 
  mutate(increased = reactive_2)

count(x1, increased)
