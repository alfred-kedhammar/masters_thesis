
########## Look at control beads across all samples

# Points
ggplot(filter(d_full, bead_id %in% 381:384), aes(x=sba_id, y = mfi)) +
  geom_point(aes(color = source_plate)) +
  facet_wrap(~as.factor(bead_id), ncol=2)

# Boxplot
ggplot(filter(d, ag_name %in% c("anti-hIgG", "EBNA1", "EMPTY", "His6ABP")), 
       aes(x=source_plate, y = mfi)) +
  geom_boxplot(aes(color = source_plate)) +
  facet_wrap(~ag_name, ncol=2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


########## General
## Look at specific bead across all samples
ggplot(filter(d_comp, ag_name == "IFNA2A"), aes(x=sba_id, y = mfi)) +
  geom_point(aes(color = simple_class))


## Look at specific sample across all beads
ggplot(filter(d_norm, sba_id %in% 149:151), aes(x=bead_id, y = mfi)) +
  geom_col(aes(color = class)) +
  geom_hline(aes(yintercept = mfi_sample_median)) +
  geom_hline(aes(yintercept = mfi_sample_median + mfi_sample_mad),
             color = "red") +
  facet_wrap(~sba_id)


########## Normalization evaluation

ggplot(filter(d_norm, sba_id %in% 149:151), aes(x=bead_id, y = mfi)) +
  geom_col(aes(color = class)) +
  geom_hline(aes(yintercept = mfi_sample_median)) +
  geom_hline(aes(yintercept = mfi_sample_median + mfi_sample_mad),
             color = "red") +
  facet_wrap(~sba_id)

ggplot(filter(d,
              sample_name %in% 17:32,
              bead_id %in% 1:380),
       aes(x=bead_id, y = mfi)) +
  geom_col(aes(color = class)) +
  facet_wrap(~sample_name, ncol=4)

# The antigen-normalized median correlates to stickyness
ggplot(filter(by_ag, bead_id %in% 1:380),
       aes(x = antigen_info$db[1:380], y = mfi_ag_mad)) +
  geom_point() +
  geom_text(aes(label=ifelse(mfi_ag_mad > 100,as.character(bead_id),'')),
            hjust=0,vjust=0)

# Normalizing by antigen decreases this tendency
ggplot(filter(d, bead_id %in% 1:380),
       aes(x = db, y = mads)) +
  geom_jitter(position = "jitter")
ggplot(filter(d, bead_id %in% 1:380),
       aes(x = db, y = madsmads)) +
  geom_point(position = "jitter")


####### Big plot
d_matrix <- filter(d, bead_id %in% 1:380, 
                   !(class %in% c("control", "empty"))) %>% 
  mutate(sample_id = if_else(is.na(phase), sample_name, 
                             paste(sample_name, phase, sep = "_"))) %>% 
  select(sample_id, ag_name, madsmads) %>% 
  pivot_wider(id_cols = c(sample_id, ag_name),
              names_from = ag_name,
              values_from = madsmads)
d_matrix_input <- as.matrix(select(d_matrix, -sample_id))
heatmap(d_matrix_input)

d_big <- d %>% 
  filter(bead_id %in% 350:380) %>%
  mutate(reactive = 
           mfi > 2000 |
           mads > 50)

ggplot(d_big %>% 
         mutate(sample_id = if_else(is.na(phase),
                                    sample_name,
                                    paste(sample_name, phase, sep = "_"))),
       aes(x = sample_id, y = ag_name)) +
  geom_tile(aes(fill = reactive)) +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.95, 
                                   vjust = 0.2),
        axis.text.y = element_text(size = 8))

########## Full length interferons
ggplot(filter(d, ag_name %in% c("IFNW1", "IFNA2A")), aes(x=sba_id, y = mads)) +
  geom_point(aes(color = class)) +
  facet_wrap(~ag_name, ncol=2)


########## Ordered fall-off of all reactivities
# mfi
ggplot(filter(d) %>% arrange(desc(mfi)), aes(x=1:nrow(d), y = mfi)) +
  geom_point(aes(color = bead_id == 383))


########## Antigen profile across samples + testing cut-off
d_plot <- d %>% 
  filter(bead_id %in% 1:46| bead_id %in% 381:384) %>%
  mutate(reactive = 
           mfi > 2000 |
           mads > 50)

# Variable grid Fall-off
ggplot(d_plot %>% 
         pivot_longer(cols = c(mfi,mads,madsmads), 
                      names_to = "var", values_to = "val") %>% 
         mutate(var_f = factor(var, levels = c("mfi", "mads", "madsmads"))) %>% 
         arrange(desc(val))) +
  geom_point(aes(x = 1:length(val), y = val, color = reactive)) +
  facet_wrap(bead_id ~ var_f, ncol = 6, scales = "free")

# Fall-off
ggplot(arrange(d_plot, desc(mfi))) +
  geom_violin(aes(x = 1:length(sba_id), y = madsmads, color = reactive)) +
  facet_wrap(~bead_id, ncol = 10, scales = "free_y")

# Boxplot
ggplot(d_plot) +
  geom_histogram(aes(y = madsmads_sd)) +
  geom_vline(xintercept = 2) +
  facet_wrap(~bead_id, ncol = 10)

# Samples
d_plot_samples <- d %>% 
  filter(sba_id %in% unique(sba_id)[1:16]) %>%
  mutate(reactive = if_else(mfi > 2000, "mfi > 2000",
                            if_else(mads > 50, "mads > 50", "None")))
ggplot(d_plot_samples) +
  geom_point(aes(x = bead_id, y = mfi, color = reactive)) +
  facet_wrap(~sba_id, ncol = 8, scales = "free_y")


########## Various interferons
ggplot(d_ifn %>% filter(class %in% c("covid", "MS", "HC"))) +
  geom_quasirandom(dodge.width=0.8, size = 0.9,
                   aes(y = madsmads, x = class, color = bead_info)) +
  facet_wrap(~bead_info, ncol = 5, scales = "free_y")


######### Empty wells
ggplot(filter(d, class == "empty"), aes(x=bead_id, y = mfi)) +
  geom_col(aes(color = class)) +
  geom_text(aes(label=ifelse(mads > 50,as.character(ag_name),'')),
            hjust=0,vjust=0) +
  facet_wrap(~sba_plate, ncol=2)


########## Longitudinal
# Based on diff and baseline
ggplot(filter(d_follow, diff > 4, baseline < 5)) +
  geom_line(aes(y = mads, x = phase,
                color = interaction(ag_name, sample_name))) +
  theme(legend.position = "none")

# Based on genes
ggplot(filter(d_follow, diff < -50, str_detect(gene, "IFN")),
       aes(y = mads, x = phase)) +
  geom_line(aes(color = interaction(ag_name, sample_name))) +
  theme(legend.position = "none") #+
  geom_text(aes(label=ifelse(diff< -1 & phase == 2, as.character(ag_name),'')),
            hjust=0,vjust=0)


ggplot(filter(d_follow, diff < -500), aes(y = mads, x = phase)) +
  geom_line(aes(color = interaction(ag_name, sample_name))) +
  theme(legend.position = "none") #+
geom_text(aes(label=ifelse(diff< -1 & phase == 2, as.character(ag_name),'')),
          hjust=0,vjust=0)

################################################################################
## Findings

# HIF3A reactivity from 42k analysis pool 2 can be traced to patient
ggplot(filter(d, PrEST == "HPRR3140817"), aes(x=sba_id, y = sds)) +
  geom_point(aes(color = class)) +
  facet_wrap(~ag_name, ncol=2)
