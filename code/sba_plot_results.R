
# Custom colorblind palette
mypal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
           "#D55E00", "#F0E442", "#56B4E9", "#999999")

# Show control beads
d_plot <- d_full %>% 
  mutate(source_plate_f = factor(source_plate, levels = c("DSPATI_MIX_1",
                                                          "DSPATI_MIX_2",
                                                          "DSPATI_9",
                                                          "DSPATI_10",
                                                          "EIMS_14",
                                                          "EIMS_16",
                                                          "control",
                                                          "empty")))

ggplot(filter(d_plot, ag_name %in% c("anti-hIgG", "EBNA1", "EMPTY", "His6ABP")),
       aes(x=source_plate_f, y = mfi)) +
  geom_boxplot(aes(color = source_plate_f)) +
  facet_wrap(~ag_name, ncol=2) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "MFI", color = "Source plate") +
  scale_color_manual(values=mypal, labels = c("Phase 1, plate 1/2",
                                              "Phase 1, plate 2/2",
                                              "Phase 2",
                                              "Phase 3",
                                              "EIMS plate 14",
                                              "EIMS plate 16",
                                              "Commercial plasma",
                                              "Empty well"
                                              ))



## Per-sample normalization

# mfi
ggplot(filter(d_trimmed, bead_id %in% c(3, 5, 377, 378, 379, 381, 382)),
       aes(x = bead_info,
           y = mfi)) +
  geom_quasirandom(aes(color = simple_class)) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust=0, vjust=0)) +
  xlab("") +
  ylab("MFI") +
  labs(color = "Sample type") +
  theme_bw() +
  scale_color_manual(values=mypal, labels = c("Covid", "No Covid"))

# mads
ggplot(filter(d_trimmed, bead_id %in% c(3, 5, 377, 378, 379, 381, 382)),
       aes(x = bead_info,
           y = mads_both)) +
  geom_quasirandom(aes(color = simple_class)) +
  geom_hline(yintercept = 1.5) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust=0, vjust=0)) +
  xlab("") +
  ylab("MADs by antigen x sample") +
  labs(color = "Sample type") +
  theme_bw() +
  scale_color_manual(values=mypal, labels = c("Covid", "No Covid"))


## Called reactive
# How many reactive per sample?
reactive_per_sample <- d_called %>% 
  group_by(sba_id) %>% 
  count(reactive) %>% 
  filter(reactive == TRUE)

ggplot(reactive_per_sample) +
  geom_histogram(aes(x = n), bins = 14) +
  xlab("n of reactivities per sample") +
  theme_bw()

reactive_per_ag <- d_called %>% 
  group_by(bead_id) %>% 
  count(reactive) %>% 
  filter(reactive == TRUE)

ggplot(reactive_per_ag) +
  geom_histogram(aes(x = n), bins = 32) +
  xlab("n of reactivities per antigen") +
  theme_bw()


########## Significant ag ############################
beads_relevant <- c(378,1,234,281,361,274,33,10)
beads_various <- c(279,254,44,236,154,25,347,321,336,271,166,177,327,68,16)
beads_excluded <- c(369,72,240,61,192,131,342,304,137,189,101)


# Facet
tag_sig_pati <- d_called %>% 
  select(sample_name, phase, bead_id, reactive) %>% 
  pivot_wider(names_from = phase, names_prefix = "phase_",
              values_from = reactive) %>% 
  # If for the current antigen, the current patient is reactive in any phase
  mutate(sig_pat = if_else(phase_1|phase_2|phase_3, TRUE, FALSE)) %>% 
  select(sample_name, bead_id, sig_pat)

plot_facet <- d_called %>% 
  filter(bead_id %in% c(beads_relevant)) %>%
  left_join(tag_sig_pati, by = c("sample_name", "bead_id"))

ggplot(plot_facet,
       aes(x = factor(phase), y = mads_sample,
           color = factor(interaction(simple_class, reactive),
                          levels = c("covid.TRUE", "covid.FALSE",
                                     "no_covid.TRUE", "no_covid.FALSE")),
           group = simple_class)) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_line(data = filter(plot_facet, sig_pat), alpha = 0.5, aes(group = sample_name)) +
  labs(x = "", y = "MADs by sample") +
  facet_wrap(~reorder(bead_info, bead_id), ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13)) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = 
                     c("Phase 1", "Phase 2", "Phase 3", "Control group")) +
  scale_color_manual(values=mypal)


######### IFN #################################################################

# Facet
tag_sig_pati <- d_called %>% 
  select(sample_name, phase, bead_id, reactive) %>% 
  pivot_wider(names_from = phase, names_prefix = "phase_",
              values_from = reactive) %>% 
  # If for the current antigen, the current patient is reactive in any phase
  mutate(sig_pat = if_else(phase_1|phase_2|phase_3, TRUE, FALSE)) %>% 
  select(sample_name, bead_id, sig_pat)

t1_ifn <- antigen_info %>%
  filter(str_detect(gene, "IFN"),
          !str_detect(gene, "IFNG"))

d_ifn <- d_called %>% 
  filter(bead_id %in% c(32,193,328,379)) %>% 
  left_join(tag_sig_pati, by = c("sample_name", "bead_id"))

ggplot(d_ifn,
       aes(x = factor(phase), y = mads_sample,
           color = factor(interaction(simple_class, reactive),
                          levels = c("covid.TRUE", "covid.FALSE",
                                     "no_covid.TRUE", "no_covid.FALSE")),
           group = simple_class)) +
  geom_quasirandom(dodge.width = 0.8) +
  geom_line(data = filter(d_ifn, sig_pat), alpha = 0.5, aes(group = sample_name)) +
  labs(x = "", y = "MADs by sample") +
  facet_wrap(~reorder(bead_info, bead_id), ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13)) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Phase 1", "Phase 2", "Phase 3", "Control group"))  +
  scale_color_manual(values=mypal)

