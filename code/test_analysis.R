library(tidyverse)

# Read well information
wells <- read_tsv("AK_wells.txt")

# Tidy count data
count <- read_tsv("bead_count.txt", col_names = F) %>% 
  mutate(well = 1:21) %>% 
  select(well, everything()) %>% 
  pivot_longer(!well, names_to = "id", values_to = "count") %>% 
  mutate(id = rep(1:384, times = 21))

# Tidy MFI data
mfi <- read_tsv("mfi.txt", col_names = F) %>% 
  mutate(well = 1:21) %>% 
  select(well, everything()) %>% 
  pivot_longer(!well, names_to = "id", values_to = "mfi") %>% 
  mutate(id = rep(1:384, times = 21))

# Join count, MFI and well data
d <- count %>% 
  mutate(mfi = mfi$mfi) %>% 
  left_join(wells, by = "well") %>% 
  # Add plate info
  arrange(id) %>% 
  mutate(plate = c(rep(1,8064/4), rep(2,8064/4),
                   rep(3,8064/4), rep(4,8064/4)))
  
bead_info <- read_tsv("bead_info.txt") %>% 
  mutate(id = as.integer(id))

################################################################################


################################## Check bead count
# Plot bead count for all id
plotdata = d
ggplot(data = plotdata) +
  geom_boxplot(aes(x = id, y = count, group = id, color = factor(plate)),
               outlier.shape = NA) +
  scale_color_discrete(name = "Plate") +
  geom_hline(yintercept = 30, color = "red") +
  theme_bw() +
  labs(y = "Bead count", x = "Bead ID")


# Summarize bead count per plate
summarise(group_by(plotdata,plate), median = median(count), mean = mean(count),
          min = min(count), max = max(count))

# Summarize bead count per id
summarise(group_by(plotdata,id), median = median(count), mean = mean(count),
          min = min(count), max = max(count)) %>% 
  arrange(median)

################################## Check coupling test

coupling_sum <- summarise(
  group_by(
    filter(d, test=="Coupling"),
    id, control),
  mean_MFI = mean(mfi)) 

coupling <- ungroup(coupling_sum) %>% 
  add_column(plate = c(rep(1,192), rep(2,192), rep(3,192), rep(4,192))) %>% 
  mutate(ctrl = if_else(id %in% c(378:379,381:384) & is.na(control),
                        id, as.integer(0)))

ggplot() +
  geom_col(data = filter(coupling, is.na(control)),
           aes(x = id, y = mean_MFI,
               fill = if_else(ctrl != 0, factor(ctrl), NULL)), width=1) +
  geom_col(data = filter(coupling, control == TRUE),
           aes(x = id, y = mean_MFI), fill = "black", width=1) +
  theme_bw() +
  labs(y = "mean(MFI)", x = "Bead ID", fill = "Antigen") +
  scale_fill_discrete(labels = c("IFNW1", "IFNA2A", "His6ABP", "EMPTY", 
                                 "anti-hIgG", "EBNA1", "PrEST"))

################################## Check sample test

samp_sum <- summarise(
  group_by(
    filter(d, test=="Sample"),
    id, control),
  mean_MFI = mean(mfi))

samp <- samp_sum %>% 
  mutate(ctrl = if_else(id %in% c(381:384) & is.na(control), id, as.integer(0)))

ggplot() +
  geom_col(data = filter(samp, is.na(control)),
           aes(x = id, y = mean_MFI,
               fill = if_else(ctrl != 0, factor(ctrl), NULL)), width=1) +
  geom_col(data = filter(samp, control == TRUE),
           aes(x = id, y = mean_MFI), fill = "black", width=1) +
  theme_bw() +
  labs(x = "Bead ID", y = "mean(MFI)", fill = "Antigen") +
  scale_fill_discrete(labels = 
                        c("His6ABP", "EMPTY", "anti-hIgG", "EBNA1", "Other"))


################################## Check antigen test

spec_control <- d %>% 
  filter(test == "Specific", control == TRUE) %>% 
  mutate(group = 1)
spec_control_rep <- spec_control
for(i in 2:6){
  spec_control <- d %>% 
    filter(test == "Specific", control == TRUE) %>% 
    mutate(group = i)
  spec_control_rep <- full_join(spec_control_rep, spec_control)
}

ag <- d %>% 
  filter(test == "Specific", is.na(control)) %>% 
  full_join(spec_control_rep) %>% 
  group_by(group, id, control) %>% 
  summarise(mean_MFI = mean(mfi)) %>% 
  ungroup() %>% 
  mutate(ctrl = if_else(id %in% c(381:384) & is.na(control),
                        id, as.integer(0))) %>% 
  left_join(distinct(select(d, group, target_id)), by = "group") %>% 
  left_join(bead_info, by = "id")

ggplot(data = filter(ag),
       aes(x = id, y = mean_MFI, fill = control)) +
  geom_col(position = position_stack(reverse = T),width=1) +
  # geom_col(data = filter(ag, control == TRUE),
  #          aes(x = id, y = mean_MFI), color = "black") +
  geom_text(aes(x = id, y = mean_MFI, 
                label = if_else(mean_MFI > 2000 & is.na(control),
                                bead_info, "")),
            size = 3, hjust = 0.8) +
  facet_wrap(~target_id, ncol = 2, scales = "free") +
  theme_bw() +
  labs(x = "Bead ID", y = "mean(MFI)", fill = "Sample (stacked)") +
  scale_fill_discrete(labels = c("Negative control", "Test duplicate"))

