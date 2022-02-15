library(tidyverse)

df <- read_tsv("MSC_6.tsv")

df_all <- df %>% 
  select("sample_name", "class", "gender", "Plate96", "Well96", "Pos96", 
         "Date of birth", "Colection Date") %>% 
  rename(date_sampled = `Colection Date`) %>% 
  mutate(date_born = as.Date(paste("19", `Date of birth`, sep=""), 
                             format = "%Y%m%d"),
         sample_age = lubridate::time_length(difftime(date_sampled, date_born), 
                                             "years"),
         half_plate = if_else(Pos96 <= 96/2,
                              paste(Plate96,1,sep="_"),
                              paste(Plate96,2,sep="_"))) %>% 
  arrange(Plate96, Pos96)

by_plate <- group_by(df_all, Plate96)
by_half_plate <- group_by(df_all, half_plate)

plate <- summarise(by_plate,
                    n_HC = sum(class == "Control", na.rm=T),
                    n_MS = sum(class == "Case", na.rm=T),
                    prop_HC = n_HC/(n_HC+n_MS),
                    n_M = sum(gender == "M", na.rm=T),
                    n_F = sum(gender == "F", na.rm=T),
                    prop_M = n_M/(n_M+n_F),
                    n_HC_M = sum(gender == "M" & class == "Control", na.rm=T),
                    n_HC_F = sum(gender == "F" & class == "Control", na.rm=T),
                    prop_HC_M = n_HC_M/(n_M+n_F),
                    med_age = median(sample_age, na.rm=T))
half_plate <- summarise(by_half_plate,
                      n_HC = sum(class == "Control", na.rm=T),
                      n_MS = sum(class == "Case", na.rm=T),
                      prop_HC = n_HC/(n_HC+n_MS),
                      n_M = sum(gender == "M", na.rm=T),
                      n_F = sum(gender == "F", na.rm=T),
                      prop_M = n_M/(n_M+n_F),
                      n_HC_M = sum(gender == "M" & class == "Control", na.rm=T),
                      n_HC_F = sum(gender == "F" & class == "Control", na.rm=T),
                      prop_HC_M = n_HC_M/(n_M+n_F),
                      med_age = median(sample_age, na.rm=T))

ggplot(plate, mapping = aes(n_HC, n_M)) +
  geom_point() +
  geom_text(aes(label=Plate96), hjust=-0.5, vjust=-0.5)

ggplot(half_plate, mapping = aes(n_HC, n_M)) +
  geom_point() +
  geom_text(aes(label=half_plate), hjust=-0.5, vjust=-0.5)

# Whole of plate 14 and second half of plate 16

d_14 <- df_all %>% 
  filter(Plate96 == 14, Pos96 < 93) %>% 
  rename(age = sample_age, source_pos = Pos96) %>% 
  add_column(source_plate = "EIMS_14") %>% 
  mutate(class = if_else(class == "Case", "MS", "HC")) %>% 
  select(sample_name, class, gender, age, source_plate, source_pos) %>% 
  arrange(source_pos)

d_16 <- df_all %>% 
  filter(Plate96 == 16, Pos96 > 49) %>% 
  rename(age = sample_age, source_pos = Pos96) %>% 
  add_column(source_plate = "EIMS_14") %>% 
  mutate(class = if_else(class == "Case", "MS", "HC")) %>% 
  select(sample_name, class, gender, age, source_plate, source_pos) %>% 
  arrange(source_pos)

# write_tsv(d_14, "p14_info.txt")
# write_tsv(d_16, "p16_info.txt")


