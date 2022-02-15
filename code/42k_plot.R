## PLOTS
# Which data structure will be used for plotting?
###############
d_plot <- d_sub
###############

# Custom colorblind palette
mypal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
           "#D55E00", "#F0E442", "#56B4E9", "#999999")

# Index pools and count occurence
pool_index <- c(1:nrow(filter(d_plot, pool == "pool_1")),
                1:nrow(filter(d_plot, pool == "pool_2")),
                1:nrow(filter(d_plot, pool == "pool_3")),
                1:nrow(filter(d_plot, pool == "pool_4")))
(count_df <- count(d_plot, pool, set))

# Plot fall-off for each pool
ggplot(arrange(d_plot, pool, desc(sds))) +
  geom_point(aes(y = sds,
                 x = pool_index,
                 color = pool,
                 shape = set),
             size = 1.5) +
  geom_hline(yintercept = 4)

# Boxplot of pools and sets
ggplot(d_plot, aes(x = pool, y = sds)) +
  geom_boxplot(outlier.shape = NA, aes(fill = set)) +
  labs(y = "SDs", x = "Pool") +
  geom_quasirandom(dodge.width=0.8, size = 0.9, aes(group = set)) +
  geom_text(inherit.aes = F,
            stat="count",
            aes(x = pool,label=..count.., group = set, color = set),
            y = 2.5, 
            position = position_dodge(width = 0.8)) +
  theme_bw() +
  scale_x_discrete(labels = factor(1:4)) +
  scale_color_manual(values=mypal) +
  scale_fill_manual(values=mypal)

# Wilcox difference between sets for every pool
ggplot(d_plot, aes(x = set, y = sds, fill = set)) +
  geom_boxplot() +
  facet_wrap(~ pool, ncol = 2) +
  geom_quasirandom(dodge.width=0.8) +
  stat_compare_means(label.y = 75)

# Plot xy-correlation
ggplot(d_plot, aes(y=sds, x=f635_median, color = pool, shape = set)) +
  geom_point(size=2) +
  geom_line(aes(group = interaction(pool,set))) +
  labs(y = "SDs", x = "MFI", shape = "Set", color = "Pool") +
  theme_bw() +
  scale_shape_discrete(labels = factor(1:2)) +
  scale_color_manual(labels = factor(1:4),values=mypal)
