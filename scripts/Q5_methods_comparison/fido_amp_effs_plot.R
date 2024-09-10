#Load data
all_fits_18s=read.csv(here("data/amp_effs/all_amp_effs_18s_all_sub.csv"))
# Plot
all_fits_18s %>%
  ggplot(., aes(x = Lambda.mean, y = Lambda.coord)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(pool_type)), height = 0.2) +
  geom_point(aes(color = as.factor(pool_type)),size = 8, alpha=0.8) + 
  labs(title="18S Family Amplification Efficiencies",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Pool") +
  facet_wrap(~size_fraction, nrow=3)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 

#By taxa
amp_effs_all_by_taxa_18s <- all_fits_18s %>%
  filter(rank_pool != "5") %>% 
  ggplot(aes(x = as.factor(size_fraction), y = Lambda.mean)) + 
  geom_hline(aes(yintercept = 0), color = "black", alpha = 1, size = 2) +
  geom_point(data = . %>% filter(pool_type != "AllandSub"), 
             aes(color = as.factor(rank_pool)), size = 8, alpha = 0.9, stroke = 2) +
  geom_point(data = . %>% filter(pool_type == "AllandSub"),
             aes(), size = 4, shape = 18, fill = "black",stroke=1, alpha=0.9) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5, color = as.factor(rank_pool)), width = 0.2) +
  labs(title = "18S Family Amplification Efficiencies",
       x = "Size Fraction",
       y = "Centered Log Ratio",
       color = "Pool") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-0.25, 0.25, by = 0.05)) +
  scale_x_discrete(labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  facet_wrap(~Lambda.coord) +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 14)) +
  scale_color_manual(values = c("#FFABAB", "#c996d4", "#93d182", "#a8bbe3", "#000000"), 
                     labels = c("All Pools","Pool 1", "Pool 2", "Pool 3", "All+Subpools"))

amp_effs_all_by_taxa_18s
#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_18s.png"),
  plot = amp_effs_all_by_taxa_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa.pdf"),
  plot = amp_effs_all_by_taxa_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)



# COI ---------------------------------------------------------------------



# Plot
all_fits_coi %>%
  ggplot(., aes(x = Lambda.mean, y = Lambda.coord)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(rank_pool)), height = 0.2) +
  geom_point(aes(color = as.factor(rank_pool)),size = 8, alpha=0.8) + 
  labs(title="18S Family Amplification Efficiencies",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Pool") +
  facet_wrap(~size_fraction, nrow=3)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 

#By taxa
amp_effs_all_by_taxa_coi=all_fits_coi %>%
  filter(rank_pool != "5") %>% 
  ggplot(., aes(x = as.factor(size_fraction), y = Lambda.mean)) + 
  geom_hline(aes(yintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_point(aes(color = as.factor(rank_pool),), size = 8, alpha=1.4, stroke=2) + 
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5,color = as.factor(rank_pool)), width=0.2) +
  labs(title = "COI Genus Amplification Efficiencies",
       x = "Size Fraction",
       y = "Centered Log Ratio",
       color = "Pool") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-0.25, 0.25, by = 0.05)) +
  scale_x_discrete(labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"))+
  facet_wrap(~Lambda.coord) +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 16),
        axis.title.x =element_text(size = 16),
        axis.title.y =element_text(size = 16),
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 14))+  # Increase legend entries size+
  scale_color_manual(values=c("#FFABAB","#c996d4", "#93d182", "#a8bbe3"),labels=c("All","Pool 1","Pool 2","Pool 3", "All and Subpools"))
# Adjust font size for y-axis tick labels
amp_effs_all_by_taxa_coi
#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_coi.png"),
  plot = amp_effs_all_by_taxa_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_coi.pdf"),
  plot = amp_effs_all_by_taxa_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)
