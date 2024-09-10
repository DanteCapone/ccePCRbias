# Script fro comparison of amplification efficiency and total relative read abundances

# Load necessary libraries
library(dplyr)
library(stringr)
library(fido)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(here)
library(gridExtra)
library(ggpubr)

# Load amplification efficiency and fido otu data

#18s
all_amp_effs_18s = read.csv(here("data/amp_effs/all_amp_effs_18s_all_sub.csv"))%>%
  # Extract the family name from Lambda.coord
  mutate(Family = str_extract(Lambda.coord, "[^_]+$"))
fido_18s_s1 = read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_18s_s2 = read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_18s_s3 = read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)




# Function to sum columns in fido datasets containing "S1", "S2", or "S3" in their names
sum_s_columns_18s <- function(df, family, suffix) {
  matching_columns <- names(df)[str_detect(names(df), suffix)]
  if (length(matching_columns) > 0) {
    normalized_sum <- sum(sapply(matching_columns, function(col) {
      column_sum_family <- sum(df[df$Family == family, col], na.rm = TRUE)
      column_sum <- sum(df[[col]], na.rm = TRUE)
      if (column_sum_family == 0) {
        return(0)
      }
      column_sum_family / column_sum
    }), na.rm = TRUE)
    return(normalized_sum)
  } else {
    return(NA)
  }
}



# 18s ---------------------------------------------------------------------


# Create the final dataframe for S1, S2, and S3
result_s1 <- all_amp_effs_18s %>%
  filter(str_detect(pool, "S1")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_18s(fido_18s_s1, Family, "S1"),
         SizeFraction = "S1") %>%
  select(Family, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

result_s2 <- all_amp_effs_18s %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_18s(fido_18s_s2, Family, "S2"),
         SizeFraction = "S2") %>%
  select(Family, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

result_s3 <- all_amp_effs_18s %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_18s(fido_18s_s3, Family, "S3"),
         SizeFraction = "S3") %>%
  select(Family, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

# Combine results into one dataframe
result_combined <- bind_rows(result_s1, result_s2, result_s3)

# Summarize the data and order Family by increasing Lambda.mean
family_order <- result_combined %>%
  group_by(Family) %>%
  summarise(mean_Lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_Lambda) %>%
  pull(Family)

summarized_result <- result_combined %>%
  group_by(Family, SizeFraction) %>%
  summarise(Lambda.mean = mean(Lambda.mean, na.rm = TRUE), Sum = mean(Sum, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Family = factor(Family, levels = family_order))

pearson_18s <- cor.test(~ Lambda.mean + Sum, data = summarized_result)

cat("18S Dataset:\n")
cat("Pearson's R:", pearson_18s$estimate, "\n")
cat("p-value:", pearson_18s$p.value, "\n\n")

# Plot with custom colors
custom_colors <- c(
  "Calanidae" = "#E41A1C",
  "Clausocalanidae" = "#377EB8",
  "Eucalanidae" = "#4DAF4A",
  "Euphausiidae" = "#984EA3",
  "Metridinidae" = "#FF7F00",
  "Oithonidae" = "#FFFF33",
  "other" = "#A65628",
  "Paracalanidae" = "#F781BF",
  "Rhincalanidae" = "#999999",
  "Salpidae" = "#66C2A5",
  "unidentified Calanoida" = "#FC8D62"
)

summarized_result %>% 
  filter(Sum > 0) %>%
  ggplot(., aes(x = Lambda.mean, y = Sum)) +
  geom_point(aes(color = Family, shape = SizeFraction), size = 10) +
  theme_classic() +
  labs(title = "Sum of Raw Relative Reads vs.\nAmplification Efficiency by Family (18S)",
       x = "Amplification Efficiency",
       y = "Sum of Relative Reads",
       color = "Family",
       shape = "Size Fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "right") +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(16, 17, 18), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")) -> amp_effs_vs_rra_18s
amp_effs_vs_rra_18s

#PNG & PDF Save
ggsave(
  filename = here("plots/methods_comparison/amp_effs_vs_rra_18s.png"),
  plot = amp_effs_vs_rra_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/amp_effs_vs_rra_18s.pdf"),
  plot = amp_effs_vs_rra_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)




# COI ---------------------------------------------------------------------

#COI
all_amp_effs_coi = read.csv(here("data/amp_effs/all_amp_effs_coi_all_sub.csv"))%>%
  # Extract the genus name from Lambda.coord
  mutate(Genus = str_extract(Lambda.coord, "[^_]+$"))
fido_coi_s1 = read.csv(here("data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_coi_s2 = read.csv(here("data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_coi_s3 = read.csv(here("data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE, row.names = 1)

# Function to sum columns in fido datasets containing "S1", "S2", or "S3" in their names
sum_s_columns_coi <- function(df, genus, suffix) {
  matching_columns <- names(df)[str_detect(names(df), suffix)]
  if (length(matching_columns) > 0) {
    normalized_sum <- sum(sapply(matching_columns, function(col) {
      column_sum_genus <- sum(df[df$Genus == genus, col], na.rm = TRUE)
      column_sum <- sum(df[[col]], na.rm = TRUE)
      if (column_sum_genus == 0) {
        return(0)
      }
      column_sum_genus / column_sum
    }), na.rm = TRUE)
    return(normalized_sum)
  } else {
    return(NA)
  }
}

# Create the final dataframe for S1, S2, and S3
result_s1 <- all_amp_effs_coi %>%
  filter(str_detect(pool, "S1")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_coi(fido_coi_s1, Genus, "S1"),
         SizeFraction = "S1") %>%
  select(Genus, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

result_s2 <- all_amp_effs_coi %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_coi(fido_coi_s2, Genus, "S2"),
         SizeFraction = "S2") %>%
  select(Genus, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

result_s3 <- all_amp_effs_coi %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Sum = sum_s_columns_coi(fido_coi_s3, Genus, "S3"),
         SizeFraction = "S3") %>%
  select(Genus, pool, Lambda.mean, Sum, SizeFraction) %>%
  ungroup()

# Combine results into one dataframe
result_combined <- bind_rows(result_s1, result_s2, result_s3)


# Summarize the data and order Genus by increasing Lambda.mean
genus_order <- result_combined %>%
  group_by(Genus) %>%
  summarise(mean_Lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_Lambda) %>%
  pull(Genus)

summarized_result <- result_combined %>%
  group_by(Genus, SizeFraction) %>%
  summarise(Lambda.mean = mean(Lambda.mean, na.rm = TRUE), Sum = mean(Sum, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genus = factor(Genus, levels = genus_order))

pearson_coi <- cor.test(~ Lambda.mean + Sum, data = summarized_result)
# Print the Pearson's R and p-values
cat("COI Dataset:\n")
cat("Pearson's R:", pearson_coi$estimate, "\n")
cat("p-value:", pearson_coi$p.value, "\n\n")

# Plot 

custom_colors <- c(
  "Pyroteuthis" = "#F8766D",
  "Muggiaea" = "#D89000",
  "Rosacea" = "#A3A500",
  "Euchaeta" = "#39B600",
  "Metridia" = "#00BF7D",
  "Ditrichocorycaeus" = "#00BFC4",
  "Sagitta" = "#00B0F6",
  "Pleuromamma" = "#9590FF",
  "Acrocalanus" = "#E76BF3",
  "Candacia" = "#FF62BC",
  "Calanus" = "#FF8B00",
  "Eucalanus" = "#00C08B",
  "unidentified Sagittidae" = "#00BA38",
  "Ctenocalanus" = "#00B4F0",
  "Euphausia" = "#E5841B",
  "Clausocalanus" = "#D39200",
  "other" = "#6A3D9A",
  "unidentified Calanoida" = "#FF6103"
)
summarized_result %>% 
  filter(Sum > 0) %>%
  filter(Genus != "unidentified Sagittidae") %>% 
ggplot(., aes(x = Lambda.mean, y = Sum)) +
  geom_point(aes(color = Genus, shape = SizeFraction), size = 10) +
  labs(title = "Sum of Raw Relative Reads vs.\n Amplification Efficiency by Genus (COI)",
       x = "Amplification Efficiency",
       y = "Sum of Relative Reads",
       color = "Genus",
       shape = "Size Fraction") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "right") +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(16, 17, 18), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")) -> amp_effs_vs_rra_coi
amp_effs_vs_rra_coi

#PNG & PDF Save
ggsave(
  filename = here("plots/methods_comparison/amp_effs_vs_rra_coi.png"),
  plot = amp_effs_vs_rra_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/amp_effs_vs_rra_coi.pdf"),
  plot = amp_effs_vs_rra_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)


grid.arrange(amp_effs_vs_rra_coi
             , amp_effs_vs_rra_18s
             , ncol = 2 )




