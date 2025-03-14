#Zooscan analysis final
librarian::shelf(tidyverse, googledrive, stringr,here,vegan,ggpubr, patchwork)
here()
source(here("PCR_bias_correction/scripts/Zooscan_analysis/zooscan_functions.R"))
source("PCR_bias_correction/scripts/helpful_functions/general_helper_functions.R")


#Load in the data

##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
  # Standardize key columns (some samples are incorerctly notated to compare with Zooscan)
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_")))%>%
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = max_size) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-2", # Assign 1-2
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  ) %>%
  select(-X,-Sample_ID_dot, -Sizefractionmm, -size_fraction_numeric) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1) 


# List all .tsv files in the folder
tsv_files <- list.files(here("PCR_bias_correction/data/Zooscan/"), pattern = "\\.tsv$", full.names = TRUE)

# Find the most recently added file
latest_ecotaxa <- tsv_files[which.max(file.info(tsv_files)$ctime)]
latest_ecotaxa

# Read the .tsv file into a data frame
zooscan_exp <- read.table(latest_ecotaxa, header=TRUE, sep="\t", encoding="latin1")


#Basic plots to show the relative abundance and biomass of zooscan results from each station
zooscan_processed=readEcotaxa(zooscan_exp)%>%
  #Cleaning up formatting issues
  mutate(sample_id = str_replace_all(sample_id, "-", "_"),
         sample_id = ifelse(sample_id == "c2_t1_h36", "ct2_t1_h36", sample_id),
         sample_id = ifelse(sample_id == "ct1_t8_h10", "c1_t8_h10", sample_id),
         sample_id = ifelse(sample_id == "ct2_t9_h19", "c2_t9_h19", sample_id),
         sample_id = ifelse(sample_id == "c3_bt6_h25", "c3_t6_h25", sample_id)) %>%
  #Need to fix hyperiids
  filter(!(object_annotation_category %in% c("Hyperiidea","part<Crustacea","dark<sphere", "multiple organisms","head<Chaetognatha",
                                             "egg<Actinopterygii","other<living","Insecta"))) %>%
  #Fix missing volume filtered 
  mutate(sample_tot_vol = case_when(
    sample_id == "ct1_t1_h28" ~ 279,
    sample_id == "ct1_t2_h29" ~ 278,
    TRUE ~ sample_tot_vol  # This line keeps the original values for all other rows
  ))


# Retain PC1 and Sample_ID_short, and calculate C-biomass, relative abundances, and biomass proportions
zooscan_combined <- zooscan_processed %>%
  # Calculate biomass
  transform_by_taxa_group(., "feret") %>%
  mutate(dryweight_C_mg = dryweight_C_ug / 1000) %>%
  # Add metadata first to get sample_conc, PC1, and Sample_ID_short
  left_join(metadata %>% select(sample_id, size_fraction, PC1, Sample_ID_short), 
            by = c("sample_id", "size_fraction")) %>%
  # Fill missing PC1 values within each Sample_ID_short
  group_by(Sample_ID_short) %>%
  fill(PC1, .direction = "downup") %>%
  ungroup() %>%
  # Compute relative abundances and biomass sums while preserving sample_conc
  group_by(sample_id,size_fraction_numeric ,size_fraction, object_annotation_category, PC1, Sample_ID_short) %>%
  summarise(count = n(),
            dryweight_C_mg_sum_taxa = sum(dryweight_C_mg, na.rm = TRUE),
            dryweight_C_ug_sum_taxa = sum(dryweight_C_ug, na.rm = TRUE),
            sample_conc = mean(sample_conc, na.rm = TRUE), 
            .groups = 'drop') %>%
  group_by(sample_id, size_fraction, size_fraction_numeric) %>%
  mutate(total = sum(count),
         relative_abundance = count / total,
         dryweight_C_mg_sum_sample = sum(dryweight_C_mg_sum_taxa, na.rm = TRUE),
         dryweight_C_ug_sum_sample = sum(dryweight_C_ug_sum_taxa, na.rm = TRUE)) %>%
  ungroup() %>%
  # Now calculate biomass/m2 and log biomass/m2
  mutate(log10_dryweight_C_ug_m2_taxa = log(dryweight_C_ug_sum_taxa),
         dryweight_C_mg_m2_taxa = dryweight_C_mg_sum_taxa * sample_conc,
         log10_dryweight_C_ug_m2_sample = log(dryweight_C_ug_sum_sample),
         dryweight_C_mg_m2_sample = dryweight_C_mg_sum_sample * sample_conc,
         biomass_prop_taxa = dryweight_C_mg_m2_taxa / dryweight_C_mg_m2_sample,
         biomass_taxa_clr = clr_convert(dryweight_C_mg_m2_taxa)) %>%
  group_by(sample_id) %>%
  mutate(total_check=sum(relative_abundance)/total) %>% 
  fill(PC1, .direction = "downup") %>% 
  fill(Sample_ID_short, .direction="downup") 

#Save
write.csv(zooscan_combined,here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))

# PC1 Labels
labels_for_PC1 <- metadata %>% 
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange(PC1)

# Identify top 10 most abundant categories for relative abundance plot
top_categories_abundance <- zooscan_combined %>%
  group_by(object_annotation_category) %>%
  summarise(total_abundance = sum(relative_abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(object_annotation_category)

# Identify top 10 most abundant categories for biomass proportion plot
top_categories_biomass <- zooscan_combined %>%
  group_by(object_annotation_category) %>%
  summarise(total_biomass = sum(biomass_prop_taxa, na.rm = TRUE)) %>%
  arrange(desc(total_biomass)) %>%
  slice_head(n = 10) %>%
  pull(object_annotation_category)

# Get unique taxa across all size fractions and plots
unique_taxa <- unique(c(top_categories_abundance, top_categories_biomass))

# Create an "Other" category for anything not in the top 10
zooscan_combined <- zooscan_combined %>%
  mutate(object_annotation_category = ifelse(object_annotation_category %in% unique_taxa,
                                             object_annotation_category,
                                             "Other"))

# Set "Other" as the last level for consistent ordering
zooscan_combined$object_annotation_category <- factor(
  zooscan_combined$object_annotation_category, 
  levels = c(unique_taxa, "Other")
)

# Filter for abundance and biomass data with the new "Other" category included
abundance_data <- zooscan_combined %>%
  group_by(sample_id, size_fraction, object_annotation_category, PC1) %>%
  summarise(relative_abundance = sum(relative_abundance), .groups = 'drop')

biomass_data <- zooscan_combined %>%
  group_by(sample_id, size_fraction, object_annotation_category, PC1) %>%
  summarise(biomass_prop_taxa = sum(biomass_prop_taxa), .groups = 'drop')

# Generate consistent color palette with grey for "Other"
exclude_colors <- c("#FFFF99", "#B3B3B3")  # Exclude yellow and grey from main palette
num_taxa <- length(unique_taxa)
color_palette <- brewer.pal(n = min(num_taxa + 1, 12), name = "Set3")
color_palette <- setdiff(color_palette, exclude_colors)
if (num_taxa > length(color_palette)) {
  color_palette <- colorRampPalette(color_palette)(num_taxa)
}
# Add grey for "Other"
color_palette <- c(color_palette, "Other" = "#B3B3B3")
taxa_colors <- setNames(color_palette, levels(zooscan_combined$object_annotation_category))

# Order size fractions explicitly
ordered_size_fractions <- c("0.2-0.5", "0.5-1", "1-2",">2")

# Function to create abundance and biomass plots for each size_fraction
create_plot <- function(data, y, ylabel, title_prefix) {
  lapply(seq_along(ordered_size_fractions), function(i) {
    size <- ordered_size_fractions[i]
    
    p <- ggplot(data %>% filter(size_fraction == size), 
                aes(x = as.factor(PC1), y = .data[[y]], fill = object_annotation_category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = taxa_colors, name="Taxa") + 
      scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
      labs(
        title = paste(title_prefix, size),
        x = expression(paste("Offshore ", PC1, " Onshore")),
        y = ylabel
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increased for better readability (O)
        axis.title.x = element_text(size = 16),  # Larger axis title (O)
        axis.title.y = element_text(size = 16),  # Larger y-axis title (O)
        plot.title = element_text(size = 18, face = "bold"),  # Increased title size (O)
        legend.text = element_text(size = 14),  # Larger legend text (L)
        legend.title = element_text(size = 16, face = "bold"),  # Larger, bold legend title (L)
        legend.position = ifelse(i %in% c(3), "bottom", "none")  # Show legend only on first plot
      )
    
    return(p)
  })
}

# Create plots
abundance_plots <- create_plot(abundance_data, "relative_abundance", "Relative Abundance", "")

abundance_plots
# Arrange in grid layout with left column for Abundance and right for Biomass
abundance_plot_save <- wrap_plots(
  c(abundance_plots), 
  ncol = 4
)

# Display combined plot
abundance_plot_save

# Define output directory
output_dir <- "PCR_bias_correction/figures/zooscan"

# Save as PNG
ggsave(
  filename = file.path(output_dir, "zooscan_abundance_stacked_bar.png"),
  plot = abundance_plot_save,  # Using the final combined plot
  dpi = 300,
  width = 18,  # Adjust for better layout
  height = 7,
  units = "in",
  limitsize = FALSE  # Ensures large plots are not cropped
  
)

# Save as PDF
ggsave(
  filename = file.path(output_dir, "zooscan_abundance_stacked_bar.pdf"),
  plot = abundance_plot_save,
  dpi = 300,
  width = 18,
  height = 7,
  units = "in",
  limitsize = FALSE  # Ensures large plots are not cropped
  
)


biomass_plots <- create_plot(biomass_data , "biomass_prop_taxa", "Biomass Proportion", "")

biomass_plots
# Arrange in grid layout with left column for Abundance and right for Biomass
biomass_plot_save <- wrap_plots(
  c(biomass_plots), 
  ncol = 4
)

# Display combined plot
biomass_plot_save



# Save as PNG
ggsave(
  filename = file.path(output_dir, "zooscan_biomass_stacked_bar.png"),
  plot = biomass_plot_save,  # Using the final combined plot
  dpi = 300,
  width = 18,  # Adjust for better layout
  height = 7,
  units = "in",
  limitsize = FALSE  # Ensures large plots are not cropped
  
)

# Save as PDF
ggsave(
  filename = file.path(output_dir, "zooscan_biomass_stacked_bar.pdf"),
  plot = biomass_plot_save,
  dpi = 300,
  width = 18,
  height = 7,
  units = "in",
  limitsize = FALSE  # Ensures large plots are not cropped
  
)

# Calanoids ---------------------------------------------------------------

#CALANOID DATA FRAME FOR PLOTTING
zoop_calanoid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Calanoida")

#Calanoid dataframe for relative abundances
zoop_calanoid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Calanoida")

#Save Biomass and relative abundances
write.csv(zoop_calanoid_by_sample,here("PCR_bias_correction/data/Zooscan/zoop_calanoid_by_sample_biomass.csv"))
write.csv(relative_abundances_map,here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))



#Zooscan and Calanoid Biomass by cycle
zoop_calanoid_by_sample %>%
  group_by(sample_id, size_fraction, cycle) %>%  # Ensure cycle is included in grouping
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%  # Summarize all numeric columns
  ggplot(aes(x = cycle, y = dryweight_C_mg_m2_sample, color=size_fraction)) +
  geom_violin() +
  geom_point() +
  facet_wrap(~size_fraction, nrow=3) +
  labs(
    title = "Dry Weight by Cycle and Size Fraction",
    x = "Cycle",
    y = "Dry Weight (mg/m²)"
  ) +
  theme_minimal()


zoop_calanoid_by_sample %>%
  ggplot(aes(x = as.factor(PC1), y = biomass_prop_taxa)) +
  geom_bar(stat = "identity", fill = "#4F86F7") +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y") +
  theme_minimal() +
  labs(title = "Biomass Proportion by Taxa",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Biomass Proportion",
       fill = "Taxa") +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))


# Analysis ----------------------------------------------------------------

# Taxa Composition by Sample ----------------------------------------------


#Biomass
# Prepare data for plotting
plot_data <- zooscan_by_sample %>%
  mutate(Sample_ID_short = sample_id,  # Assuming sample_id is the short version
         biomass_prop = dryweight_C_mg_sum_taxa / dryweight_C_mg_sum_sample) %>%
  select(Sample_ID_short, object_annotation_category, biomass_prop)

# Create the stacked bar plot
ggplot(plot_data, aes(x = Sample_ID_short, 
                      y = biomass_prop, 
                      fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample ID", 
       y = "Proportion of Biomass", 
       fill = "Taxa Category",
       title = "Proportion of Biomass by Sample and Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Abundance
# Custom palette of 10 soft, high-contrast colors
custom_colors <- c(
  "#66C2A5", # Soft green
  "#FC8D62", # Soft orange
  "#8DA0CB", # Soft blue
  "#E78AC3", # Soft pink
  "#A6D854", # Soft lime green
  "#FFD92F", # Soft yellow
  "#E5C494", # Soft beige
  "#B3B3B3", # Soft grey
  "#A6CEE3", # Soft light blue
  "#1F78B4"  # Soft dark blue
)

# Compute relative abundances and filter to top 10 categories
relative_abundances <- zooscan_processed %>%
  group_by(sample_id, size_fraction, object_annotation_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(sample_id) %>%
  mutate(total = sum(count),
         relative_abundance = count / total) %>%
  ungroup() %>%
  mutate(Sample_ID_short = sample_id)  # Assuming sample_id is the short version

# Identify top 10 most abundant categories
top_categories <- relative_abundances %>%
  group_by(object_annotation_category) %>%
  summarise(total_abundance = sum(relative_abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(object_annotation_category)

# Filter to top 10 categories
relative_abundances_top10 <- relative_abundances %>%
  filter(object_annotation_category %in% top_categories)

# Create the stacked bar plot
ggplot(relative_abundances_top10, aes(x = Sample_ID_short, 
                                      y = relative_abundance, 
                                      fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) + # Apply custom color palette
  labs(
    title = "Relative Abundance by Sample and Top 10 Taxa Categories",
    x = "Sample ID",
    y = "Relative Abundance",
    fill = "Taxa Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Rotate x-axis labels
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 10), 
    title = element_text(size = 12), 
    legend.position = "bottom" # Position legend at the bottom
  )

# -------------------------------------------------------------------------



# Compute Shannon Diversity Index by site
# Count the number of individuals for each taxa within each site
count_data <- zooscan_processed %>%
  group_by(sample_id,object_annotation_category,size_fraction) %>%
  summarise(count = n())

# Compute Shannon Diversity Index by site
shannon_index <- count_data %>%
  group_by(sample_id,size_fraction) %>%
  summarise(shannon_index = diversity(count, index = "shannon"))

shannon_index_plot=shannon_index %>%
  #Make size class a factor
  mutate(size_fraction=factor(size_fraction, levels = c("0.2-1", "1-2", ">5"))) %>%
  left_join(metadata, by=c("sample_id","size_fraction")) %>% 
  distinct(.) 


#Shannon Index Plot
my_palette=custom_pallete()
shannon_index_plot %>%
  ggplot(aes(x = PC1, y = shannon_index, color = as.factor(size_fraction))) +
  geom_point(size=8, aes(shape=cycle, fill=cycle))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  scale_x_continuous(breaks = seq(-6, 7, by = 2))+
  scale_fill_manual(values = my_palette) +
  facet_wrap(~size_fraction, nrow=3)+
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="COI") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 4, label.y = 3.5)+
  theme_classic()





# Scale Comparison: Zooscan Biomass vs. Volume Filtered -------------------
#Plot grouped bar plot
labels_for_map=zooscan_by_sample %>% 
  ungroup()%>%
  select(sample_id,PC1) %>%
  unique(.) %>%
  arrange((PC1))


library(ggplot2)
library(dplyr)
library(ggpubr)

# Filter data and remove rows with NA or infinite values
cleaned_data <- zooscan_by_sample %>%
  filter(size_fraction != ">2") %>%
  filter(object_annotation_category == "Calanoida") %>%
  filter(is.finite(sample_tot_vol) & is.finite(dryweight_C_mg_sum_sample * acq_sub_part))

# Check if there are enough finite observations for correlation calculation
if (nrow(cleaned_data) < 3) {
  stop("Not enough finite observations for correlation calculation.")
}

# Plotting
cleaned_data %>%
  ggplot(aes(x = sample_tot_vol, y = log10(dryweight_C_mg_sum_sample * acq_sub_part), fill = as.factor(size_fraction), color = as.factor(size_fraction), shape = cycle)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
  scale_fill_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_color_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  stat_cor(method = "spearman", label.x = max(cleaned_data$sample_tot_vol) * 0.9, 
           label.y = max(log10(cleaned_data$dryweight_C_mg_sum_sample * cleaned_data$acq_sub_part)) * 0.9) +
  labs(x = "Volume Filtered (m^3)", y = "log10(Carbon Biomass) [g/m^3]", title = "Zooscan Sample Carbon Biomass vs. Volume Filtered", shape = "Cycle", color = "Size Fraction") + 
  theme_minimal() +
  facet_wrap(~size_fraction, nrow=3) +
  guides(size = FALSE, fill = FALSE) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



zooscan_by_sample %>% 
  filter(object_annotation_category=="Calanoida") %>% 
  ggplot(aes(x=sample_tot_vol, y=dryweight_C_mg_sum_taxa*acq_sub_part, color=size_fraction))+
  facet_wrap(~size_fraction )+
  stat_cor(method = "pearson", label.x = 300, label.y = 0.5)+
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +  # Add linear regression line  
  geom_point()+
  theme_minimal()








# Other TAXA --------------------------------------------------------------


# Euphausiids -------------------------------------------------------------
#Euphausiid DATA FRAME FOR PLOTTING
zoop_euphausiid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Euphausiacea")

#euphausiid dataframe for relative abundances
zoop_euphausiid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Euphausiacea")

#Save Biomass and relative abundances
write.csv(zoop_euphausiid_by_sample,here("PCR_bias_correction/data/Zooscan/zoop_euphausiid_by_sample_biomass.csv"))
write.csv(zoop_euphausiid_by_sample_relative_abundance,here("PCR_bias_correction/data/Zooscan/zoop_euphausiid_by_sample_relative_abundance.csv"))




# Oithonidae --------------------------------------------------------------

#Euphausiid DATA FRAME FOR PLOTTING
zoop_euphausiid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Euphausiacea")

#euphausiid dataframe for relative abundances
zoop_euphausiid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Euphausiacea")

#Save Biomass and relative abundances
write.csv(zoop_euphausiid_by_sample,here("PCR_bias_correction/data/Zooscan/zoop_euphausiid_by_sample_biomass.csv"))
write.csv(zoop_euphausiid_by_sample_relative_abundance,here("PCR_bias_correction/data/Zooscan/zoop_euphausiid_by_sample_relative_abundance.csv"))





####PLOTTING

#Color pallete
custom_palette <-  c("#FF6F61", "#FFA07A", "#7FB3D5", "#77DD77", "#B19CD9")  # Add more colors if needed


#==Proportions==#
relative_abundances_map %>%
  filter(object_annotation_category=="Calanoida") %>%
  filter(size_fraction!=">2") %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = cycle)) +
  geom_bar(stat = "identity") +
  labs(title = "Calanoid Copepod Zooscan Relative Abundances",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Relative Abundance"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  ylim(0, 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_props_zooscan
calanoid_props_zooscan





## ==== biomass === #
labels_for_map=biomass_map %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

zoop_calanoid_by_sample %>%
# biomass_map %>%
#   filter(object_annotation_category=="Calanoida") %>%
  # filter(Sample_ID_short %in% unique(dryweights$Sample_ID_short))
  filter(size_fraction!=">2") %>%
  ggplot(aes(x = as.factor(PC1), y = dryweight_C_mg_m2_taxa, fill = cycle)) +
  geom_bar(stat = "identity", width=0.8) +
  labs(title = "Calanoid Copepod Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = expression("Biomass (mg C " ~ m^-2 * ")"),
       fill = "Cycle") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm",">2 mm"))), scales = "free_y") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))+
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->calanoid_biomass_zooscan
calanoid_biomass_zooscan



#Continuous PC1
biomass_map %>%
  filter(object_annotation_category=="Calanoida") %>%
  ggplot(aes(x = PC1, y = dryweight_C, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack", width=0.2) +
  labs(title = "Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Biomass (ug C)",
       fill = "Taxa Group") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c(">2 mm","0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)









#With continuous PC1
###Zooscan relative abundances vs. pc1
##
labels_for_map=relative_abundances_map %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange(desc(PC1))

relative_abundances_map %>%
  filter(relative_abundance > 0.05, !(object_annotation_category %in% c("multiple organisms", "darksphere"))) %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = object_annotation_category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Zooscan Biomass",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Relative Abundance",
       fill = "Taxa Group") +
  facet_wrap(~size_fraction, nrow = 4, labeller = label_bquote(rows = .(c(">2 mm","0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)
