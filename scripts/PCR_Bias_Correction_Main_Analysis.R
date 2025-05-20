#Analysis Script for PCR Bias Correction Paper



# Read in the Data --------------------------------------------------------

# Packages and Functions --------------------------------------------------
librarian::shelf(tidyverse, googledrive, stringr,here,phyloseq,
                 extrafont, RColorBrewer, fido, compositions, ggpubr, patchwork, colorspace,
                 purrr, broom)

#Add functions for myself
source(("PCR_bias_correction/scripts/helpful_functions/treemap_funs_Capone.R"))
source("PCR_bias_correction/scripts/helpful_functions/phyloseq_mapping_funs.R")
source("PCR_bias_correction/scripts/helpful_functions/general_helper_functions.R")

saving=0


# Load in the data -----------------------------------------------------

#Metadata
metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

#Add depths df
depths=read.csv(here("PCR_bias_correction/data/physical_environmental_data/sample_depths.csv")) %>%
  select(-X) %>%
  mutate(Sample_ID_short=Sample_ID)


#Volume filtered (add to metadata)
volume_filtered=read.csv(here("PCR_bias_correction/data/raw_data/biomass/p2107_bt_volume_filtered.csv"))

#Dryweights
dryweights=read.csv("PCR_bias_correction/data/raw_data/biomass/dryweights_forzoopmetab.csv") %>%
  mutate(biomass_dry=8/3*biomass_dry) %>% 
  left_join(.,volume_filtered, by = c("Sample_ID_short"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>% 
  #Replace missing tow depths with 210
  mutate(depth = replace(depth, which(is.na(depth)), 210)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID) %>% 
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = size_fraction)  # Rename the existing column


env_metadata=metadata %>% 
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x)



# 18s ---------------------------------------------------------------------
#Using family

#Taxa file from pre-processed fido families for 18S
zhan_taxa=read.csv(here("PCR_bias_correction/data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))  %>% 
  select(-X) %>% 
  distinct() %>% 
  column_to_rownames("Family")



#Proportions

# Function to dynamically load the most recent file based on the naming pattern
load_most_recent_file <- function(suffix,primer) {
  # Directory path
  target_dir <- here("PCR_bias_correction/data/predicted_og")
  
  # List all files matching the pattern for the given suffix
  files <- list.files(target_dir, 
                      pattern = paste0("predicted_og_", primer, "_\\d{2}_\\d{2}_\\d{4}_", suffix, "_phy_all_and_subpools\\_nocollodaria.csv"), 
                      full.names = TRUE)
  # List all files matching the pattern for the given suffix
  # files <- list.files(target_dir, 
  #                     pattern = paste0("predicted_og_", primer, "_\\d{2}_\\d{2}_\\d{4}_", suffix, "_phy_all_and_subpools\\_nofilt.csv"), 
  #                     full.names = TRUE)
  # List all files matching the pattern for the given suffix
  # files <- list.files(target_dir, 
  #                     pattern = paste0("predicted_og_", primer, "_\\d{2}_\\d{2}_\\d{4}_", suffix, "_phy_all_and_subpools\\.csv"), 
  #                     full.names = TRUE)

  
  if (length(files) == 0) {
    stop(paste("No files found for suffix:", suffix))
  }
  
  # Extract dates from filenames
  extract_date <- function(file) {
    match <- regmatches(file, regexpr("\\d{2}_\\d{2}_\\d{4}", file))
    as.Date(match, format = "%m_%d_%Y")
  }
  
  # Find the most recent file
  files_with_dates <- data.frame(
    file = files,
    date = sapply(files, extract_date)
  )
  
  latest_file <- files_with_dates %>%
    arrange(desc(date)) %>%
    slice(1) %>%
    pull(file)
  
  # Read the file and apply transformations
  read.csv(latest_file) %>%
    select(-X) %>%
    mutate(coord = str_remove(coord, "^clr_"))
}

# Load files dynamically for s1, s2, and s3
fido_s1_18s <- load_most_recent_file("s1","18s")
fido_s2_18s <- load_most_recent_file("s2","18s")
fido_s3_18s <- load_most_recent_file("s3","18s")

#Merge
final_data_all_sizes_18s=rbind(fido_s1_18s,fido_s2_18s,fido_s3_18s) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#Make final dataframe
phy_taxa_pcr_18s= final_data_all_sizes_18s %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr_18s=phy_taxa_pcr_18s %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
  filter(Order=="Calanoida") 

#All taxa
pcr_join_prop_18s=phy_taxa_pcr_18s %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction_numeric,PC1,cycle,taxa) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))


# #Plotting ---------------------------------------------------------------

#Colors
taxa_colors_18s <- c(
  "other" = "grey",
  "Calanidae" = "#77DD77",       # Pastel green
  "Clausocalanidae" = "#72872d", # Bright pastel green
  "Eucalanidae" = "#89CFF0",     # Baby blue
  "Metridinidae" = "#9370DB",    # Pastel purple
  "Rhincalanidae" = "#4682B4",   # Steel blue
  "Paracalanidae" = "#eb6098",   # Lime pastel green
  "Oithonidae" = "#f0ca62",      # Light blue
  "Euphausiidae" = "#FF6961",    # Pastel red
  "Salpidae" = "#FFB6C1",        # Pastel pink
  "unidentified Calanoida" = "#2d8087",   
  "unidentified Collodaria" = "#912330",    # Pastel orange
  "unidentified Siphonophorae" = "#CDA4DE",  # Pastel lavender
  "Doliolidae" = "#27416b",
  "Sphaerozoidae"= "#7d0180"
)


# 18s ---------------------------------------------------------------------


# =========== Raw Reads Relative Abundance Data using taxa that went into fido model


# 18S ---------------------------------------------------------------------
fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools_nocollodaria.csv")) %>%
# fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools_nofilt.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)

fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools_nocollodaria.csv")) %>%
# fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools_nofilt.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) 

fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools_nocollodaria.csv")) %>%
# fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools_nofilt.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Family", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Family", all = TRUE)%>%
  column_to_rownames("Family") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw



#Metadata
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")

#Make phyloseq objects
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)

#Raw in counts
phy_18s_raw_counts=phyloseq(OTU, TAX, meta) %>% 
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>% 
  select(-asv_code)


#Raw in Proportions
phy_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  select(-asv_code)



#Filter to calanoid copepods
taxa_sel="Calanoida"


#2. Join with PCR Proportions
#All taxa
pcr_and_raw_18s_all=
  phy_18s %>% 
  mutate(taxa=Family) %>% 
  group_by(Sample_ID,taxa) %>%
  summarise(n_reads_raw=sum(n_reads)) %>% 
  left_join(pcr_join_prop_18s, by=c("Sample_ID","taxa"))

#3. Join in counts
phy_18s_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join_prop_18s, by="Sample_ID") %>%
  select(-size_fraction_numeric.x) %>%
  rename(size_fraction_numeric=size_fraction_numeric.y)->pcr_and_raw_18s_counts

#Add PCR-bias corrected counts
counts_raw_all_18s=phy_18s_raw_counts %>% 
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw_sum=sum(n_reads)) %>% 
  select(Sample_ID,size_fraction_numeric,n_reads_raw_sum)


pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
  left_join(.,counts_raw_all_18s,by=c("Sample_ID","size_fraction_numeric")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)




#Make phyloseq objects


# 18s ---------------------------------------------------------------------
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)

#Raw in counts
phy_18s_raw_counts=phyloseq(OTU, TAX, meta) %>% 
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>% 
  select(-asv_code)


#Raw in Proportions
phy_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  select(-asv_code)



#Filter to calanoid copepods
taxa_sel="Calanoida"

#2. Join with PCR Proportions
#All taxa
pcr_and_raw_18s_all=
  phy_18s %>% 
  mutate(taxa=Family) %>% 
  group_by(Sample_ID,taxa) %>%
  summarise(n_reads_raw=sum(n_reads)) %>% 
  left_join(pcr_join_prop_18s, by=c("Sample_ID","taxa"))

#3. Join in counts
phy_18s_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join_prop_18s, by="Sample_ID") %>%
  select(-size_fraction_numeric.x) %>%
  rename(size_fraction_numeric=size_fraction_numeric.y)->pcr_and_raw_18s_counts

#Add PCR-bias corrected counts
counts_raw_all=phy_18s_raw_counts %>% 
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw_sum=sum(n_reads)) %>% 
  select(Sample_ID,size_fraction_numeric,n_reads_raw_sum)


pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
  left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction_numeric")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)







# Q1. Raw Reads vs. PCR Bias Corrected ------------------------------------

#18s
plot_data_18s=pcr_and_raw_18s_all %>%
  # Create a new variable Sample_ID_short (everything before the first '_')
  mutate(Sample_ID_short = sub("_.*", "", Sample_ID)) %>%
  # Ensure Sample_ID is ordered by increasing PC1
  mutate(Sample_ID_short = reorder(Sample_ID_short, PC1)) %>%
  group_by(Sample_ID_short, taxa, size_fraction_numeric) %>%
  # Summarize raw and PCR reads
  summarize(
    n_reads_raw = sum(n_reads_raw, na.rm = TRUE),
    n_reads_pcr = sum(n_reads_pcr, na.rm = TRUE),
    PC1=mean(PC1),
    .groups = "drop"
  ) %>%
  # Pivot longer for faceting
  pivot_longer(
    cols = c(n_reads_raw, n_reads_pcr),
    names_to = "Metric",  # Create a column to identify raw and PCR
    values_to = "Reads"   # Move the values into a single column
  ) %>%
  group_by(Sample_ID_short, Metric, size_fraction_numeric) %>%
  # Normalize reads within each Sample_ID for both metrics
  mutate(Proportion = Reads / sum(Reads, na.rm = TRUE))

# Create custom labels for size_fraction_numeric and Metric
facet_labels <- list(
  size_fraction_numeric = c(
    "0.2-0.5mm" = "0.5",
    "0.5-1mm" = "1",
    "1-2mm" = "2"
  ),
  Metric = c(
    "n_reads_pcr" = "PCR corrected",
    "n_reads_raw" = "Raw"
  )
)

# Define labeller function
custom_labeller <- labeller(
  size_fraction_numeric = c(
    "0.2" = "0.2-0.5mm",
    "0.5" = "0.5-1mm",
    "1" = "1-2mm"
  ),
  Metric = c(
    "n_reads_pcr" = "PCR corrected",
    "n_reads_raw" = "Raw"
  )
)


# PC1 Labels
labels_for_PC1 <- metadata %>% 
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange(PC1)


fig1 <- ggplot(plot_data_18s, aes(x = as.factor(PC1), y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) +
  scale_fill_manual(values = taxa_colors_18s) +
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
  labs(
    title = "",
    x = expression("Offshore  \u2190   PC1   \u2192  Onshore"),
    y = "Proportion Reads",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),         # ~6–8 pt
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8, face = "bold"),                # ~8–9 pt
    axis.title.y = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.8, "lines")
  )  # Reduce spacing for compact fit)  # Show legend only on first plot
  
fig1

# Define output directory
output_dir <- "PCR_bias_correction/figures/v0"

# Save the plot as PNG
ggsave(
  filename = file.path(output_dir, "figure_2_v0.png"),
  plot = fig1,
  dpi = 600,
  width = 7,
  height = 4.5,  # adjust if needed
  units = "in"
)

# ---- Save as PDF (same size in mm) ----
ggsave(
  filename = file.path(output_dir, "figure_2_v0.pdf"),
  plot = fig1,
  dpi = 600,
  width = 7,
  height = 4.5,
  units = "mm",
  device = cairo_pdf  # ensures font embedding
)

# 
# ## By Order to better visualize alognside zooscan
# # Define a unique color palette for each Order
# order_colors <- c(
#   "other" = "grey",
#   "Calanoida" = "#8DD3C7",       # Pastel green
#   "Cyclopoida" = "#89CFF0",      # Baby blue
#   "Euphausiacea" = "#BC80BD",    # Pastel red
#   "Collodaria" = "#2d8087"      # Deep blue-green
# )
# 
# # Modify plot_data to mutate `taxa` into `Family`, then join with `zhan_taxa` by Family
# plot_data_18s <- plot_data_18s %>%
#   rename(Family = taxa) %>%  # Rename taxa to Family
#   left_join(zhan_taxa %>% rownames_to_column("Family"), by = "Family") %>%  # Join with zhan_taxa
#   mutate(Order = ifelse(is.na(Order), "other", Order))  # Assign "other" for missing values
# 
# # Create the stacked bar plot by Order
# fig1_order <- ggplot(plot_data, aes(x = Sample_ID_short, y = Proportion, fill = Order)) +
#   geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
#   facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) +  # Apply custom labels
#   scale_fill_manual(values = order_colors) +  # Apply custom Order colors
#   scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
#   labs(
#     title = "",
#     x = "Sample ID",
#     y = "Proportion Reads",
#     fill = "Order"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
#     axis.title.x = element_text(size = 16),  # Larger axis title (O)
#     axis.title.y = element_text(size = 16),  # Larger y-axis title (O)
#     plot.title = element_text(size = 18, face = "bold"),  # Increased title size (O)
#     legend.text = element_text(size = 16),  # Larger legend text (L)
#     strip.text = element_text(face = "bold", size = 16), # Facet text size
#     legend.title = element_text(size = 16, face = "bold"),  # Larger, bold legend title (L)
#     legend.position ="bottom")  # Show legend only on first plot
# 
# # Print the plot
# fig1_order
# 
# # Define output directory
# output_dir <- "PCR_bias_correction/figures/v0"
# 
# # Save the plot as PNG
# ggsave(
#   filename = file.path(output_dir, "figure_2_v0_order.png"),
#   plot = fig1_order,
#   dpi = 600,
#   width = 18,
#   height = 7,
#   units = "in",
# )
# 
# # Save the plot as PDF
# ggsave(
#   filename = file.path(output_dir, "figure_2_v0_order.pdf"),
#   plot = fig1_order,
#   dpi = 600,
#   width = 18,
#   height = 7,
#   units = "in",
# )



# Q2. Taxa Effects of PCR Bias --------------------------------------------
#Use all taxa
# Load data
all_amp_effs_18s <- read.csv(here("PCR_bias_correction/data/amp_effs/all_amp_effs_18s_all_sub_nofilt.csv")) %>%
  mutate(Family = str_extract(Lambda.coord, "[^_]+$"))

# Calculate mean Lambda per Family and order them ascending
ordered_taxa <- all_amp_effs_18s %>%
  group_by(Family) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_lambda) %>%
  pull(Family)

# Create x-axis group (Family + size_fraction) and apply order
all_amp_effs_18s <- all_amp_effs_18s %>%
  mutate(
    Family = factor(Family, levels = ordered_taxa),
    size_fraction = factor(size_fraction, levels = c(0.2, 0.5, 1)),
    group = interaction(Family, size_fraction, lex.order = TRUE),
    x_label = ifelse(size_fraction == 0.5, as.character(Family), " ")
  )

# Create named vector of labels
x_labels_named <- all_amp_effs_18s %>%
  distinct(group, x_label) %>%
  deframe()
x_labels_named

# Order 'group' by Family order and within-Family size fraction
ordered_groups <- all_amp_effs_18s %>%
  arrange(Family, size_fraction) %>%
  pull(group) %>%
  unique()

# Apply ordered factor
all_amp_effs_18s <- all_amp_effs_18s %>%
  mutate(group = factor(group, levels = ordered_groups))

# Compute dynamic limits and breaks
lambda_range <- range(all_amp_effs_18s$Lambda.mean, na.rm = TRUE)
max_abs <- max(abs(lambda_range)) * 1.5
y_limits <- c(-max_abs, max_abs)
y_breaks <- seq(from = floor(y_limits[1] * 10) / 10,
                to   = ceiling(y_limits[2] * 10) / 10,
                by = 0.1)

# Final plot
amp_effs_all_and_subpools_by_taxa_18s_no_filt <- ggplot(all_amp_effs_18s, aes(x = group, y = Lambda.mean,
                                                                              color = Family,
                                                                              shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5), width = 0.2, size = 1.1) +
  geom_point(size = 4, stroke = 1.2) +
  scale_color_manual(values = taxa_colors_18s) +
  scale_shape_manual(name = "Size Fraction",
                     values = c("0.2" = 16, "0.5" = 17, "1" = 18),
                     labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_x_discrete(labels = x_labels_named)+
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 0.1, trim = TRUE)
  ) +
  labs(
    x = "",
    y = "Amplification Efficiency (CLR)",
    title = "",
    color = "Family"
  ) +
  theme_minimal() +
  guides(color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  )



# View the plot
amp_effs_all_and_subpools_by_taxa_18s_no_filt



# Save the plot as PNG
output_dir <- "PCR_bias_correction/figures/supporting_info/"
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_nofilt.png"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_no_filt,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_nofilt.pdf"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_no_filt,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)



# Figure 4b. Amplification efficiencies vs. Proportion reads -----------------
detach("package:fido", unload = TRUE)

# #NO filter --------------------------------------------------------------
# Function to sum columns in fido datasets containing "S1", "S2", or "S3" in their names
mean_s_columns_18s_family <- function(df, family, suffix) {
  matching_columns <- names(df)[str_detect(names(df), suffix)]
  if (length(matching_columns) > 0) {
    normalized_values <- sapply(matching_columns, function(col) {
      column_sum_hash <- sum(df[df$Family == family, col], na.rm = TRUE)
      column_sum <- sum(df[[col]], na.rm = TRUE)
      if (column_sum_hash == 0 || column_sum == 0) {
        return(0)
      }
      column_sum_hash / column_sum
    })
    return(mean(normalized_values, na.rm = TRUE))
  } else {
    return(NA)
  }
}

all_amp_effs_18s_nofilt <- read.csv(here("PCR_bias_correction/data/amp_effs/all_amp_effs_18s_all_sub_nofilt.csv")) %>%
  mutate(Family = str_extract(Lambda.coord, "[^_]+$"))



fido_18s_s1_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools_nofilt.csv"), header = TRUE, check.names = FALSE)
fido_18s_s2_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools_nofilt.csv"), header = TRUE, check.names = FALSE)
fido_18s_s3_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools_nofilt.csv"), header = TRUE, check.names = FALSE)


# 18s ---------------------------------------------------------------------


# Create the final dataframe for S1, S2, and S3
result_s1_nofilt <- all_amp_effs_18s_nofilt %>%
  filter(str_detect(pool, "S1")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s1_nofilt, Family, "S1"),
         SizeFraction = "S1") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s2_nofilt <- all_amp_effs_18s_nofilt %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s2_nofilt, Family, "S2"),
         SizeFraction = "S2") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s3_nofilt <- all_amp_effs_18s_nofilt %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s3_nofilt, Family, "S3"),
         SizeFraction = "S3") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

# Combine S1, S2, S3 with updated mean calc and SizeFraction tag
result_combined_nofilt <- bind_rows(
  all_amp_effs_18s_nofilt %>% filter(str_detect(pool, "S1")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_18s_family(fido_18s_s1_nofilt, Family, "S1"),
           SizeFraction = "S1") %>%
    select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup(),
  all_amp_effs_18s_nofilt %>% filter(str_detect(pool, "S2")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_18s_family(fido_18s_s2_nofilt, Family, "S2"),
           SizeFraction = "S2") %>%
    select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup(),
  all_amp_effs_18s_nofilt %>% filter(str_detect(pool, "S3")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_18s_family(fido_18s_s3_nofilt, Family, "S3"),
           SizeFraction = "S3") %>%
    select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup()
)

# Reorder Family by increasing Lambda.mean
family_order <- result_combined_nofilt %>%
  group_by(Family) %>%
  summarise(mean_Lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_Lambda) %>%
  pull(Family)


# Compute correlation stats per size fraction
lm_stats_by_size <- summarized_result %>%
  filter(Mean > 0) %>%
  group_by(SizeFraction) %>%
  summarise(
    fit = list(lm(log2(Mean) ~ Lambda.mean)),
    .groups = "drop"
  ) %>%
  mutate(
    r2 = map_dbl(fit, ~ summary(.x)$r.squared),
    p = map_dbl(fit, ~ summary(.x)$coefficients[2, 4]),
    annotation = paste0(
      "R² = ", formatC(r2, digits = 2, format = "f"), "\n",
      "P = ", ifelse(p < 0.001, "<0.001", formatC(p, digits = 2, format = "f"))
    )
  )

# Summarize and set factor order
summarized_result <- result_combined_nofilt %>%
  group_by(Family, SizeFraction) %>%
  summarise(
    Lambda.mean = mean(Lambda.mean, na.rm = TRUE),
    Mean = mean(Mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Family = factor(Family, levels = family_order))%>%
  mutate(Family = factor(Family, levels = names(taxa_colors_18s)))

# Correlation
pearson_18s <- cor.test(~ Lambda.mean + Mean, data = summarized_result)
cat("Pearson's R:", pearson_18s$estimate, "\n")
cat("p-value:", pearson_18s$p.value, "\n\n")


size_fraction_labels <- c(
  "S1" = "0.2–0.5 mm",
  "S2" = "0.5–1 mm",
  "S3" = "1–2 mm"
)

# Plot
amp_effs_vs_rra_18s_nofilt <- summarized_result %>%
  filter(Mean > 0) %>%
  ggplot(aes(x = Lambda.mean, y = log2(Mean))) +
  geom_point(
    aes(fill = Family, shape = SizeFraction),
    size = 4,
    stroke = 0.8,
    color = "black",
    alpha = 0.85
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8)+
  geom_text(
    data = lm_stats_by_size,
    aes(x = -Inf, y = Inf, label = annotation),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 4, fontface = "italic"
  ) +
  scale_fill_manual(values = taxa_colors_18s) +
  scale_shape_manual(
    values = c(21, 22, 23),
    labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"),
    name = "Size Fraction"
  ) +
  labs(
    x = "Amplification Efficiency",
    y = expression(log[2]*"(mean relative read abundance)"),
    fill = "Order",  # <- update to 'fill' instead of 'color' since you're using fill
    shape = "Size Fraction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face="bold"),
    legend.key.size = unit(0.4, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1, "lines"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.spacing.x = unit(0.1, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")
  ) +
  guides(
    shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5),
    fill = guide_legend(
      override.aes = list(
        fill = unname(taxa_colors_18s),  # strip names from vector
        color = "black",
        shape = 21,
        stroke = 0.8,
        size = 4
      ),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 4
    )
  ) +
  facet_wrap(~SizeFraction, nrow = 1,
             labeller = labeller(SizeFraction = size_fraction_labels))

amp_effs_vs_rra_18s_nofilt

#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/v0/fig3b_amp_effs_vs_rra_18s_nofilt.png"),
  plot = amp_effs_vs_rra_18s_nofilt,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/v0/amp_effs_vs_rra_18s_nofilt.pdf"),
  plot = amp_effs_vs_rra_18s_nofilt,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)


#COmbine for figure 4
library(patchwork)

# Remove x-axis labels and ticks from the top plot
amp_effs_all_and_subpools_by_taxa_18s_no_filt <- amp_effs_all_and_subpools_by_taxa_18s_no_filt +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"  # this removes the legend entirely from the top plot
  )

amp_effs_vs_rra_18s_nofilt <- amp_effs_vs_rra_18s_nofilt +
  guides(
    shape = guide_legend(
      ncol = 1,  # this forces one column
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    ),
    fill = guide_legend(
      override.aes = list(
        fill = unname(taxa_colors_18s),
        color = "black",
        shape = 21,
        stroke = 0.8,
        size = 4
      ),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 4
    )
  )
  

# Combine the plots
figure_s4_combined <-  wrap_elements(amp_effs_all_and_subpools_by_taxa_18s_no_filt) /
  amp_effs_vs_rra_18s_nofilt +
  plot_layout(heights = c(1.1, 1), guides = "keep") +
  plot_annotation(tag_levels = 'A')  # This will label them "A", "B"

figure_s4_combined

# Save the combined figure
ggsave(
  filename = file.path(output_dir, "figure_4_combined.png"),
  plot = figure_s4_combined,
  dpi = 600,
  width = 7,
  height = 8,  # Adjust height to accommodate both plots vertically
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "figure_4_combined.pdf"),
  plot = figure_s4_combined,
  dpi = 600,
  width = 7,
  height = 8,
  units = "in"
)


# Q3. Comparison with Zooscan for calanoids-------------------------------------------------
# Zooscan -----------------------------------------------------------------


#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read in processed Zooscan data and look at biomass proortion
zoo_metric="biomass_prop"



#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Euphausiacea"


#Zooscan biomass proportion
# All taxa
zooscan_all=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass_esd.csv"), row.names = 1)%>%
  # select(-X, acq_min_mesh) %>% 
  # mutate(Sample_ID=sample_id) %>%  
  distinct(.)


#Calanoida
zooscan_taxa=zooscan_all%>%
  filter(object_annotation_category==taxa_sel)  



#18S
pcr_and_raw_18s_calanoids=pcr_and_raw_18s_all %>% 
  mutate(Family=taxa) %>% 
  left_join(.,zhan_taxa %>% rownames_to_column("Family") %>% 
              select(Family, Order), by="Family") %>% 
  filter(Order==taxa_sel) %>% 
  select(-Family,-taxa)%>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric,cycle) %>% 
  summarise(n_reads_raw=sum(n_reads_raw), 
            n_reads_pcr=sum(n_reads_pcr),
            PC1=mean(PC1)) 
  
#Propotions
pcr_raw_zoo_18s=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(pcr_and_raw_18s_calanoids, by=c("PC1","size_fraction_numeric")) %>% 
  unique(.)





# Correlations. Make a grid plot of Y axis is:  ----------------------------
# Raw reads or PCR corrected
# X is: Zooscan biomass, abundance or relative of each

facet_x_labels <- c(
  "biomass_prop_taxa" = "Log10(Biomass Proportion)",
  "dryweight_C_mg_m2_taxa_mean" = "Log10 (Dry Weight (mg Cm²))",
  "relative_abundance" = "Log10(Relative Abundance)",
  "abundance_m2_mean" = "Log10(Abundance (m²))"
)

facet_y_labels <- c(
  "n_reads_raw" = "Raw Reads",
  "n_reads_pcr" = "PCR Bias-Corrected Reads"
)

size_fraction_labels <- c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")

# Convert data to long format and filter for log transformations
pcr_raw_zoo_18s_long <- pcr_raw_zoo_18s %>%
  filter(!is.na(cycle)) %>%
  pivot_longer(cols = c(biomass_prop_taxa, dryweight_C_mg_m2_taxa_mean, relative_abundance, abundance_m2_mean), 
               names_to = "x_variable", values_to = "x_value") %>%
  pivot_longer(cols = c(n_reads_raw, n_reads_pcr), 
               names_to = "y_variable", values_to = "y_value") %>%
  filter(x_value > 0, y_value >= 0) %>%  # Ensure valid sqrt transform (y_value >= 0)
  mutate(
    x_variable = factor(x_variable, levels = names(facet_x_labels), labels = facet_x_labels),
    y_variable = factor(y_variable, levels = names(facet_y_labels), labels = facet_y_labels),
    size_fraction = factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2"))
  )


# Total number of unique comparisons
n_tests <- pcr_raw_zoo_18s_long %>%
  distinct(size_fraction, y_variable, x_variable) %>%
  nrow()

# Compute correlation per facet combination
cor_results_grid <- pcr_raw_zoo_18s_long %>%
  group_by(size_fraction, y_variable, x_variable) %>%
  summarise(
    cor_test = list(cor.test(log10(x_value), asin(sqrt(y_value)))),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(cor_test, ~ .x$estimate),
    p_uncorrected = map_dbl(cor_test, ~ .x$p.value),
    p_bonferroni = pmin(p_uncorrected * n_tests, 1),
    label = paste0(
      "r = ", formatC(r, digits = 2, format = "f"), "\n",
      "P = ", ifelse(p_bonferroni < 0.001, "<0.001", formatC(p_bonferroni, digits = 2, format = "f"))
    )
  )


# Create 6-row, 4-column facet grid plot
zoo_vs_pcr_grid <- ggplot(pcr_raw_zoo_18s_long, 
                          aes(x = ((x_value)), 
                              y = asin(sqrt(y_value)))) +  # Apply asin sqrt transformation
  geom_point(aes(shape = cycle, size = 3, color = as.factor(size_fraction), fill = as.factor(size_fraction))) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                    labels = size_fraction_labels) +
  scale_color_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                     labels = size_fraction_labels) +
  facet_grid(rows = vars(size_fraction, y_variable), cols = vars(x_variable), scales = "free") +  # 6x4 layout
  labs(
    x = "Log10 Transformed X-Axis Variables",
    y = "asin(sqrt(Read Counts))",
    shape = "Cycle",
    color = "Size Fraction"
  ) +
  geom_text(
    data = cor_results_grid,
    aes(label = label),
    x = 0.25, y = 1.2,
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "italic"
  )  +  # Pearson correlation
  guides(size = FALSE, fill = FALSE) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Major grid lines for readability
    panel.grid.minor = element_line(color = "gray90", size = 0.2),  # Minor grid lines for subtle separation
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around each facet
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Light gray background for facet labels
    strip.text = element_text(size = 16, face = "bold"),  # Larger and bold facet labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


# Print plot
zoo_vs_pcr_grid


#Just Absolute biomass
# Filter and reshape for plotting
pcr_vs_dryweight <- pcr_raw_zoo_18s %>%
  filter(!is.na(cycle)) %>%
  select(sample_id, cycle, size_fraction, dryweight_C_mg_m2_taxa_mean, n_reads_raw, n_reads_pcr) %>%
  pivot_longer(cols = c(n_reads_raw, n_reads_pcr),
               names_to = "y_variable", values_to = "y_value") %>%
  filter(dryweight_C_mg_m2_taxa_mean > 0, y_value >= 0) %>%
  mutate(
    y_variable = factor(y_variable, levels = c("n_reads_raw", "n_reads_pcr"),
                        labels = c("Raw Reads", "PCR-Corrected Reads")),
    size_fraction = factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2"))
    
  )

# Number of comparisons = number of facet combinations
n_tests <- pcr_vs_dryweight %>%
  distinct(y_variable, size_fraction) %>%
  nrow()

# Compute Pearson correlation per facet
cor_results <- pcr_vs_dryweight %>%
  group_by(y_variable, size_fraction) %>%
  summarise(
    cor_test = list(cor.test(
      ~ log10(dryweight_C_mg_m2_taxa_mean) + asin(sqrt(y_value))
    )),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(cor_test, ~ .x$estimate),
    p_uncorrected = map_dbl(cor_test, ~ .x$p.value),
    p_bonferroni = pmin(p_uncorrected * n_tests, 1),  # Bonferroni correction
    label = paste0("r = ", formatC(r, digits = 2, format = "f"), 
                   "\nP = ", ifelse(p_bonferroni < 0.001, "<0.001", formatC(p_bonferroni, digits = 2, format = "f")))
  )

#Define pastel color palette
pastel_colors <- c("0.2-0.5" = "#AEDFF7", "0.5-1" = "#B6E3B6", "1-2" = "#FFB3AB")

# Plot
ggplot(pcr_vs_dryweight,
       aes(x = log10(dryweight_C_mg_m2_taxa_mean),
           y = asin(sqrt(y_value)),
           fill = size_fraction)) +
  geom_point(aes(shape = cycle), size = 6, stroke = 1.2, color = "black") +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
  scale_fill_manual(values = pastel_colors, labels = size_fraction_labels, name = "Size Fraction") +
  coord_cartesian(xlim = c(0, 3)) +
  facet_grid(rows = vars(y_variable), cols = vars(size_fraction), scales = "free_y") +
  labs(
    x = "Log10 Dry Weight C (mg)",
    y = "arcsine-sqrt Relative Read Abundance",
    shape = "Cycle",
    fill = "Size Fraction"
  ) +
  geom_text(
    data = cor_results,
    aes(label = label),
    x = 0.25, y = 0.5,
    inherit.aes = FALSE,
    size = 4, fontface = "italic"
  )  +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6, stroke = 1.2)),
    shape = guide_legend(override.aes = list(fill = "gray90", color = "black")),
    color = "none"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(size = 16, face = "bold", family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold", family = "serif"),
    axis.text.y = element_text(size = 14, face = "bold", family = "serif"),
    axis.title = element_text(size = 16, face = "bold", family = "serif"),
    legend.title = element_text(size = 16, face = "bold", family = "serif"),
    legend.text = element_text(size = 14, family = "serif")
  )-> methods_correlation_plot

methods_correlation_plot
  #PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_18s_scatter.png"),
  plot = methods_correlation_plot,
  dpi = 600,
  width = 12, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_18s_scatter.pdf"),
  plot = methods_correlation_plot,
  dpi = 600,
  width = 12, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)



#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions

#Format Long and add difference metrics
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s, 
                                     cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance_metric") %>%
  # filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance_metric[Method == "n_reads_pcr"] - relative_abundance_metric),
         diff_raw = abs(relative_abundance_metric[Method == "n_reads_raw"] - relative_abundance_metric))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance_metric, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance_metric, "xx",NA))


# Plot with secondary y-axis
ggplot(pcr_raw_zoo_18s_long, aes(x = as.factor(PC1), fill = Method)) +
  geom_bar(data = filter(pcr_raw_zoo_18s_long, Method != "dryweight_C_mg_sum_taxa"),
           aes(y = relative_abundance_metric),
           stat = "identity",
           position = position_dodge(width = 0.9)) +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y") +
  scale_y_continuous(
    name = "Proportion Reads or Biomass"
  ) +
  scale_fill_manual(values = c("#70BF41", "#4F86F7", "#F78D4F", "grey50"),
                    labels = c("Zooscan Biomass", "PCR-corrected", "Raw Reads", "Dryweight C")) +
  theme_minimal() +
  labs(x = "",
       fill = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "right") +
  geom_text(aes(y = relative_abundance_metric, label = is_closer),
            position = position_dodge(width = 0.9), vjust = 0, size = 4) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)->grouped_bar_all_18s




grouped_bar_all_18s
#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/calanoid_methods_compare_18s.png"),
  plot = grouped_bar_all_18s,
  dpi = 600,
  width =6, # Adjust width (inches) based on journal requirements
  height = 4.5, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/calanoid_methods_compare_18s.pdf"),
  plot = grouped_bar_all_18s,
  dpi = 600,
  width =6, # Adjust width (inches) based on journal requirements
  height = 4.5, # Adjust height (inches) based on journal requirements
  units = "in"
)


#Scatter Plot
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle))%>%
  ggplot(.,aes(x=((log10(dryweight_C_mg_m2_taxa))), y=asin(sqrt(n_reads_pcr))))+
  geom_point(aes(shape=cycle, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3, scale="free_y") +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))->zoo_vs_pcr
zoo_vs_pcr



# Mean squared error ---------------------------------------------------------------------
pcr_raw_zoo_18s_mse=pcr_raw_zoo_18s%>% 
  # filter(!is.na(n_reads_raw)) %>% 
  group_by(size_fraction) %>% 
  mutate(se_zoo_pcr = ((biomass_prop_taxa) - (n_reads_pcr))^2,
         se_zoo_raw = ((biomass_prop_taxa) - (n_reads_raw))^2,
         average_se = (se_zoo_pcr + se_zoo_raw) / 2)%>%
  mutate(difference = se_zoo_raw - se_zoo_pcr,
         colorr = ifelse(difference >= 0, "#f5776e", "#8d9af2")) %>% 
  ungroup()

hist(pcr_raw_zoo_18s_mse$se_zoo_pcr)
hist(pcr_raw_zoo_18s_mse$se_zoo_raw)

# Calculate average MSE per size
sum_mse <- pcr_raw_zoo_18s_mse %>%
  group_by(size_fraction) %>%
  summarise(
    average_mse_by_size_pcr = median(se_zoo_pcr, na.rm = TRUE),
    average_mse_by_size_raw = median(se_zoo_raw, na.rm = TRUE)
  )

sum_mse

# Error barplot
ggplot(pcr_raw_zoo_18s_mse, aes(x = as.factor(PC1), y = difference, fill = colorr)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("#f5776e" = "#f5776e", "#8d9af2" = "#8d9af2")) +
  theme_minimal() +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y", 
             labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  labs(title = "",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Difference in Square-Error (RRA-PCR-RA)",
       fill="") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none")  +
  geom_text(data = sum_mse, aes(x = 0.5, y = -0.2, label = sprintf("MedSE PCR-RA = %.3f\nMedSE RRA = %.3f", average_mse_by_size_pcr, average_mse_by_size_raw)),
            hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE)->mse_plot

mse_plot
#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/mse_plot_methods_compare_18s.png"),
  plot = mse_plot,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 4.5, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/mse_plot_methods_compare_18s.pdf"),
  plot = mse_plot,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 4.5, # Adjust height (inches) based on journal requirements
  units = "in"
)


#Combine Figure S5
# Remove x-axis labels and ticks from the top plot
grouped_bar_all_18s=grouped_bar_all_18s +
  theme(
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.position = "right"  # this removes the legend entirely from the top plot
  )

# Combine the plots
figure_s5_combined <-  grouped_bar_all_18s +
  mse_plot +
  plot_layout(heights = c(1.1, 1), guides = "keep") +
  plot_annotation(tag_levels = 'A')  # This will label them "A", "B"

figure_s5_combined

# Save the combined figure
ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/figure_s5_combined.png"),
  plot = figure_s5_combined,
  dpi = 600,
  width = 7,
  height = 8,  # Adjust height to accommodate both plots vertically
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/figure_s5_combined.pdf"),
  plot = figure_s5_combined,
  dpi = 600,
  width = 7,
  height = 8,
  units = "in"
)


# Correlations ------------------------------------------------------------

pcr_raw_zoo_18s %>%
  filter(!is.na(cycle))%>%
  # filter(size_fraction_numeric !=0.2) %>% 
  ggplot(.,aes(x=log10((abundance_m2)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  # facet_wrap(~size_fraction, nrow=3, scales="free_x") +
  labs(x = "Zooscan Biomass", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))


# COI ---------------------------------------------------------------------




pcr_and_raw_coi_calanoids=pcr_and_raw_coi_all %>% 
  mutate(Genus=taxa) %>% 
  left_join(.,coi_taxa %>% rownames_to_column("Genus") %>% 
              select(Genus, Order), by="Genus") %>% 
  filter(Order=="Calanoida") %>% 
  select(-Genus,-taxa)%>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric,cycle) %>% 
  summarise(n_reads_raw=sum(n_reads_raw), 
            n_reads_pcr=sum(n_reads_pcr),
            PC1=mean(PC1)) 

#Propotions
pcr_raw_zoo_coi=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(pcr_and_raw_coi_calanoids, by=c("PC1","size_fraction_numeric")) %>% 
  unique(.)


#Format Long and add difference metrics
pcr_raw_zoo_coi_long <- pivot_longer(pcr_raw_zoo_coi, 
                                     cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance_metric") %>%
  # filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance_metric[Method == "n_reads_pcr"] - relative_abundance_metric),
         diff_raw = abs(relative_abundance_metric[Method == "n_reads_raw"] - relative_abundance_metric))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance_metric, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance_metric, "xx",NA))



# Correlations. Make a grid plot of Y axis is:  ----------------------------
# Raw reads or PCR corrected
# X is: Zooscan biomass, abundance or relative of each
# Define facet labels
facet_x_labels <- c(
  "biomass_prop_taxa" = "Log10(Biomass Proportion)",
  "dryweight_C_mg_m2_taxa" = "Log10 (Dry Weight (mg C m²))",
  "count" = "Log10(Relative Abundance)",
  "abundance_m2" = "Log10(Abundance (m²))"
)

facet_y_labels <- c(
  "n_reads_raw" = "Raw Reads",
  "n_reads_pcr" = "PCR Bias-Corrected Reads"
)

size_fraction_labels <- c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")

# Convert data to long format and filter for log transformations
pcr_raw_zoo_coi_long <- pcr_raw_zoo_coi %>%
  filter(!is.na(cycle)) %>%
  pivot_longer(cols = c(biomass_prop_taxa, dryweight_C_mg_m2_taxa, count, abundance_m2), 
               names_to = "x_variable", values_to = "x_value") %>%
  pivot_longer(cols = c(n_reads_raw, n_reads_pcr), 
               names_to = "y_variable", values_to = "y_value") %>%
  filter(x_value > 0, y_value >= 0) %>%  # Ensure valid sqrt transform (y_value >= 0)
  mutate(
    x_variable = factor(x_variable, levels = names(facet_x_labels), labels = facet_x_labels),
    y_variable = factor(y_variable, levels = names(facet_y_labels), labels = facet_y_labels),
    size_fraction = factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2"))
  )

# Create 6-row, 4-column facet grid plot
zoo_vs_pcr_grid <- ggplot(pcr_raw_zoo_coi_long, 
                          aes(x = log10(x_value), 
                              y = asin(sqrt(y_value)))) +  # Apply asin sqrt transformation
  geom_point(aes(shape = cycle, size = 3, color = as.factor(size_fraction), fill = as.factor(size_fraction))) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                    labels = size_fraction_labels) +
  scale_color_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                     labels = size_fraction_labels) +
  facet_grid(rows = vars(size_fraction, y_variable), cols = vars(x_variable), scales = "free") +  # 6x4 layout
  labs(
    x = "Log10 Transformed X-Axis Variables",
    y = "asin(sqrt(Read Counts))",
    shape = "Cycle",
    color = "Size Fraction"
  ) +
  stat_cor(method = "pearson", label.x = -1, label.y = 1) +  # Pearson correlation
  guides(size = FALSE, fill = FALSE) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Major grid lines for readability
    panel.grid.minor = element_line(color = "gray90", size = 0.2),  # Minor grid lines for subtle separation
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around each facet
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Light gray background for facet labels
    strip.text = element_text(size = 16, face = "bold"),  # Larger and bold facet labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


# Print plot
(zoo_vs_pcr_grid)



#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions


pcr_raw_zoo_coi_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance_metric, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y") +
  theme_minimal() +
  labs(title = "",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Biomass ",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "right") +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = 0, size = 9) +
  scale_fill_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"),
                    labels = c("Zooscan Biomass", "PCR-corrected", "Raw Reads")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  ylim(0, 1.1)+
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+# Define linetypes for each Method
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F")) -> grouped_bar_all_coi  # Define colors for each Method

grouped_bar_all_coi
#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_coi.png"),
  plot = grouped_bar_all_coi,
  dpi = 600,
  width = 10, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_coi.pdf"),
  plot = grouped_bar_all_coi,
  dpi = 600,
  width = 10, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)

#Scatter Plot
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle))%>%
  ggplot(.,aes(x=log10(((dryweight_C_mg_m2_taxa))), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  # facet_wrap(~size_fraction, nrow=3, scale="free_y") +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))->zoo_vs_pcr
zoo_vs_pcr







# Supporting Information --------------------------------------------------



#Figure S


# Load data
all_amp_effs_18s <- read.csv(here("PCR_bias_correction/data/amp_effs/all_amp_effs_18s_all_sub.csv")) %>%
  mutate(Family = str_extract(Lambda.coord, "[^_]+$"))

# Calculate mean Lambda per Family and order them ascending
ordered_taxa <- all_amp_effs_18s %>%
  group_by(Family) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_lambda) %>%
  pull(Family)

# Create x-axis group (Family + size_fraction) and apply order
all_amp_effs_18s <- all_amp_effs_18s %>%
  mutate(
    Family = factor(Family, levels = ordered_taxa),
    size_fraction = factor(size_fraction, levels = c(0.2, 0.5, 1)),
    group = interaction(Family, size_fraction, lex.order = TRUE),
    x_label = ifelse(size_fraction == 0.5, as.character(Family), " ")
  )

# Create named vector of labels
x_labels_named <- all_amp_effs_18s %>%
  distinct(group, x_label) %>%
  deframe()
x_labels_named

# Order 'group' by Family order and within-Family size fraction
ordered_groups <- all_amp_effs_18s %>%
  arrange(Family, size_fraction) %>%
  pull(group) %>%
  unique()

# Apply ordered factor
all_amp_effs_18s <- all_amp_effs_18s %>%
  mutate(group = factor(group, levels = ordered_groups))

# Compute dynamic limits and breaks
lambda_range <- range(all_amp_effs_18s$Lambda.mean, na.rm = TRUE)
max_abs <- max(abs(lambda_range)) * 3
y_limits <- c(-max_abs, max_abs)
y_breaks <- seq(from = floor(y_limits[1] * 10) / 10,
                to   = ceiling(y_limits[2] * 10) / 10,
                by = 0.05)

# Final plot
amp_effs_all_and_subpools_by_taxa_18s <- ggplot(all_amp_effs_18s, aes(x = group, y = Lambda.mean,
                                                                      color = Family,
                                                                      shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5), width = 0.2, size = 1.1) +
  geom_point(size = 4, stroke = 1.2) +
  scale_color_manual(values = taxa_colors_18s) +
  scale_shape_manual(name = "Size Fraction",
                     values = c("0.2" = 16, "0.5" = 17, "1" = 18),
                     labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_x_discrete(labels = x_labels_named)+
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 0.1, trim = TRUE)
  ) +
  labs(
    x = "",
    y = "Amplification Efficiency (CLR)",
    title = "",
    color = "Family"
  ) +
  theme_minimal() +
  guides(color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  )

# View the plot
amp_effs_all_and_subpools_by_taxa_18s



# Save the plot as PNG
output_dir <- "PCR_bias_correction/figures/v0"
ggsave(
  filename = file.path(output_dir, "amp_effs_18s.png"),
  plot = amp_effs_all_and_subpools_by_taxa_18s,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "amp_effs_18s.pdf"),
  plot = amp_effs_all_and_subpools_by_taxa_18s,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)


fido_18s_s1 = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)
fido_18s_s2 = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)
fido_18s_s3 = read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)




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
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s1, Family, "S1"),
         SizeFraction = "S1") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s2 <- all_amp_effs_18s %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s2, Family, "S2"),
         SizeFraction = "S2") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s3 <- all_amp_effs_18s %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s_family(fido_18s_s3, Family, "S3"),
         SizeFraction = "S3") %>%
  select(Family, pool, Lambda.mean, Mean, SizeFraction) %>%
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
  summarise(Lambda.mean = mean(Lambda.mean, na.rm = TRUE), Mean= mean(Mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Family = factor(Family, levels = family_order))

pearson_18s <- cor.test(~ Lambda.mean + Mean, data = summarized_result)

cat("18S Dataset:\n")
cat("Pearson's R:", pearson_18s$estimate, "\n")
cat("p-value:", pearson_18s$p.value, "\n\n")

summarized_result %>% 
  filter(Mean > 0) %>%
  ggplot(., aes(x = Lambda.mean, y = log2(Mean))) +
  geom_point(aes(color = Family, shape = SizeFraction), size = 10) +
  geom_smooth(method="lm")+
  theme_classic() +
  labs(title = "",
       x = "Amplification Efficiency",
       y = "log2(Mean Relative Read Abundance)",
       color = "Family",
       shape = "Size Fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "right") +
  scale_color_manual(values = taxa_colors_18s) +
  scale_shape_manual(values = c(16, 17, 18), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")) -> amp_effs_vs_rra_18s
amp_effs_vs_rra_18s

#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/amp_effs_vs_rra_18s_nocollodaria.png"),
  plot = amp_effs_vs_rra_18s,
  dpi = 600,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/amp_effs_vs_rra_18s_nocollodaria.pdf"),
  plot = amp_effs_vs_rra_18s,
  dpi = 600,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)





# #Figure S5 ASV Level ----------------------------------------------------

#Taxa file for matching Hash
taxa_18s=read.csv(here("PCR_bias_correction/data/taxa_files/blast_metazoo_18s.csv")) %>% 
  select(-X) %>% 
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>% 
  #Replace Orders that are empty with 'other'
  mutate_all(~replace_na(., "other")) %>% 
  distinct() %>%
  mutate(Species = case_when(
    Species != "other" & Species != "" ~ Species,
    Genus != "other" & Genus != "" ~ paste("unidentified", Genus),
    Family != "other" & Family != "" ~ paste("unidentified", Family),
    Order != "other" & Order != "" ~ paste("unidentified", Order),
    Class != "other" & Class != "" ~ paste("unidentified", Class),
    Phylum != "other" & Phylum != "" ~ paste("unidentified", Phylum),
    TRUE ~ "unidentified"
  ))


# Load and merge amplification efficiency data
all_amp_effs_18s_asv <- read.csv(here("PCR_bias_correction/data/amp_effs/amp_effs_18s_asv.csv"), row.names = 1) %>%
  mutate(Hash = str_extract(Lambda.coord, "[^_]+$")) %>%
  left_join(taxa_18s %>% select(Hash, Species), by = "Hash") %>%
  mutate(Species = replace_na(Species, "other"),
         Species_simple = str_remove(Species, "\\sASV\\s\\d+$"))

# Determine species order by mean Lambda.mean
ordered_species <- all_amp_effs_18s_asv %>%
  group_by(Species_simple) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_lambda) %>%
  pull(Species_simple)

# Determine hash order within species by mean Lambda.mean
hash_order_df <- all_amp_effs_18s_asv %>%
  group_by(Species_simple, Hash) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species_simple) %>%
  arrange(mean_lambda, .by_group = TRUE) %>%
  mutate(hash_order = row_number())

# Join hash order and define group label
all_amp_effs_18s_asv <- all_amp_effs_18s_asv %>%
  left_join(hash_order_df, by = c("Species_simple", "Hash")) %>%
  mutate(
    size_fraction = factor(size_fraction, levels = c(0.2, 0.5, 1)),
    Species_simple = factor(Species_simple, levels = ordered_species),
    group = paste(Species_simple, sprintf("h%02d", hash_order), size_fraction, sep = "_"),
    x_label = ifelse(size_fraction == 0.5, as.character(Species_simple), "")
  )

# Create ordered group factor
ordered_groups <- all_amp_effs_18s_asv %>%
  arrange(Species_simple, hash_order, size_fraction) %>%
  pull(group) %>%
  unique()

all_amp_effs_18s_asv <- all_amp_effs_18s_asv %>%
  mutate(group = factor(group, levels = ordered_groups))%>%
  mutate(group_chr = as.character(group))

# Define y-axis limits and breaks
lambda_range <- range(all_amp_effs_18s_asv$Lambda.mean, na.rm = TRUE)
max_abs <- max(abs(lambda_range)) * 1.5
y_limits <- c(-max_abs, max_abs)
y_breaks <- seq(from = floor(y_limits[1] * 10) / 10,
                to   = ceiling(y_limits[2] * 10) / 10,
                by = 0.1)

species_colors <- c(
  "Calanus helgolandicus" = "#77DD77",
  "Metridia pacifica" = "#9370DB",
  "unidentified Collodaria" = "#912330",
  "unidentified Euphausiidae" = "#FF6961",
  "unidentified Calanus" = "#2d8087",
  "other" = "grey",
  "Paracalanus parvus" = "#eb6098"
)

x_labels_named <- all_amp_effs_18s_asv %>%
  mutate(group_chr = as.character(group)) %>%
  group_by(Species_simple) %>%
  arrange(hash_order, size_fraction) %>%
  slice(floor(n() / 2) + 1) %>%
  ungroup() %>%
  select(group_chr, Species_simple) %>%
  mutate(
    group_chr = as.character(group_chr),
    Species_simple = as.character(Species_simple)
  ) %>%
  deframe()

# Create final plot
amp_effs_all_and_subpools_by_taxa_18s_asv <- ggplot(all_amp_effs_18s_asv, aes(x = group, y = Lambda.mean,
                                                                              color = Species,
                                                                              shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5), width = 0.2, size = 1.1) +
  geom_point(size = 4, stroke = 1.2) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_shape_manual(name = "Size Fraction",
                     values = c("0.2" = 16, "0.5" = 17, "1" = 18),
                     labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(x_labels_named), x_labels_named[x], ""))+
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::number_format(accuracy = 0.1)
  )+
  labs(
    x = "",
    y = "Amplification Efficiency (CLR)"
  ) +
  theme_minimal() +
  guides(color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  )

# View the plot
amp_effs_all_and_subpools_by_taxa_18s_asv



# Save the plot as PNG
output_dir <- "PCR_bias_correction/figures/supporting_info/"
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_asv.png"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_asv,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_asv.pdf"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_asv,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)

# Loosened Filtering Criteria ASV---------------------------------------------
#Taxa file for matching Hash
taxa_18s=read.csv(here("PCR_bias_correction/data/taxa_files/blast_metazoo_18s.csv")) %>% 
  select(-X) %>% 
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>% 
  #Replace Orders that are empty with 'other'
  mutate_all(~replace_na(., "other")) %>% 
  distinct() %>%
  mutate(Species = case_when(
    Species != "other" & Species != "" ~ Species,
    Genus != "other" & Genus != "" ~ paste("unidentified", Genus),
    Family != "other" & Family != "" ~ paste("unidentified", Family),
    Order != "other" & Order != "" ~ paste("unidentified", Order),
    Class != "other" & Class != "" ~ paste("unidentified", Class),
    Phylum != "other" & Phylum != "" ~ paste("unidentified", Phylum),
    TRUE ~ "unidentified"
  ))


# Load and merge amplification efficiency data
all_amp_effs_18s_asv_nofilt <- read.csv(here("PCR_bias_correction/data/amp_effs/amp_effs_18s_asv_nofilt.csv"), row.names = 1) %>%
  mutate(Hash = str_extract(Lambda.coord, "[^_]+$")) %>%
  left_join(taxa_18s %>% select(Hash, Species), by = "Hash") %>%
  mutate(Species = replace_na(Species, "other"),
         Species_simple = str_remove(Species, "\\sASV\\s\\d+$"))

# Determine species order by mean Lambda.mean
ordered_species <- all_amp_effs_18s_asv_nofilt %>%
  group_by(Species_simple) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_lambda) %>%
  pull(Species_simple)

# Determine hash order within species by mean Lambda.mean
hash_order_df <- all_amp_effs_18s_asv_nofilt %>%
  group_by(Species_simple, Hash) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species_simple) %>%
  arrange(mean_lambda, .by_group = TRUE) %>%
  mutate(hash_order = row_number())

# Join hash order and define group label
all_amp_effs_18s_asv_nofilt <- all_amp_effs_18s_asv_nofilt %>%
  left_join(hash_order_df, by = c("Species_simple", "Hash")) %>%
  mutate(
    size_fraction = factor(size_fraction, levels = c(0.2, 0.5, 1)),
    Species_simple = factor(Species_simple, levels = ordered_species),
    group = paste(Species_simple, sprintf("h%02d", hash_order), size_fraction, sep = "_"),
    x_label = ifelse(size_fraction == 0.5, as.character(Species_simple), "")
  )

# Create ordered group factor
ordered_groups <- all_amp_effs_18s_asv_nofilt %>%
  arrange(Species_simple, hash_order, size_fraction) %>%
  pull(group) %>%
  unique()

all_amp_effs_18s_asv_nofilt <- all_amp_effs_18s_asv_nofilt %>%
  mutate(group = factor(group, levels = ordered_groups))%>%
  mutate(group_chr = as.character(group))

# Define y-axis limits and breaks
lambda_range <- range(all_amp_effs_18s_asv_nofilt$Lambda.mean, na.rm = TRUE)
max_abs <- max(abs(lambda_range)) * 1.5
y_limits <- c(-max_abs, max_abs)
y_breaks <- seq(from = floor(y_limits[1] * 10) / 10,
                to   = ceiling(y_limits[2] * 10) / 10,
                by = 0.1)

species_colors <- c(
  "Calanus helgolandicus" = "#77DD77",
  "Metridia pacifica" = "#9370DB",
  "unidentified Collodaria" = "#912330",
  "unidentified Euphausiidae" = "#FF6961",
  "unidentified Calanus" = "#2d8087",
  "other" = "grey",
  "Paracalanus parvus" = "#eb6098"
)

x_labels_named <- all_amp_effs_18s_asv_nofilt %>%
  mutate(group_chr = as.character(group)) %>%
  group_by(Species_simple) %>%
  arrange(hash_order, size_fraction) %>%
  slice(floor(n() / 2) + 1) %>%
  ungroup() %>%
  select(group_chr, Species_simple) %>%
  mutate(
    group_chr = as.character(group_chr),
    Species_simple = as.character(Species_simple)
  ) %>%
  deframe()

# Add padding to the levels of the factor
group_levels_padded <- c("pad_left", levels(all_amp_effs_18s_asv_nofilt$group), "pad_right")

# Apply the new levels
all_amp_effs_18s_asv_nofilt <- all_amp_effs_18s_asv_nofilt %>%
  mutate(group = factor(group, levels = group_levels_padded))


# Create final plot
amp_effs_all_and_subpools_by_taxa_18s_asv_nofilt <- ggplot(all_amp_effs_18s_asv_nofilt, aes(x = group, y = Lambda.mean,
                                                                                            color = Species,
                                                                                            shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5), width = 0., size = 0.8) +
  geom_point(size = 2, stroke = 1) +
  scale_shape_manual(name = "Size Fraction",
                     values = c("0.2" = 16, "0.5" = 17, "1" = 18),
                     labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_x_discrete(labels = NULL)+
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    x = "",
    y = "Amplification Efficiency (CLR)"
  ) +
  theme_minimal() +
  # guides(color = "none") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    color = guide_legend(ncol = 4, title.position = "top", title.hjust = 0.5),
    shape = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)
  )+
  facet_wrap(~size_fraction, ncol=1)


# View the plot
amp_effs_all_and_subpools_by_taxa_18s_asv_nofilt



# Save the plot as PNG
output_dir <- "PCR_bias_correction/figures/supporting_info/"
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_asv_nofilt.png"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_asv_nofilt,
  dpi = 600,
  width = 7,
  height = 6,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_s4_aes_asv_nofilt.pdf"),
  plot = amp_effs_all_and_subpools_by_taxa_18s_asv_nofilt,
  dpi = 600,
  width = 7,
  height = 6,
  units = "in",
)


# Use a violin plot
ggplot(all_amp_effs_18s_asv_nofilt, aes(x = Species, y = Lambda.mean,
                                        color = Species,
                                        fill = Species,
                                        shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_violin(width = 1, alpha = 0.4, scale = "width", color = NA) +
  geom_point(
    size = 2, stroke = 1,
    position = position_jitter(width = 0.2, height = 0)
  ) +
  scale_shape_manual(
    name = "Size Fraction",
    values = c("0.2" = 16, "0.5" = 17, "1" = 18),
    labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")
  ) +

  scale_x_discrete(labels = NULL) +
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::number_format(accuracy = 0.1)
  ) +
  labs(
    x = NULL,
    y = "Amplification Efficiency (CLR)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    color = guide_legend(ncol = 4, title.position = "top", title.hjust = 0.5),
    fill = "none",
    shape = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)
  )


# AE vs Reads for ASVs ----------------------------------------------------


fido_18s_s1 = read.csv(here("PCR_bias_correction/data/fido/asv_level/fido_18s_s1_asv_nofilt.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_18s_s2 = read.csv(here("PCR_bias_correction/data/fido/asv_level/fido_18s_s2_asv_nofilt.csv"), header = TRUE, check.names = FALSE, row.names = 1)
fido_18s_s3 = read.csv(here("PCR_bias_correction/data/fido/asv_level/fido_18s_s3_asv_nofilt.csv"), header = TRUE, check.names = FALSE, row.names = 1)


# Function to sum columns in fido datasets containing "S1", "S2", or "S3" in their names
mean_s_columns_18s <- function(df, hash, suffix) {
  matching_columns <- names(df)[str_detect(names(df), suffix)]
  if (length(matching_columns) > 0) {
    normalized_values <- sapply(matching_columns, function(col) {
      column_sum_hash <- sum(df[df$Hash == hash, col], na.rm = TRUE)
      column_sum <- sum(df[[col]], na.rm = TRUE)
      if (column_sum_hash == 0 || column_sum == 0) {
        return(0)
      }
      column_sum_hash / column_sum
    })
    return(mean(normalized_values, na.rm = TRUE))
  } else {
    return(NA)
  }
}




# 18s ---------------------------------------------------------------------


# Create the final dataframe for S1, S2, and S3
result_s1 <- all_amp_effs_18s_asv_nofilt %>%
  filter(str_detect(pool, "S1")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s(fido_18s_s1, Hash, "S1"),
         SizeFraction = "S1") %>%
  select(Hash, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s2 <- all_amp_effs_18s_asv_nofilt %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s(fido_18s_s2, Hash, "S2"),
         SizeFraction = "S2") %>%
  select(Hash, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s3 <- all_amp_effs_18s_asv_nofilt %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_18s(fido_18s_s3, Hash, "S3"),
         SizeFraction = "S3") %>%
  select(Hash, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

# Combine results into one dataframe
result_combined <- bind_rows(result_s1, result_s2, result_s3)

# Summarize the data and order Hash by increasing Lambda.mean
hash_order <- result_combined %>%
  group_by(Hash) %>%
  summarise(mean_Lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_Lambda) %>%
  pull(Hash)

summarized_result <- result_combined %>%
  group_by(Hash, SizeFraction) %>%
  summarise(Lambda.mean = mean(Lambda.mean, na.rm = TRUE), Mean= mean(Mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Hash = factor(Hash, levels = hash_order))

pearson_18s <- cor.test(~ Lambda.mean + Mean, data = summarized_result)

cat("18S Dataset:\n")
cat("Pearson's R:", pearson_18s$estimate, "\n")
cat("p-value:", pearson_18s$p.value, "\n\n")


# Filter non-zero rows, rename NA
filtered_data <- summarized_result %>% filter(Mean > 0) %>% 
  left_join(., taxa_18s %>% select(Hash, Order), by="Hash")%>%
  mutate(Order = replace_na(Order, "unidentified Order"))  %>%
  mutate(Order = factor(Order, levels = names(taxa_colors_18s_order)))

# Fit linear models by SizeFraction and extract R² and p-value
lm_stats_by_size <- filtered_data %>%
  group_by(SizeFraction) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(log2(Mean) ~ Lambda.mean, data = .x)),
    model_summary = map(model, glance)
  ) %>%
  unnest(model_summary) %>%
  select(SizeFraction, r.squared, p.value)

# Format the text for annotation
lm_stats_by_size <- lm_stats_by_size %>%
  mutate(
    annotation = paste0("R² = ", round(r.squared, 2), 
                        "\nP = ", format.pval(p.value, digits = 2, eps = 0.001))
  )


taxa_colors_18s_order <- c(
  "Calanoida"             = "#89CFF0",  # Soft red-pink
  "Cyclopoida"            = "#00CED1",  # Light turquoise
  "Euphausiacea"          = "#FF6961",  # Baby blue
  "Salpida"               = "#FFB6C1",  # Pastel pink
  "Collodaria"            = "#912330",  # Light brown/gold
  "Doliolida"             = "#9ACD32",  # Yellow-green
  "Siphonophorae"         = "#CDA4DE",  # Medium orchid
  "unidentified Order"    = "#A9A9A9",  # Dark grey
  "other"                 = "#D3D3D3"   # Light grey
)




amp_effs_vs_rra_18s_asv <- filtered_data %>%
  ggplot(aes(x = Lambda.mean, y = log2(Mean))) +
  geom_point(
    aes(fill = Order, shape = SizeFraction),
    color = "black",
    size = 4,
    stroke = 0.8,
    alpha = 0.8
  ) +
  geom_text(
    data = lm_stats_by_size,
    aes(x = -Inf, y = Inf, label = annotation),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 4, fontface = "italic"
  ) +
  theme_classic() +
  labs(
    x = "Amplification Efficiency",
    y = "log2(mean relative read abundance)",
    fill = "Order",  # <- update to 'fill' instead of 'color' since you're using fill
    shape = "Size Fraction"
  ) +
  scale_fill_manual(values = taxa_colors_18s_order, name = "Order") +
  scale_shape_manual(
    values = c(21, 22, 23),
    labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"),
    name = "Size Fraction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    shape = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5),
    fill = guide_legend(
      override.aes = list(
        fill = unname(taxa_colors_18s_order),  # strip names from vector
        color = "black",
        shape = 21,
        stroke = 0.8,
        size = 4
      ),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 3
    )
  ) +
  facet_wrap(~SizeFraction, nrow = 1)

 amp_effs_vs_rra_18s_asv  
    

#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/figures_amp_effs_vs_rra_18s_asv.png"),
  plot = amp_effs_vs_rra_18s_asv,
  dpi = 600,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 4, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/supporting_info/figures_amp_effs_vs_rra_18s_nocollodaria.pdf"),
  plot = amp_effs_vs_rra_18s_asv,
  dpi = 600,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 4, # Adjust height (inches) based on journal requirements
  units = "in"
)


# Log2 Bias from Silverman ------------------------------------------------

#Plot Bias as the log ratio of the proportions
# ---- 1️⃣ Calculate log2 Fold Change and Error Bars ----
# Step 1: Calculate log2 bias and CI bounds
bias_df <- pcr_and_raw_18s_all %>%
  mutate(
    log2_bias = log2(n_reads_raw / n_reads_pcr),
    log2_low  = log2(n_reads_raw / p.97.5),
    log2_high = log2(n_reads_raw / p.2.5)
  ) %>%
  filter(is.finite(log2_bias), is.finite(log2_low), is.finite(log2_high)) %>% 
  left_join(metadata %>% select(Sample_ID, offshore_onshore), by="Sample_ID")

# Step 2: Summarize by taxa, size fraction, and cycle
bias_summary <- bias_df %>%
  group_by(taxa, size_fraction_numeric, offshore_onshore) %>%
  summarise(
    log2_bias = mean(log2_bias, na.rm = TRUE),
    log2_low = mean(log2_low, na.rm = TRUE),
    log2_high = mean(log2_high, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Order taxa by average bias across all groups
ordered_taxa <- bias_summary %>%
  group_by(taxa) %>%
  summarise(mean_bias = mean(log2_bias, na.rm = TRUE)) %>%
  arrange(mean_bias) %>%
  pull(taxa)

bias_summary <- bias_summary %>%
  mutate(taxa = factor(taxa, levels = ordered_taxa))

# Step 4: Plot
jitter_pos <- position_jitter(width = 0.15, height = 0)
ggplot(bias_summary, aes(x = taxa, y = log2_bias, color = taxa, shape = as.factor(size_fraction_numeric))) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
  geom_errorbar(aes(ymin = log2_low, ymax = log2_high),
                width = 0.6, linewidth = 0.9,
                position = jitter_pos) +
  geom_point(size = 3, stroke = 1, position = jitter_pos) +
  facet_wrap(~offshore_onshore, nrow = 2, labeller = label_both) +
  scale_shape_manual(name = "Size Fraction (mm)", values = c(16, 17, 18, 15)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)) +
  labs(
    x = "Taxa (ordered by mean PCR bias)",
    y = expression(log[2]~"(Raw / PCR Proportion)"),
    title = "Mean PCR Bias by Taxon, Size Fraction, and Cycle",
    color = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )





# COI Analysis --------------------------------------------------------------------


# COI ---------------------------------------------------------------------
plot_data_coi=pcr_and_raw_coi_all %>%
  # Create a new variable Sample_ID_short (everything before the first '_')
  mutate(Sample_ID_short = sub("_.*", "", Sample_ID)) %>%
  # Ensure Sample_ID is ordered by increasing PC1
  mutate(Sample_ID_short = reorder(Sample_ID_short, PC1)) %>%
  group_by(Sample_ID_short, taxa, size_fraction_numeric) %>%
  # Summarize raw and PCR reads
  summarize(
    n_reads_raw = sum(n_reads_raw, na.rm = TRUE),
    n_reads_pcr = sum(n_reads_pcr, na.rm = TRUE),
    PC1=mean(PC1),
    .groups = "drop"
  ) %>%
  # Pivot longer for faceting
  pivot_longer(
    cols = c(n_reads_raw, n_reads_pcr),
    names_to = "Metric",  # Create a column to identify raw and PCR
    values_to = "Reads"   # Move the values into a single column
  ) %>%
  group_by(Sample_ID_short, Metric, size_fraction_numeric) %>%
  # Normalize reads within each Sample_ID for both metrics
  mutate(Proportion = Reads / sum(Reads, na.rm = TRUE))


# Create custom labels for size_fraction_numeric and Metric
facet_labels <- list(
  size_fraction_numeric = c(
    "0.2-0.5mm" = "0.5",
    "0.5-1mm" = "1",
    "1-2mm" = "2"
  ),
  Metric = c(
    "n_reads_pcr" = "PCR corrected",
    "n_reads_raw" = "Raw"
  )
)

# Define labeller function
custom_labeller <- labeller(
  size_fraction_numeric = c(
    "0.2" = "0.2-0.5mm",
    "0.5" = "0.5-1mm",
    "1" = "1-2mm"
  ),
  Metric = c(
    "n_reads_pcr" = "PCR corrected",
    "n_reads_raw" = "Raw"
  )
)




# PC1 Labels
labels_for_PC1 <- metadata %>% 
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange(PC1)

#Colors
# Define the new taxa color palette
taxa_colors_coi <- c(
  "Euphausia" = "#FF6961",          # Pastel red (from Euphausiidae)
  "Ctenocalanus" = "#77DD77",       # Pastel green (from Calanidae)
  "Calanus" = "#77DD77",            # Pastel green (from Calanidae)
  "Eucalanus" = "#89CFF0",          # Baby blue (from Eucalanidae)
  "Clausocalanus" = "#72872d",      # Bright pastel green (from Clausocalanidae)
  "Pleuromamma" = "#9370DB",        # Pastel purple (from Metridinidae)
  "Rosacea" = "#912330",            # Pastel orange (from unidentified Collodaria)
  "Sagitta" = "#CDA4DE",            # Pastel lavender (from unidentified Siphonophorae)
  "Metridia" = "#9370DB",           # Pastel purple (from Metridinidae)
  "Acrocalanus" = "#eb6098",        # Lime pastel green (from Paracalanidae)
  "unidentified Halocyprididae" = "#2d8087",  # Matching unidentified Calanoida
  "Acartia" = "#f0ca62",            # Light blue (from Oithonidae)
  "Ditrichocorycaeus" = "#FFB6C1"   # Pastel pink (from Salpidae)
)

fig1 <- ggplot(plot_data_coi, aes(x = Sample_ID_short, y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) + # Apply custom labels
  scale_fill_manual(values = taxa_colors_coi) + # Apply custom taxa colors
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)+
  labs(
    title = "",
    x = expression(paste("Offshore ", PC1, " Onshore")),
    y = "Proportion Reads",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
    axis.title.x = element_text(size = 16),  # Larger axis title (O)
    axis.title.y = element_text(size = 16),  # Larger y-axis title (O)
    plot.title = element_text(size = 18, face = "bold"),  # Increased title size (O)
    strip.text = element_text(face = "bold", size = 16), # Facet text size
    legend.text = element_text(size = 16),  # Larger legend text (L)
    legend.title = element_text(size = 16, face = "bold"),  # Larger, bold legend title (L)
    legend.position = "bottom",
    panel.spacing = unit(0.8, "lines"))  # Reduce spacing for compact fit)  # Show legend only on first plot

fig1

# Define output directory
output_dir <- "PCR_bias_correction/figures/v0"

# Save the plot as PNG
ggsave(
  filename = file.path(output_dir, "figure_2_v0.png"),
  plot = fig1,
  dpi = 600,
  width = 18,
  height = 7,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_2_v0.pdf"),
  plot = fig1,
  dpi = 600,
  width = 18,
  height = 7,
  units = "in",
)


## By Order to better visualize alognside zooscan
# Define a unique color palette for each Order
order_colors <- c(
  "other" = "grey",
  "Calanoida" = "#8DD3C7",       # Pastel green
  "Cyclopoida" = "#89CFF0",      # Baby blue
  "Euphausiacea" = "#BC80BD",    # Pastel red
  "Collodaria" = "#2d8087"      # Deep blue-green
)

# Modify plot_data to mutate `taxa` into `Genus`, then join with `zhan_taxa` by Genus
plot_data_coi <- plot_data_coi %>%
  rename(Genus = taxa) %>%  # Rename taxa to Genus
  left_join(coi_taxa %>% rownames_to_column("Genus"), by = "Genus") %>%  # Join with zhan_taxa
  mutate(Order = ifelse(is.na(Order), "other", Order))  # Assign "other" for missing values

# Create the stacked bar plot by Order
fig1_order <- ggplot(plot_data_coi, aes(x = Sample_ID_short, y = Proportion, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) +  # Apply custom labels
  # scale_fill_manual(values = order_colors) +  # Apply custom Order colors
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short) +
  labs(
    title = "",
    x = "Sample ID",
    y = "Proportion Reads",
    fill = "Order"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increased for better readability (O)
    axis.title.x = element_text(size = 16),  # Larger axis title (O)
    axis.title.y = element_text(size = 16),  # Larger y-axis title (O)
    plot.title = element_text(size = 18, face = "bold"),  # Increased title size (O)
    legend.text = element_text(size = 16),  # Larger legend text (L)
    strip.text = element_text(face = "bold", size = 16), # Facet text size
    legend.title = element_text(size = 16, face = "bold"),  # Larger, bold legend title (L)
    legend.position ="bottom")  # Show legend only on first plot

# Print the plot
fig1_order

# Define output directory
output_dir <- "PCR_bias_correction/figures/v0"

# Save the plot as PNG
ggsave(
  filename = file.path(output_dir, "figure_2_v0_order.png"),
  plot = fig1_order,
  dpi = 600,
  width = 18,
  height = 7,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_2_v0_order.pdf"),
  plot = fig1_order,
  dpi = 600,
  width = 18,
  height = 7,
  units = "in",
)


# Define the new taxa color palette
taxa_colors_coi <- c(
  "other" = "grey",
  "Euphausia" = "#FF6961",          # Pastel red
  "Ctenocalanus" = "#4682B4",       # Steel blue
  "Calanus" = "#0e9c38",            # Pastel green
  "Eucalanus" = "#89CFF0",          # Baby blue
  "Clausocalanus" = "#72872d",      # Olive green
  "Pleuromamma" = "#628be3",        # Light periwinkle
  "Rosacea" = "#CDA4DE",            # Lavender
  "Sagitta" = "#89f57f",            # Mint green
  "Metridia" = "#9370DB",           # Purple
  "Acrocalanus" = "#eb6098",        # Pink-lilac
  "unidentified Halocyprididae" = "#f5971d",  # Orange-gold
  "Acartia" = "#f0ca62",            # Pale yellow
  "Ditrichocorycaeus" = "#FFB6C1",  # Light pink
  "Membranipora" = "#8B008B",       # Dark magenta (distinct for bryozoan)
  "Hematodinium" = "#DC143C",       # Crimson red (flagged parasitic)
  "Paraeuchaeta" = "#20B2AA",        # Light sea green (copepod, carnivore)
  "unidentified" = "#7c8087"
)


# COI: AE Analysis --------------------------------------------------------

#Use all taxa
# Load data
all_amp_effs_coi <- read.csv(here("PCR_bias_correction/data/amp_effs/all_amp_effs_coi_all_sub_nofilt.csv")) %>%
  mutate(Genus = str_extract(Lambda.coord, "[^_]+$"))

# Calculate mean Lambda per Genus and order them ascending
ordered_taxa <- all_amp_effs_coi %>%
  group_by(Genus) %>%
  summarise(mean_lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_lambda) %>%
  pull(Genus)

# Determine which genera include a 0.5 size fraction
has_05 <- all_amp_effs_coi %>%
  group_by(Genus) %>%
  summarise(has_0_5 = any(size_fraction == 0.5), .groups = "drop")

# Join back and apply conditional labeling
all_amp_effs_coi <- all_amp_effs_coi %>%
  left_join(has_05, by = "Genus") %>%
  mutate(
    Genus = factor(Genus, levels = ordered_taxa),
    size_fraction = factor(size_fraction, levels = c(0.2, 0.5, 1)),
    group = interaction(Genus, size_fraction, lex.order = TRUE),
    x_label = case_when(
      has_0_5 ~ ifelse(size_fraction == 0.5, as.character(Genus), " "),
      !has_0_5 ~ ifelse(size_fraction == 1, as.character(Genus), " "),
      TRUE ~ " "  # fallback
    )
  )
# Create named vector of labels
x_labels_named <- all_amp_effs_coi %>%
  distinct(group, x_label) %>%
  deframe()
x_labels_named

# Order 'group' by Genus order and within-Genus size fraction
ordered_groups <- all_amp_effs_coi %>%
  arrange(Genus, size_fraction) %>%
  pull(group) %>%
  unique()

# Apply ordered factor
all_amp_effs_coi <- all_amp_effs_coi %>%
  mutate(group = factor(group, levels = ordered_groups))

# Compute dynamic limits and breaks
lambda_range <- range(all_amp_effs_coi$Lambda.mean, na.rm = TRUE)
max_abs <- max(abs(lambda_range)) * 1.5
y_limits <- c(-max_abs, max_abs)
y_breaks <- seq(from = floor(y_limits[1] * 10) / 10,
                to   = ceiling(y_limits[2] * 10) / 10,
                by = 0.1)

# Final plot
amp_effs_all_and_subpools_by_taxa_coi_no_filt <- ggplot(all_amp_effs_coi, aes(x = group, y = Lambda.mean,
                                                                              color = Genus,
                                                                              shape = as.factor(size_fraction))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5), width = 0.2, size = 1.1) +
  geom_point(size = 4, stroke = 1.2) +
  scale_color_manual(values = taxa_colors_coi) +
  scale_shape_manual(name = "Size Fraction",
                     values = c("0.2" = 16, "0.5" = 17, "1" = 18),
                     labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  scale_x_discrete(labels = x_labels_named)+
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 0.1, trim = TRUE)
  ) +
  labs(
    x = "",
    y = "Amplification Efficiency (CLR)",
    title = "",
    color = "Genus"
  ) +
  theme_minimal() +
  guides(color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(1, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.01, "lines"),
    legend.spacing.y = unit(0.01, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  )



# View the plot
amp_effs_all_and_subpools_by_taxa_coi_no_filt



# Save the plot as PNG
output_dir <- "PCR_bias_correction/figures/supporting_info/coi/"
ggsave(
  filename = file.path(output_dir, "figure_s_aes_nofilt_coi.png"),
  plot = amp_effs_all_and_subpools_by_taxa_coi_no_filt,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_s_aes_nofilt_coi.pdf"),
  plot = amp_effs_all_and_subpools_by_taxa_coi_no_filt,
  dpi = 600,
  width = 7,
  height = 5.5,
  units = "in",
)


detach("package:fido", unload = TRUE)

# Function to sum columns in fido datasets containing "S1", "S2", or "S3" in their names
mean_s_columns_coi_genus <- function(df, genus, suffix) {
  matching_columns <- names(df)[str_detect(names(df), suffix)]
  if (length(matching_columns) > 0) {
    normalized_values <- sapply(matching_columns, function(col) {
      column_sum_hash <- sum(df[df$Genus == genus, col], na.rm = TRUE)
      column_sum <- sum(df[[col]], na.rm = TRUE)
      if (column_sum_hash == 0 || column_sum == 0) {
        return(0)
      }
      column_sum_hash / column_sum
    })
    return(mean(normalized_values, na.rm = TRUE))
  } else {
    return(NA)
  }
}

all_amp_effs_coi_nofilt <- read.csv(here("PCR_bias_correction/data/amp_effs/all_amp_effs_coi_all_sub_nofilt.csv"), row.names = 1) %>%
  mutate(Genus = str_extract(Lambda.coord, "[^_]+$"))



fido_coi_s1_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s1_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)
fido_coi_s2_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s2_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)
fido_coi_s3_nofilt = read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s3_genus_phy_all_subpools.csv"), header = TRUE, check.names = FALSE)


# coi ---------------------------------------------------------------------


# Create the final dataframe for S1, S2, and S3
result_s1_nofilt <- all_amp_effs_coi_nofilt %>%
  filter(str_detect(pool, "S1")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_coi_genus(fido_coi_s1_nofilt, Genus, "S1"),
         SizeFraction = "S1") %>%
  select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s2_nofilt <- all_amp_effs_coi_nofilt %>%
  filter(str_detect(pool, "S2")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_coi_genus(fido_coi_s2_nofilt, Genus, "S2"),
         SizeFraction = "S2") %>%
  select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

result_s3_nofilt <- all_amp_effs_coi_nofilt %>%
  filter(str_detect(pool, "S3")) %>%
  rowwise() %>%
  mutate(Mean= mean_s_columns_coi_genus(fido_coi_s3_nofilt, Genus, "S3"),
         SizeFraction = "S3") %>%
  select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
  ungroup()

# Combine S1, S2, S3 with updated mean calc and SizeFraction tag
result_combined_nofilt <- bind_rows(
  all_amp_effs_coi_nofilt %>% filter(str_detect(pool, "S1")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_coi_genus(fido_coi_s1_nofilt, Genus, "S1"),
           SizeFraction = "S1") %>%
    select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup(),
  all_amp_effs_coi_nofilt %>% filter(str_detect(pool, "S2")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_coi_genus(fido_coi_s2_nofilt, Genus, "S2"),
           SizeFraction = "S2") %>%
    select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup(),
  all_amp_effs_coi_nofilt %>% filter(str_detect(pool, "S3")) %>%
    rowwise() %>%
    mutate(Mean = mean_s_columns_coi_genus(fido_coi_s3_nofilt, Genus, "S3"),
           SizeFraction = "S3") %>%
    select(Genus, pool, Lambda.mean, Mean, SizeFraction) %>%
    ungroup()
)

# Reorder Genus by increasing Lambda.mean
genus_order <- result_combined_nofilt %>%
  group_by(Genus) %>%
  summarise(mean_Lambda = mean(Lambda.mean, na.rm = TRUE)) %>%
  arrange(mean_Lambda) %>%
  pull(Genus)


# Compute correlation stats per size fraction
lm_stats_by_size <- summarized_result %>%
  filter(Mean > 0) %>%
  group_by(SizeFraction) %>%
  summarise(
    fit = list(lm(log2(Mean) ~ Lambda.mean)),
    .groups = "drop"
  ) %>%
  mutate(
    r2 = map_dbl(fit, ~ summary(.x)$r.squared),
    p = map_dbl(fit, ~ summary(.x)$coefficients[2, 4]),
    annotation = paste0(
      "R² = ", formatC(r2, digits = 2, format = "f"), "\n",
      "P = ", ifelse(p < 0.001, "<0.001", formatC(p, digits = 2, format = "f"))
    )
  )

# Summarize and set factor order
summarized_result <- result_combined_nofilt %>%
  group_by(Genus, SizeFraction) %>%
  summarise(
    Lambda.mean = mean(Lambda.mean, na.rm = TRUE),
    Mean = mean(Mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Genus = factor(Genus, levels = genus_order))%>%
  mutate(Genus = factor(Genus, levels = names(taxa_colors_coi)))

# Correlation
pearson_coi <- cor.test(~ Lambda.mean + Mean, data = summarized_result)
cat("Pearson's R:", pearson_coi$estimate, "\n")
cat("p-value:", pearson_coi$p.value, "\n\n")


size_fraction_labels <- c(
  "S1" = "0.2–0.5 mm",
  "S2" = "0.5–1 mm",
  "S3" = "1–2 mm"
)

# Plot
amp_effs_vs_rra_coi_nofilt <- summarized_result %>%
  # filter(Mean > 0) %>%
  # filter(!is.na(cycle))
  ggplot(aes(x = Lambda.mean, y = log2(Mean))) +
  geom_point(
    aes(fill = Genus, shape = SizeFraction),
    size = 4,
    stroke = 0.8,
    color = "black",
    alpha = 0.85
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8)+
  geom_text(
    data = lm_stats_by_size,
    aes(x = -Inf, y = Inf, label = annotation),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 4, fontface = "italic"
  ) +
  scale_fill_manual(values = taxa_colors_coi) +
  scale_shape_manual(
    values = c(21, 22, 23),
    labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"),
    name = "Size Fraction"
  ) +
  labs(
    x = "Amplification Efficiency",
    y = expression(log[2]*"(mean relative read abundance)"),
    fill = "Order",  # <- update to 'fill' instead of 'color' since you're using fill
    shape = "Size Fraction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face="bold"),
    legend.key.size = unit(0.4, "lines"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1, "lines"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.spacing.x = unit(0.1, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90")
  ) +
  guides(
    shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5),
    fill = guide_legend(
      override.aes = list(
        fill = unname(taxa_colors_coi),  # strip names from vector
        color = "black",
        shape = 21,
        stroke = 0.8,
        size = 4
      ),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 4
    )
  ) +
  facet_wrap(~SizeFraction, nrow = 1,
             labeller = labeller(SizeFraction = size_fraction_labels))

amp_effs_vs_rra_coi_nofilt

#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/v0/figs_amp_effs_vs_rra_coi_nofilt.png"),
  plot = amp_effs_vs_rra_coi_nofilt,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/v0/figs_amp_effs_vs_rra_coi_nofilt.png"),
  plot = amp_effs_vs_rra_coi_nofilt,
  dpi = 600,
  width = 6, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)


#COmbine for figure 4
library(patchwork)

# Remove x-axis labels and ticks from the top plot
amp_effs_all_and_subpools_by_taxa_coi_no_filt <- amp_effs_all_and_subpools_by_taxa_coi_no_filt +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"  # this removes the legend entirely from the top plot
  )

amp_effs_vs_rra_coi_nofilt <- amp_effs_vs_rra_coi_nofilt +
  guides(
    shape = guide_legend(
      ncol = 1,  # this forces one column
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5
    ),
    fill = guide_legend(
      override.aes = list(
        fill = unname(taxa_colors_coi),
        color = "black",
        shape = 21,
        stroke = 0.8,
        size = 4
      ),
      title.position = "top",
      title.hjust = 0.5,
      nrow = 4
    )
  )


# Combine the plots
figure_scoi_combined <-  wrap_elements(amp_effs_all_and_subpools_by_taxa_coi_no_filt) /
  amp_effs_vs_rra_coi_nofilt +
  plot_layout(heights = c(1.1, 1), guides = "keep") +
  plot_annotation(tag_levels = 'A')  # This will label them "A", "B"

figure_scoi_combined

# Save the combined figure
ggsave(
  filename = file.path(output_dir, "figure_scoi_combined.png"),
  plot = figure_scoi_combined,
  dpi = 600,
  width = 7,
  height = 8,  # Adjust height to accommodate both plots vertically
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "figure_scoi_combined.pdf"),
  plot = figure_s4_combined,
  dpi = 600,
  width = 7,
  height = 8,
  units = "in"
)



# COI: Zooscan Comparison -------------------------------------------------


#Taxa file from pre-processed fido families for 18S
coi_taxa=read.csv(here("PCR_bias_correction/data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv"))  %>% 
  select(-X) %>% 
  distinct() %>% 
  filter(Family != "Cliidae") %>% 
  column_to_rownames("Genus")



#COI
# Load files dynamically for s1, s2, and s3
load_most_recent_file_coi <- function(suffix,primer) {
  # Directory path
  target_dir <- here("PCR_bias_correction/data/predicted_og")
  
  # List all files matching the pattern for the given suffix
  files <- list.files(target_dir, 
                      pattern = paste0("predicted_og_", primer, "_\\d{2}_\\d{2}_\\d{4}_", suffix, "\\_nofilt.csv"), 
                      full.names = TRUE)
  
  
  if (length(files) == 0) {
    stop(paste("No files found for suffix:", suffix))
  }
  
  # Extract dates from filenames
  extract_date <- function(file) {
    match <- regmatches(file, regexpr("\\d{2}_\\d{2}_\\d{4}", file))
    as.Date(match, format = "%m_%d_%Y")
  }
  
  # Find the most recent file
  files_with_dates <- data.frame(
    file = files,
    date = sapply(files, extract_date)
  )
  
  latest_file <- files_with_dates %>%
    arrange(desc(date)) %>%
    slice(1) %>%
    pull(file)
  
  # Read the file and apply transformations
  read.csv(latest_file) %>%
    select(-X) %>%
    mutate(coord = str_remove(coord, "^clr_"))
}


fido_s1_coi <- load_most_recent_file_coi("s1","coi")
fido_s2_coi <- load_most_recent_file_coi("s2","coi")
fido_s3_coi <- load_most_recent_file_coi("s3","coi")


#Merge
final_data_all_sizes_coi=rbind(fido_s1_coi,fido_s2_coi,fido_s3_coi) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 


#Make final dataframe
phy_taxa_pcr_coi= final_data_all_sizes_coi %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr_coi=phy_taxa_pcr_coi %>% mutate(Genus=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Genus"), by="Genus") %>%
  filter(Order=="Calanoida") 


#PCR-RA df ready to join with RRA
pcr_join_prop_coi=phy_taxa_pcr_coi %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction_numeric,PC1,cycle,taxa) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))


fido_s1_raw_coi=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s1_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)



fido_s2_raw_coi=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s2_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) 


fido_s3_raw_coi=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s3_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)




merge(fido_s1_raw_coi, fido_s2_raw_coi, by = "Genus", all = TRUE) %>%
  merge(.,fido_s3_raw_coi, by = "Genus", all = TRUE)%>%
  column_to_rownames("Genus") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_coi_merged_raw



#Metadata
env_metadata_phy=env_metadata %>%
  column_to_rownames("Sample_ID_dot")


OTU = otu_table(as.matrix(fido_coi_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(env_metadata_phy)

#Raw in counts
phy_coi_raw_counts=phyloseq(OTU, TAX, meta) %>% 
  phyloseq_transform_to_long(.) %>%
  mutate(Genus=asv_code) %>% 
  select(-asv_code)


#Raw in Proportions
phy_coi=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Genus=asv_code) %>%
  select(-asv_code)



#Filter to calanoid copepods
taxa_sel="Calanoida"

#2. Join with PCR Proportions
#All taxa
pcr_and_raw_coi_all=
  phy_coi %>% 
  mutate(taxa=Genus) %>% 
  group_by(Sample_ID,taxa) %>%
  summarise(n_reads_raw=sum(n_reads)) %>% 
  left_join(pcr_join_prop_coi, by=c("Sample_ID","taxa"))

#3. Join in counts
phy_coi_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join_prop_coi, by="Sample_ID") %>%
  select(-size_fraction_numeric.x) %>%
  rename(size_fraction_numeric=size_fraction_numeric.y)->pcr_and_raw_coi_counts

#Add PCR-bias corrected counts
counts_raw_all_coi=phy_coi_raw_counts %>% 
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw_sum=sum(n_reads)) %>% 
  select(Sample_ID,size_fraction_numeric,n_reads_raw_sum)


pcr_and_raw_coi_counts=pcr_and_raw_coi_counts %>% 
  left_join(.,counts_raw_all_coi,by=c("Sample_ID","size_fraction_numeric")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



pcr_and_raw_coi_calanoids=pcr_and_raw_coi_all %>% 
  mutate(Genus=taxa) %>% 
  left_join(.,coi_taxa %>% rownames_to_column("Genus") %>% 
              select(Genus, Order), by="Genus") %>% 
  filter(Order==taxa_sel) %>% 
  select(-taxa)%>% 
  filter(Genus=="Calanus") %>%
  group_by(Sample_ID,size_fraction_numeric,cycle) %>% 
  summarise(n_reads_raw=sum(n_reads_raw), 
            n_reads_pcr=sum(n_reads_pcr),
            PC1=mean(PC1)) 

#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Calanoida"


#Zooscan biomass proportion
# All taxa
zooscan_all=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass_esd.csv"), row.names = 1)%>%
  # select(-X, acq_min_mesh) %>% 
  # mutate(Sample_ID=sample_id) %>%  
  distinct(.)


#Calanoida
zooscan_taxa=zooscan_all%>%
  filter(object_annotation_category==taxa_sel)  



#Propotions
pcr_raw_zoo_coi=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(pcr_and_raw_coi_calanoids, by=c("PC1","size_fraction_numeric")) %>% 
  unique(.)





# Correlations. Make a grid plot of Y axis is:  ----------------------------
# Raw reads or PCR corrected
# X is: Zooscan biomass, abundance or relative of each

facet_x_labels <- c(
  "biomass_prop_taxa" = "Log10(Biomass Proportion)",
  "dryweight_C_mg_m2_taxa_mean" = "Log10 (Dry Weight (mg Cm²))",
  "relative_abundance" = "Log10(Relative Abundance)",
  "abundance_m2_mean" = "Log10(Abundance (m²))"
)

facet_y_labels <- c(
  "n_reads_raw" = "Raw Reads",
  "n_reads_pcr" = "PCR Bias-Corrected Reads"
)

size_fraction_labels <- c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")

# Convert data to long format and filter for log transformations
pcr_raw_zoo_coi_long <- pcr_raw_zoo_coi %>%
  filter(!is.na(cycle)) %>%
  pivot_longer(cols = c(biomass_prop_taxa, dryweight_C_mg_m2_taxa_mean, relative_abundance, abundance_m2_mean), 
               names_to = "x_variable", values_to = "x_value") %>%
  pivot_longer(cols = c(n_reads_raw, n_reads_pcr), 
               names_to = "y_variable", values_to = "y_value") %>%
  filter(x_value > 0, y_value >= 0) %>%  # Ensure valid sqrt transform (y_value >= 0)
  mutate(
    x_variable = factor(x_variable, levels = names(facet_x_labels), labels = facet_x_labels),
    y_variable = factor(y_variable, levels = names(facet_y_labels), labels = facet_y_labels),
    size_fraction = factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2"))
  )


# Total number of unique comparisons
n_tests <- pcr_raw_zoo_coi_long %>%
  distinct(size_fraction, y_variable, x_variable) %>%
  nrow()

# Compute correlation per facet combination
cor_results_grid <- pcr_raw_zoo_coi_long %>%
  group_by(size_fraction, y_variable, x_variable) %>%
  summarise(
    cor_test = list(cor.test(log10(x_value), asin(sqrt(y_value)))),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(cor_test, ~ .x$estimate),
    p_uncorrected = map_dbl(cor_test, ~ .x$p.value),
    p_bonferroni = pmin(p_uncorrected * n_tests, 1),
    label = paste0(
      "r = ", formatC(r, digits = 2, format = "f"), "\n",
      "P = ", ifelse(p_bonferroni < 0.001, "<0.001", formatC(p_bonferroni, digits = 2, format = "f"))
    )
  )

# Compute per-facet positions using quantiles
label_positions <- pcr_raw_zoo_coi_long %>%
  group_by(x_variable, y_variable, size_fraction) %>%
  summarise(
    x = quantile(log(x_value), 0.25, na.rm = TRUE),
    y = quantile(asin(sqrt(y_value)), 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Merge with your correlation label table
cor_results_grid_pos <- cor_results_grid %>%
  left_join(label_positions, by = c("x_variable", "y_variable", "size_fraction"))


# Create 6-row, 4-column facet grid plot
zoo_vs_pcr_grid <- ggplot(pcr_raw_zoo_coi_long, 
                          aes(x = log((x_value)), 
                              y = asin(sqrt(y_value)))) +  # Apply asin sqrt transformation
  geom_point(aes(shape = cycle, size = 3, color = as.factor(size_fraction), fill = as.factor(size_fraction))) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                    labels = size_fraction_labels) +
  scale_color_manual(values = c("#5BA3D5", "#66CC66", "#FF4C38"), 
                     labels = size_fraction_labels) +
  facet_grid(rows = vars(size_fraction, y_variable), cols = vars(x_variable), scales = "free") +  # 6x4 layout
  labs(
    x = "Log10 Transformed X-Axis Variables",
    y = "asin(sqrt(Read Counts))",
    shape = "Cycle",
    color = "Size Fraction"
  ) +
  geom_text(
    data = cor_results_grid_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "italic"
  )+  # Pearson correlation
  guides(size = FALSE, fill = FALSE) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Major grid lines for readability
    panel.grid.minor = element_line(color = "gray90", size = 0.2),  # Minor grid lines for subtle separation
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around each facet
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Light gray background for facet labels
    strip.text = element_text(size = 16, face = "bold"),  # Larger and bold facet labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


# Print plot
zoo_vs_pcr_grid


#Just Absolute biomass
# Filter and reshape for plotting
pcr_vs_dryweight <- pcr_raw_zoo_coi %>%
  filter(!is.na(cycle)) %>%
  select(sample_id, cycle, size_fraction, dryweight_C_mg_m2_taxa_mean, n_reads_raw, n_reads_pcr) %>%
  pivot_longer(cols = c(n_reads_raw, n_reads_pcr),
               names_to = "y_variable", values_to = "y_value") %>%
  filter(dryweight_C_mg_m2_taxa_mean > 0, y_value >= 0) %>%
  mutate(
    y_variable = factor(y_variable, levels = c("n_reads_raw", "n_reads_pcr"),
                        labels = c("Raw Reads", "PCR-Corrected Reads")),
    size_fraction = factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2"))
    
  )

# Number of comparisons = number of facet combinations
n_tests <- pcr_vs_dryweight %>%
  distinct(y_variable, size_fraction) %>%
  nrow()

# Compute Pearson correlation per facet
cor_results <- pcr_vs_dryweight %>%
  group_by(y_variable, size_fraction) %>%
  summarise(
    cor_test = list(cor.test(
      ~ log10(dryweight_C_mg_m2_taxa_mean) + asin(sqrt(y_value))
    )),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(cor_test, ~ .x$estimate),
    p_uncorrected = map_dbl(cor_test, ~ .x$p.value),
    p_bonferroni = pmin(p_uncorrected * n_tests, 1),  # Bonferroni correction
    label = paste0("r = ", formatC(r, digits = 2, format = "f"), 
                   "\nP = ", ifelse(p_bonferroni < 0.001, "<0.001", formatC(p_bonferroni, digits = 2, format = "f")))
  )

#Define pastel color palette
pastel_colors <- c("0.2-0.5" = "#AEDFF7", "0.5-1" = "#B6E3B6", "1-2" = "#FFB3AB")

# Plot
ggplot(pcr_vs_dryweight,
       aes(x = log10(dryweight_C_mg_m2_taxa_mean),
           y = asin(sqrt(y_value)),
           fill = size_fraction)) +
  geom_point(aes(shape = cycle), size = 6, stroke = 1.2, color = "black") +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
  scale_fill_manual(values = pastel_colors, labels = size_fraction_labels, name = "Size Fraction") +
  coord_cartesian(xlim = c(0, 3)) +
  facet_grid(rows = vars(y_variable), cols = vars(size_fraction), scales = "free_y") +
  labs(
    x = "Log10 Dry Weight C (mg)",
    y = "arcsine-sqrt Relative Read Abundance",
    shape = "Cycle",
    fill = "Size Fraction"
  ) +
  geom_text(
    data = cor_results,
    aes(label = label),
    x = 0.25, y = 0.5,
    inherit.aes = FALSE,
    size = 4, fontface = "italic"
  )  +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6, stroke = 1.2)),
    shape = guide_legend(override.aes = list(fill = "gray90", color = "black")),
    color = "none"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(size = 16, face = "bold", genus = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold", genus = "serif"),
    axis.text.y = element_text(size = 14, face = "bold", genus = "serif"),
    axis.title = element_text(size = 16, face = "bold", genus = "serif"),
    legend.title = element_text(size = 16, face = "bold", genus = "serif"),
    legend.text = element_text(size = 14, genus = "serif")
  )-> methods_correlation_plot

methods_correlation_plot
#PNG & PDF Save
ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_coi_scatter.png"),
  plot = methods_correlation_plot,
  dpi = 600,
  width = 12, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)

ggsave(
  filename = here("PCR_bias_correction/figures/methods_comparison/calanoid_methods_compare_coi_scatter.pdf"),
  plot = methods_correlation_plot,
  dpi = 600,
  width = 12, # Adjust width (inches) based on journal requirements
  height = 7, # Adjust height (inches) based on journal requirements
  units = "in"
)


# SCRAP -------------------------------------------------------------------


#Proportions
filtered_data_18s <- final_data_all_sizes_18s %>%
  filter(p2.5 > 0 | p97.5 < 0) %>%  # Keep only rows where both values are either >0 or <0
  mutate(Sample_ID = as.factor(Sample_ID),  # Ensure Sample_ID is categorical
         coord = as.factor(coord), 
         size = as.factor(size)) %>%   # Ensure size is categorical for shape mapping
  left_join(env_metadata %>% select(Sample_ID,Sample_ID_short,PC1), by="Sample_ID")

labels_for_map=filtered_data_18s %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

#CLR
filtered_data_clr_18s <- final_data_all_sizes_clr %>%
  filter(p2.5 > 0 | p97.5 < 0) %>%  # Keep only rows where both values are either >0 or <0
  mutate(Sample_ID = as.factor(Sample_ID),  # Ensure Sample_ID is categorical
         coord = as.factor(coord), 
         size = as.factor(size)) %>%   # Ensure size is categorical for shape mapping
  left_join(env_metadata %>% select(Sample_ID,Sample_ID_short,PC1), by="Sample_ID")

# Create the plot: CLR
ggplot(filtered_data_clr, aes(x = as.factor(PC1), y = n_reads, color = coord, shape = size)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Larger points
  geom_errorbar(aes(ymin = p2.5, ymax = p97.5), width = 0.3, position = position_dodge(width = 0.5), linewidth = 1) +  # Larger error bars
  labs(title = "Filtered Error Bar Plot of n_reads by Sample_ID",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Number of Reads",
       color = "Coord",
       shape = "Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_color_manual(values = taxa_colors)+  # Apply custom color palette
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)

# Create the plot: Proportions
ggplot(filtered_data_clr_18s %>% filter(size=="0.2-0.5mm"), aes(x = as.factor(PC1), y = n_reads, color = coord, shape = size)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Larger points
  geom_errorbar(aes(ymin = p2.5, ymax = p97.5), width = 0.3, position = position_dodge(width = 0.5), linewidth = 1) +  # Larger error bars
  labs(title = "Filtered Error Bar Plot of n_reads by Sample_ID",
       x = "Offfshore \u2190 PC1 \u2192 Onshore",
       y = "Number of Reads",
       color = "Coord",
       shape = "Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  # scale_color_manual(values = taxa_colors)+  # Apply custom color palette
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  facet_wrap(~size, nrow=3)



#Vis 1. View proportion taxa in each sample
plot_data <- pcr_join_prop_18s %>%
  group_by(Sample_ID, taxa, size_fraction_numeric, PC1) %>%
  # Calculate the sum of n_reads_pcr for each taxa within each Sample_ID
  summarize(n_reads_pcr = sum(n_reads_pcr, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  # Normalize n_reads_pcr to ensure they sum to 1 within each Sample_ID
  mutate(Proportion = n_reads_pcr / sum(n_reads_pcr, na.rm = TRUE))

# Check if each Sample_ID sums to 1
check_sums <- plot_data %>%
  group_by(Sample_ID) %>%
  summarize(Mean= sum(Proportion, na.rm = TRUE)) %>%
  filter(abs(Mean - 1) > 1e-6) # Filter samples where the sum deviates significantly from 1

# Print the check results
if (nrow(check_sums) > 0) {
  print("Warning: The following samples do not sum to 1 after normalization:")
  print(check_sums)
} else {
  print("All samples sum to 1 after normalization.")
}

# Create the stacked bar plot
# PC1 Labels
labels_for_PC1 <- metadata %>% 
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange(PC1)

ggplot(plot_data, aes(x = as.factor(PC1), y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  labs(
    title = "Proportional Stacked Bar Plot of Reads by Taxa",
    x = "Sample ID",
    y = "Relative Reads",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )+
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)+
  facet_wrap(~size_fraction_numeric)
# scale_fill_manual(values = taxa_colors)


