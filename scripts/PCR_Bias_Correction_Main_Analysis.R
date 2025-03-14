#Analysis Script for PCR Bias Correction Paper



# Read in the Data --------------------------------------------------------

# Packages and Functions --------------------------------------------------
librarian::shelf(tidyverse, googledrive, stringr,here,phyloseq,
                 extrafont, RColorBrewer, fido, compositions, ggpubr, patchwork)

#Add functions for myself
source(("PCR_bias_correction/scripts/helpful_functions/treemap_funs_Capone.R"))
source("PCR_bias_correction/scripts/helpful_functions/phyloseq_mapping_funs.R")
source("PCR_bias_correction/scripts/helpful_functions/general_helper_functions.R")

saving=0


# Load in the data -----------------------------------------------------

#Metadata
metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

#Add depth
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
  # mutate(
  #   size_fraction = case_when(
  #     size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
  #     size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
  #     TRUE ~ NA_character_ # Handle any unexpected values
  #   )
  # )

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
load_most_recent_file <- function(suffix) {
  # Directory path
  target_dir <- here("PCR_bias_correction/data/predicted_og")
  
  # List all files matching the pattern for the given suffix
  files <- list.files(target_dir, 
                      pattern = paste0("predicted_og_18s_\\d{2}_\\d{2}_\\d{4}_", suffix, "_phy_all_and_subpools\\.csv"), 
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

# Load files dynamically for s1, s2, and s3
fido_s1 <- load_most_recent_file("s1")
fido_s2 <- load_most_recent_file("s2")
fido_s3 <- load_most_recent_file("s3")


#Proportions (Using Data with all taxa)
# fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03__2025_s1_phy_all_and_subpools.csv")) %>%
#   select(-X) %>% 
#   mutate(coord=str_remove(coord, "^clr_"))
# fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03_06_2025_s2_phy_all_and_subpools.csv")) %>%
#   select(-X)%>% 
#   mutate(coord=str_remove(coord, "^clr_"))
# fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03_06_2025_s3_phy_all_and_subpools.csv")) %>%
#   select(-X)%>% 
#   mutate(coord=str_remove(coord, "^clr_"))

#Merge
final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 


#Make final dataframe
phy_taxa_pcr= final_data_all_sizes %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr=phy_taxa_pcr %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
  filter(Order=="Calanoida") 


#CLR
fido_s1_clr=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03_11_2025_s1_phy_all_and_subpools_clr.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2_clr=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03_11_2025_s2_phy_all_and_subpools_clr.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3_clr=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_03_11_2025_s3_phy_all_and_subpools_clr.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))

#Merge
final_data_all_sizes_clr=rbind(fido_s1_clr,fido_s2_clr,fido_s3_clr) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 



#PCR-RA df ready to join with RRA

#All taxa
pcr_join_prop=phy_taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction_numeric,PC1,cycle,taxa) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))
  #Add size fraction that will match with Zooscan
  # rename(size_fraction_numeric = size_fraction)  # Rename the existing column
  # mutate(
  #   size_fraction = case_when(
  #     size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
  #     size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
  #     TRUE ~ NA_character_ # Handle any unexpected values
  #   )
  # )


#Colors
taxa_colors <- c(
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
  "unidentified Siphonophorae" = "#CDA4DE"  # Pastel lavender
)


# Filter the data to exclude cases where zero is within the error bars

#Proportions
filtered_data <- final_data_all_sizes %>%
  filter(p2.5 > 0 | p97.5 < 0) %>%  # Keep only rows where both values are either >0 or <0
  mutate(Sample_ID = as.factor(Sample_ID),  # Ensure Sample_ID is categorical
         coord = as.factor(coord), 
         size = as.factor(size)) %>%   # Ensure size is categorical for shape mapping
  left_join(env_metadata %>% select(Sample_ID,Sample_ID_short,PC1), by="Sample_ID")

labels_for_map=filtered_data %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

#CLR
filtered_data_clr <- final_data_all_sizes_clr %>%
  filter(p2.5 > 0 | p97.5 < 0) %>%  # Keep only rows where both values are either >0 or <0
  mutate(Sample_ID = as.factor(Sample_ID),  # Ensure Sample_ID is categorical
         coord = as.factor(coord), 
         size = as.factor(size)) %>%   # Ensure size is categorical for shape mapping
  left_join(env_metadata %>% select(Sample_ID,Sample_ID_short,PC1), by="Sample_ID")

labels_for_map=filtered_data_clr %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

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
ggplot(filtered_data %>% filter(size=="0.2-0.5mm"), aes(x = as.factor(PC1), y = n_reads, color = coord, shape = size)) +
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
  scale_x_discrete(labels = labels_for_map$Sample_ID_short)+
  facet_wrap(~size, nrow=3)



#Vis 1. View proportion taxa in each sample
plot_data <- pcr_join_prop %>%
  group_by(Sample_ID, taxa, size_fraction_numeric) %>%
  # Calculate the sum of n_reads_pcr for each taxa within each Sample_ID
  summarize(n_reads_pcr = sum(n_reads_pcr, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  # Normalize n_reads_pcr to ensure they sum to 1 within each Sample_ID
  mutate(Proportion = n_reads_pcr / sum(n_reads_pcr, na.rm = TRUE))

# Check if each Sample_ID sums to 1
check_sums <- plot_data %>%
  group_by(Sample_ID) %>%
  summarize(Sum = sum(Proportion, na.rm = TRUE)) %>%
  filter(abs(Sum - 1) > 1e-6) # Filter samples where the sum deviates significantly from 1

# Print the check results
if (nrow(check_sums) > 0) {
  print("Warning: The following samples do not sum to 1 after normalization:")
  print(check_sums)
} else {
  print("All samples sum to 1 after normalization.")
}

# Create the stacked bar plot
ggplot(plot_data, aes(x = Sample_ID, y = Proportion, fill = taxa)) +
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
  scale_fill_manual(values = taxa_colors)+
  facet_wrap(~size_fraction_numeric,nrow=3)



# =========== Raw Reads Relative Abundance Data using taxa that went into fido model

fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


  #11/2024 sum in Salpidae
  # mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  # group_by(Family) %>%
  # summarise(across(everything(), sum, na.rm = TRUE)) %>%
  # ungroup() 



fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) 


  #11/2024 sum in Salpidae
  # mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  # group_by(Family) %>%
  # summarise(across(everything(), sum, na.rm = TRUE)) %>%
  # ungroup() 


fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


  #11/2024 sum in Salpidae
  # mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  # group_by(Family) %>%
  # summarise(across(everything(), sum, na.rm = TRUE)) %>%
  # ungroup() 


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


#Raw in CLR
# phy_18s_clr=phyloseq(OTU, TAX, meta) %>% 
#   tax_glom(taxrank="Order", NArm=TRUE) %>%
#   transform_sample_counts(., clr_convert)%>%
#   phyloseq_transform_to_long(.)


#Filter to calanoid copepods
taxa_sel="Calanoida"


# #1. Join with PCR CLR
# phy_18s_clr %>% 
#   filter(Order==taxa_sel) %>%
#   left_join(pcr_join_clr, by=c("Sample_ID","PC1"), keep = FALSE) %>%
#   select(-size_fraction.x) %>%
#   mutate(n_reads_raw=n_reads) %>% 
#   mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s_clr

#2. Join with PCR Proportions
#All taxa
pcr_and_raw_18s_all=
  phy_18s %>% 
  mutate(taxa=Family) %>% 
  group_by(Sample_ID,taxa) %>%
  summarise(n_reads_raw=sum(n_reads)) %>% 
  left_join(pcr_join_prop, by=c("Sample_ID","taxa"))

#3. Join in counts
phy_18s_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction_numeric) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
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



# Zooscan -----------------------------------------------------------------


#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read in processed Zooscan data and look at biomass proortion
zoo_metric="biomass_prop"



#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Calanoida"


#Zooscan biomass proportion
# All taxa
zooscan_all=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) 

#Calanoida
zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  filter(object_annotation_category==taxa_sel)  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) 

# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category==taxa_sel)  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) %>%
  group_by(Sample_ID) %>% 
  left_join(.,metadata, by=c("PC1","size_fraction_numeric")) 






# Q1. Raw Reads vs. PCR Bias Corrected ------------------------------------

plot_data=pcr_and_raw_18s_all %>%
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
    "0.5" = "0.2-0.5mm",
    "1" = "0.5-1mm",
    "2" = "1-2mm"
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


fig1 <- ggplot(plot_data, aes(x = as.factor(PC1), y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) + # Apply custom labels
  scale_fill_manual(values = taxa_colors) + # Apply custom taxa colors
  scale_x_discrete(labels = labels_for_PC1$Sample_ID_short)+
  labs(
    title = "Stacked Bar Plot of Relative Read Abundances by Taxa",
    x = "Sample ID",
    y = "Proportion Reads",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # X-axis text size
    axis.text.y = element_text(size = 10), # Y-axis text size
    axis.title.x = element_text(size = 12), # X-axis label size
    axis.title.y = element_text(size = 12), # Y-axis label size
    strip.text = element_text(face = "bold", size = 10), # Facet text size
    legend.text = element_text(size = 10), # Legend text size
    legend.title = element_text(size = 10), # Legend title size
    title = element_text(size = 12), # Plot title size
    legend.position = "bottom" # Legend position
  )
fig1

# Define output directory
output_dir <- "PCR_bias_correction/figures/v0"

# Save the plot as PNG
ggsave(
  filename = file.path(output_dir, "figure_2_v0.png"),
  plot = fig1,
  dpi = 300,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)

# Save the plot as PDF
ggsave(
  filename = file.path(output_dir, "figure_2_v0.pdf"),
  plot = fig1,
  dpi = 300,
  width = 7, # Adjust width (inches) based on journal requirements
  height = 5, # Adjust height (inches) based on journal requirements
  units = "in"
)

# Q2. Taxa Effects of PCR Bias --------------------------------------------


#Figure. Amplification efficiencies by taxa


#Figure. Amplification efficiencies vs. Proportion reads


# Q3. Comparison with Zooscan for calanoids-------------------------------------------------
pcr_and_raw_18s_calanoids=pcr_and_raw_18s_all %>% 
  mutate(Family=taxa) %>% 
  left_join(.,zhan_taxa %>% rownames_to_column("Family") %>% 
              select(Family, Order), by="Family") %>% 
  filter(Order=="Calanoida") %>% 
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


#Format Long and add difference metrics
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s, 
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





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions


pcr_raw_zoo_18s_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance_metric, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y") +
  theme_minimal() +
  labs(title = "Methods Differences 18S",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Biomass ",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"),
                    labels = c("Zooscan Biomass", "PCR-corrected", "Raw Reads")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+# Define linetypes for each Method
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F")) -> grouped_bar_all_18s  # Define colors for each Method

grouped_bar_all_18s


#Scatter Plot
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle))%>%
  ggplot(.,aes(x=((biomass_prop_taxa)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  # facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr
zoo_vs_pcr
