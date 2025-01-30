#Analysis Script for PCR Bias Correction Paper



# Read in the Data --------------------------------------------------------

# Packages and Functions --------------------------------------------------
librarian::shelf(tidyverse, googledrive, stringr,here,phyloseq,
                 extrafont, RColorBrewer, fido, compositions, ggpubr)

#Add functions for myself
source(("PCR_bias_correction/scripts/helpful_functions/treemap_funs_Capone.R"))
source("PCR_bias_correction/scripts/helpful_functions/phyloseq_mapping_funs.R")
source("PCR_bias_correction/scripts/helpful_functions/general_helper_functions.R")

saving=0


# Load in the data -----------------------------------------------------

#Metadata
env_metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  mutate(size_fraction=max_size) %>% 
  select(-X, -Sizefractionmm,-max_size) %>%
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
  rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )





# Methods Comparison. ------------------------------------------------------
# Use either proportions or CLR


# 18s ---------------------------------------------------------------------
#Using family

#Taxa file from pre-processed fido families for 18S
zhan_taxa=read.csv(here("PCR_bias_correction/data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))  %>% 
  select(-X) %>% 
  distinct() %>% 
  column_to_rownames("Family")

#CLR
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/clr/predicted_og_18s_04_29_2024_s1_phy_all_and_subpools_clr.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/clr/predicted_og_18s_04_29_2024_s2_phy_all_and_subpools_clr.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/clr/predicted_og_18s_04_29_2024_s3_phy_all_and_subpools_clr.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))

#Merge
final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#Make final dataframe
phy_taxa_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr=phy_taxa_pcr %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
  filter(Order=="Calanoida") 


#PCR-RA df ready to join with RRA
pcr_join_clr=taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))


#Proportions
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_04_17_2024_s1_phy_all_and_subpools.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_04_17_2024_s2_phy_all_and_subpools.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_04_17_2024_s3_phy_all_and_subpools.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))

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


#PCR-RA df ready to join with RRA

#All taxa
pcr_join_prop=phy_taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle,taxa) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))%>%
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )


#Vis 1. View proportion taxa in each sample
plot_data <- pcr_join_prop %>%
  group_by(Sample_ID, taxa) %>%
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
    y = "Proportion (Sum to 1)",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )



# =========== Raw Reads Relative Abundance Data using taxa that went into fido model

#Predicted proportions
fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv")) %>% 
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
  column_to_rownames("Sample_ID_dot")%>%
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )

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
  # #Add size fraction that will match with Zooscan
  # rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  # mutate(
  #   size_fraction = case_when(
  #     size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
  #     size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
  #     TRUE ~ NA_character_ # Handle any unexpected values
  #   )
  # )


#Raw in CLR
phy_18s_clr=phyloseq(OTU, TAX, meta) %>% 
  tax_glom(taxrank="Order", NArm=TRUE) %>%
  transform_sample_counts(., clr_convert)%>%
  phyloseq_transform_to_long(.)


#Filter to calanoid copepods
taxa_sel="Calanoida"


#1. Join with PCR CLR
phy_18s_clr %>% 
  filter(Order==taxa_sel) %>%
  left_join(pcr_join_clr, by=c("Sample_ID","PC1"), keep = FALSE) %>%
  select(-size_fraction.x) %>%
  mutate(n_reads_raw=n_reads) %>% 
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s_clr

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
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  rename(size_fraction=size_fraction.y)->pcr_and_raw_18s_counts

#Add PCR-bias corrected counts
counts_raw_all=phy_18s_raw_counts %>% 
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw_sum=sum(n_reads)) %>% 
  select(Sample_ID,size_fraction,n_reads_raw_sum)


pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
  left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



# Zooscan -----------------------------------------------------------------


#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read in processed Zooscan data and look at biomass proortion
zoo_metric="biomass_prop"


#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Calanoida"


#Zooscan biomass proportion
zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  filter(object_annotation_category==taxa_sel)  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) 
# mutate(size_fraction = case_when(
#   size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
# #   TRUE ~ NA_real_)) %>%
# group_by(Sample_ID) %>% 
# left_join(.,env_metadata, by=c("PC1","size_fraction")) 

# Zooscan biomass proportion in CLR
zooscan_taxa_clr= read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>% 
  group_by(sample_id,size_fraction) %>% 
  ungroup() %>% 
  filter(object_annotation_category=="Calanoida")  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) %>%
  # mutate(size_fraction = case_when(
  #   size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
  #   TRUE ~ NA_real_)) %>%
  group_by(Sample_ID) %>% 
  left_join(.,env_metadata, by=c("PC1","size_fraction"))  

# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category=="Calanoida")  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) %>%
  # mutate(size_fraction = case_when(
  #   size_fraction %in% names(size_mapping) ~ size_mapping[size_fraction],
  #   TRUE ~ NA_real_)) %>%
  group_by(Sample_ID) %>% 
  left_join(.,env_metadata, by=c("PC1","size_fraction")) 


#Zooscan relative abundance CLR
zooscan_relative_clr=zooscan_relative %>% 
  group_by(sample_id,size_fraction) %>% 
  mutate(abundance_clr=clr_convert(count))%>%
  ungroup() %>% 
  group_by(sample_id) %>% 
  left_join(.,env_metadata, by=c("PC1","size_fraction"))




# Q1. Raw Reads vs. PCR Bias Corrected ------------------------------------

plot_data=pcr_and_raw_18s_all %>%
  # Create a new variable Sample_ID_short (everything before the first '_')
  mutate(Sample_ID_short = sub("_.*", "", Sample_ID)) %>%
  # Ensure Sample_ID is ordered by increasing PC1
  mutate(Sample_ID = reorder(Sample_ID, PC1)) %>%
  group_by(Sample_ID_short, taxa, size_fraction_numeric) %>%
  # Summarize raw and PCR reads
  summarize(
    n_reads_raw = sum(n_reads_raw, na.rm = TRUE),
    n_reads_pcr = sum(n_reads_pcr, na.rm = TRUE),
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

# Create the stacked bar plot with custom labels
fig1=ggplot(plot_data, aes(x = Sample_ID_short, y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) + # Apply custom labels
  scale_fill_manual(values = taxa_colors) + # Apply custom taxa colors
  labs(
    title = "Proportional Stacked Bar Plot of Reads by Taxa",
    x = "Sample ID (Shortened)",
    y = "Proportion Reads",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

#Save
# Generate the plot
plot <- ggplot(plot_data, aes(x = Sample_ID_short, y = Proportion, fill = taxa)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  facet_grid(Metric ~ size_fraction_numeric, scales = "free_y", labeller = custom_labeller) + # Apply custom labels
  scale_fill_manual(values = taxa_colors) + # Apply custom taxa colors
  labs(
    title = "Proportional Stacked Bar Plot of Reads by Taxa",
    x = "Sample ID (Shortened)",
    y = "Proportion (Sum to 1)",
    fill = "Taxa"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

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


# Q3. Comparison with Zooscan -------------------------------------------------


