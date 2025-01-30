library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(patchwork)
detach("package:fido", unload = TRUE)

# Make dfs ----------------------------------------------------------------

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
metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)


dryweights=read.csv(here("Zoop_Patterns/data/zoop_other/biomass_processed.csv"))

  #Add size fraction that will match with Zooscan
  # rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  # mutate(
  #   size_fraction = case_when(
  #     size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
  #     size_fraction_numeric >= 1 & size_fraction_numeric <= 2 ~ "1-2", # Assign 1-5
  #     size_fraction_numeric > 2 ~ ">2", # Assign 1-5
  #     TRUE ~ NA_character_ # Handle any unexpected values
  #   )
  # )



env_metadata = metadata %>% 
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))




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
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_11_22_2024_s1_phy_all_and_subpools.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_11_22_2024_s2_phy_all_and_subpools.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_11_22_2024_s3_phy_all_and_subpools.csv")) %>%
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
pcr_join_prop=taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,max_size,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))%>%
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = max_size) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-2", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
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
#Filter to onshore offshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID=Sample_ID_dot) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1) %>% unique()


env_metadata_phy=env_metadata%>% 
  left_join(.,clusters, by="Sample_ID") %>%
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
phy_18s_clr=phyloseq(OTU, TAX, meta) %>% 
  tax_glom(taxrank="Order", NArm=TRUE) %>%
  transform_sample_counts(., clr_convert)%>%
  phyloseq_transform_to_long(.)
  # #Add size fraction that will match with Zooscan
  # rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  # mutate(
  #   size_fraction = case_when(
  #     size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
  #     size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
  #     TRUE ~ NA_character_ # Handle any unexpected values
  #   )
  # )

#Filter to calanoid copepods
taxa_sel="Calanoida"


#Join with PCR CLR
phy_18s_clr %>% 
  filter(Order==taxa_sel) %>%
  left_join(pcr_join_clr, by=c("Sample_ID","PC1"), keep = FALSE) %>%
  select(-size_fraction.x) %>%
  mutate(n_reads_raw=n_reads) %>% 
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s_clr

#Join with PCR Proportions
pcr_and_raw_18s=phy_18s %>% 
  # filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction, Family,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Family, and PC1
  group_by(size_fraction, Family, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}")) %>% 
  rename(Sample_ID=Sample_ID_short)

calanoida_dna=phy_18s %>% 
  filter(Order == "Calanoida", !Family %in% c("Eucalanidae", "Rhincalanidae")) %>% 
  group_by(Sample_ID,size_fraction, Order,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Family, and PC1
  group_by(size_fraction, Order, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

eucalanidae_dna=phy_18s %>% 
  filter(Family %in% c("Eucalanidae","Rhincalanidae")) %>%
  group_by(Sample_ID,size_fraction, Family,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Family, and PC1
  group_by(size_fraction, Family, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

oithona_dna=phy_18s %>% 
  filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Family,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Family, and PC1
  group_by(size_fraction, Family, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

euphausiid_dna=phy_18s %>% 
  filter(Order=="Euphausiacea") %>%
  # filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction,Order,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Family, and PC1
  group_by(size_fraction, Order, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)
  

#Join in counts
phy_18s_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction, Family) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  rename(size_fraction=size_fraction.y)->pcr_and_raw_18s_counts

#Add PCR-bias corrected counts
counts_raw_all=phy_18s_raw_counts %>% 
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw_sum=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>% 
  select(Sample_ID,size_fraction,n_reads_raw_sum)


# pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
#   left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction")) %>% 
#   mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



# Zooscan -----------------------------------------------------------------

#Zooscan biomass proportion: Filter to offshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID_short=Sample_ID) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1,Sample_ID) %>% unique()


zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  select(-X)  %>%
  mutate(
    Sample_ID = sample_id, 
    size_fraction = case_when(
      size_fraction=="1-5" ~ "1-2",  # Map values 1-5 to "1-2"
      TRUE ~ as.character(size_fraction)  # Retain other values as is
    )
  ) %>%
  left_join(., clusters, by = "Sample_ID_short")


## Create a taxa mapping between Zooscan and DNA

#First modify the DNA df
# Define the mapping
family_to_taxa <- c(
  "Calanidae" = "Calanoida",
  "Clausocalanidae" = "Calanoida",
  "Eucalanidae" = "Eucalanidae",
  "Euphausiidae" = "Euphausiacea",
  "Metridinidae" = "Calanoida",
  "Oithonidae" = "Oithonidae",
  "Paracalanidae" = "Calanoida",
  "Rhincalanidae" = "Calanoida",
  "Salpidae" = "Salpida",
  "other" = "other",
  "unidentified Calanoida" = "Calanoida"
)

# Add the taxa_zooscan column using mutate and recode
pcr_and_raw_18s <- pcr_and_raw_18s %>%
  mutate(taxa_join = recode(Family, !!!family_to_taxa))

#Next Zooscan
# Reverse mapping definition
reverse_mapping <- c(
  "Calanoida" = "Calanoida",
  "Euphausiacea" = "Euphausiacea",
  "Eucalanidae" = "Eucalanidae",
  "Oithonidae" = "Oithonidae",
  "Salpida" = "Salpida",
  "other" = "other"
)

# Perform reverse mapping
zooscan_taxa <- zooscan_taxa %>%
  mutate(taxa_join = case_when(
    object_annotation_category %in% names(reverse_mapping) ~ reverse_mapping[object_annotation_category],
    TRUE ~ "other" # Map all other categories to "other"
  ))




zooscan_calanoida=zooscan_taxa %>% 
  filter(object_annotation_category=="Calanoida")


zooscan_eucalanus=zooscan_taxa %>% 
  filter(object_annotation_category=="Eucalanidae")


zooscan_euphausiid=zooscan_taxa %>%  
  filter(object_annotation_category %in% c("Euphausiacea", "Eumalacostraca"))

zooscan_oithona=zooscan_taxa %>% 
  filter(object_annotation_category=="Oithonidae")

# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category==taxa_sel)  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id)


#========== COMPARE: Make combined dataframe for comparing all 3 methods

#Propotions

#All
pcr_raw_zoo_18s_all=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction","taxa_join")) %>% 
  unique(.) 
  filter(taxa_join != "other")

#Save
# write.csv(pcr_raw_zoo_18s_all,here("PCR_bias_correction/data/methods_compare/pcr_raw_zoo_18s_all_taxa.csv"))

#Calanoida
pcr_raw_zoo_18s_calanoida=zooscan_calanoida %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(calanoida_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Euphausiids
pcr_raw_zoo_18s_euphausiid=zooscan_euphausiid %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(eucalanidae_dna, by=c("PC1","size_fraction")) %>% 
  unique(.) %>% 
  filter(biomass_prop_taxa>1e-3)


#Eucalanidae
pcr_raw_zoo_18s_eucalanidae=zooscan_eucalanus %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(eucalanidae_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Oithonidae
pcr_raw_zoo_18s_oithonidae=zooscan_oithona %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(oithona_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Counts
pcr_raw_zoo_18s_counts=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">5") %>%
  left_join(pcr_and_raw_18s_counts, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Abundance
pcr_raw_zoo_18s_relab=zooscan_relative %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction"))




# Plots --------------------------------------------------------------------
#Calanoida first

  my_palette=custom_pallete()
# 
# # Function to calculate R² and p-value for each group
# assess_correlation_by_group <- function(data, x_var, y_var, group_var) {
#   data %>%
#     group_by(.data[[group_var]]) %>%
#     summarise(
#       r_squared = {
#         model <- lm(as.formula(paste(y_var, "~", x_var)), data = cur_data())
#         summary(model)$r.squared
#       },
#       p_value = {
#         model <- lm(as.formula(paste(y_var, "~", x_var)), data = cur_data())
#         coef(summary(model))["biomass_prop_taxa", "Pr(>|t|)"]
#       },
#       .groups = "drop"  # Avoids grouped dataframes in the result
#     )
# }
# 
# # Add R² and p-value annotations for each facet
# add_r2_p_annotations <- function(data, x_var, y_var, group_var) {
#   correlation_results <- assess_correlation_by_group(data, x_var, y_var, group_var)
#   
#   # Add annotations to the data for plotting
#   data <- data %>%
#     left_join(correlation_results, by = c(offshore_onshore = group_var)) %>%
#     mutate(
#       annotation = paste0("R²: ", round(r_squared, 2), "\nP: ", signif(p_value, 3))
#     )
#   return(data)
# }
# 
# # Prepare data with annotations for PCR
# data_pcr <- add_r2_p_annotations(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_pcr", "offshore_onshore")
# 
# # Prepare data with annotations for Raw Reads
# data_raw <- add_r2_p_annotations(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_raw", "offshore_onshore")
# 
# 
# data_pcr_size <- add_r2_p_annotations(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_pcr", "size_fraction")
# 
# # Prepare data with annotations for Raw Reads
# data_raw_size <- add_r2_p_annotations(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_raw", "size_fraction")



# Plotting function
create_facet_plot <- function(data, x_var, y_var, title, x_transform = "none", y_transform = "none", group_var = "offshore_onshore") {
  
  # Apply transformations to x_var and y_var
  data <- data %>%
    mutate(
      transformed_x_var = case_when(
        x_transform == "log10" ~ log10(.data[[x_var]] + 1), # Add 1 to avoid log10(0)
        x_transform == "arcsinsqrt" ~ asin(sqrt(.data[[x_var]])),
        TRUE ~ .data[[x_var]] # No transformation
      ),
      transformed_y_var = case_when(
        y_transform == "log10" ~ log10(.data[[y_var]] + 1), # Add 1 to avoid log10(0)
        y_transform == "arcsinsqrt" ~ asin(sqrt(.data[[y_var]])),
        TRUE ~ .data[[y_var]] # No transformation
      )
    )
  
  # Nested function to calculate R² and p-value for each group or the entire dataset
  assess_correlation_by_group <- function(data, x_var, y_var, group_var) {
    if (is.null(group_var) || group_var == "none") {
      # No grouping: calculate correlation on the entire dataset
      tibble(
        r_squared = {
          model <- tryCatch(
            lm(as.formula(paste(y_var, "~", x_var)), data = data),
            error = function(e) NULL
          )
          if (!is.null(model)) summary(model)$r.squared else NA
        },
        p_value = {
          model <- tryCatch(
            lm(as.formula(paste(y_var, "~", x_var)), data = data),
            error = function(e) NULL
          )
          if (!is.null(model) && "Pr(>|t|)" %in% colnames(coef(summary(model)))) {
            coef(summary(model))[[2, "Pr(>|t|)"]]
          } else {
            NA
          }
        }
      )
    } else {
      # Grouped calculation
      data %>%
        group_by(.data[[group_var]]) %>%
        summarise(
          r_squared = {
            model <- tryCatch(
              lm(as.formula(paste(y_var, "~", x_var)), data = cur_data()),
              error = function(e) NULL
            )
            if (!is.null(model)) summary(model)$r.squared else NA
          },
          p_value = {
            model <- tryCatch(
              lm(as.formula(paste(y_var, "~", x_var)), data = cur_data()),
              error = function(e) NULL
            )
            if (!is.null(model) && "Pr(>|t|)" %in% colnames(coef(summary(model)))) {
              coef(summary(model))[[2, "Pr(>|t|)"]]
            } else {
              NA
            }
          },
          .groups = "drop"  # Avoids grouped dataframes in the result
        )
    }
  }
  
  # Assess correlation and add annotations to the data
  if (is.null(group_var) || group_var == "none") {
    correlation_results <- assess_correlation_by_group(data, "transformed_x_var", "transformed_y_var", NULL)
    data <- data %>%
      mutate(
        annotation = paste0("R²: ", ifelse(is.na(correlation_results$r_squared), "NA", round(correlation_results$r_squared, 2)), 
                            "\nP: ", ifelse(is.na(correlation_results$p_value), "NA", signif(correlation_results$p_value, 3)))
      )
  } else {
    correlation_results <- assess_correlation_by_group(data, "transformed_x_var", "transformed_y_var", group_var)
    data <- data %>%
      left_join(correlation_results, by = c(group_var)) %>%
      mutate(
        annotation = paste0("R²: ", ifelse(is.na(r_squared), "NA", round(r_squared, 2)), 
                            "\nP: ", ifelse(is.na(p_value), "NA", signif(p_value, 3)))
      )
  }
  
  # Plot using the transformed variables and annotations
  # Create the ggplot object
  plot <- ggplot(data, aes(x = transformed_x_var, y = transformed_y_var)) +
    geom_point(size = 4, aes(shape = cycle.x, fill = cycle.x)) +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
    scale_fill_manual(values = my_palette) +
    scale_color_manual(name = "Cycle", values = my_palette) +  # Change legend label to "Cycle"
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
    labs(
      title = title,
      x = ifelse(x_transform == "log10", "log10 Transformed",
                 ifelse(x_transform == "arcsinsqrt", "Arcsin-Sqrt Transformed", "Proportion Zooscan Biomass")),
      y = ifelse(y_transform == "log10", "log10 Transformed",
                 ifelse(y_transform == "arcsinsqrt", "Arcsin-Sqrt Transformed", "Proportion Reads"))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center main title
      strip.text = element_text(size = 10, face = "bold", hjust = 0.5),  # Center and bold facet titles
      axis.text = element_text(size = 10),  # Tick label size set to 10
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.spacing = unit(1, "lines")
    )
  
  # Add faceting and annotations dynamically
  if (!is.null(group_var) && group_var != "none") {
    plot <- plot +
      facet_wrap(as.formula(paste0("~", group_var)), labeller = labeller(.default = toupper)) +
      geom_text(data = distinct(data, !!sym(group_var), annotation), 
                aes(x = Inf, y = Inf, label = annotation), 
                hjust = 2.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
  } else {
    plot <- plot +
      geom_text(aes(x = Inf, y = Inf, label = annotation), 
                hjust = 2.1, vjust = 1.1, size = 4, inherit.aes = FALSE)
  }
  
  # Return the plot
  return(plot)
}



# Create individual plots


#1. Size Fraction

#No transformation
plot_pcr <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var="size_fraction"
)
plot_raw <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var="size_fraction"
)
# Combine plots into a grid
combined_plot <- plot_pcr / plot_raw +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot



#Arcsin
plot_pcr_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="size_fraction"
)

plot_raw_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="size_fraction"
)

# Combine plots into a grid
combined_plot_arcsin <- plot_pcr_arcsin / plot_raw_arcsin +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_arcsin




#Log
plot_pcr_log <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                     x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="size_fraction"
)

plot_raw_log <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                     x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="size_fraction"
)

# Combine plots into a grid
combined_plot_log <- plot_pcr_log / plot_raw_log +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_log


# 2. Offshore Onshore
# Create individual plots

#No transformation
plot_pcr <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var="offshore_onshore"
)
plot_raw <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var="offshore_onshore"
)
# Combine plots into a grid
combined_plot <- plot_pcr / plot_raw +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot



#Arcsin
plot_pcr_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="offshore_onshore"
)

plot_raw_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="offshore_onshore"
)

# Combine plots into a grid
combined_plot_arcsin <- plot_pcr_arcsin / plot_raw_arcsin +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_arcsin




#Log
plot_pcr_log <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                  x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                  y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                  group_var="offshore_onshore"
)

plot_raw_log <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                  x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                  y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                  group_var="offshore_onshore"
)

# Combine plots into a grid
combined_plot_log <- plot_pcr_log / plot_raw_log +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_log


# 3. No Grouping
#No transformation
plot_pcr <- create_facet_plot(pcr_raw_zoo_18s_euphausiid, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var=NULL
)
plot_raw <- create_facet_plot(pcr_raw_zoo_18s_euphausiid, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                              x_transform = "none",    # Options: "log10", "arcsinsqrt", "none"
                              y_transform = "none", # Options: "log10", "arcsinsqrt", "none"
                              group_var=NULL
)
# Combine plots into a grid
combined_plot <- plot_pcr / plot_raw +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot


#Arcsin
plot_pcr_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="none"
)

plot_raw_arcsin <- create_facet_plot(pcr_raw_zoo_18s_eucalanidae, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                     x_transform = "arcsinsqrt",    # Options: "log10", "arcsinsqrt", "none"
                                     y_transform = "arcsinsqrt", # Options: "log10", "arcsinsqrt", "none"
                                     group_var="none"
)

# Combine plots into a grid
combined_plot_arcsin <- plot_pcr_arcsin / plot_raw_arcsin +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_arcsin




#Log
plot_pcr_log <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected",
                                  x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                  y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                  group_var="none"
)

plot_raw_log <- create_facet_plot(pcr_raw_zoo_18s_calanoida, "biomass_prop_taxa", "n_reads_raw", "Calanoida - PCR Bias Corrected",
                                  x_transform = "log10",    # Options: "log10", "arcsinsqrt", "none"
                                  y_transform = "log10", # Options: "log10", "arcsinsqrt", "none"
                                  group_var="none"
)

# Combine plots into a grid
combined_plot_log <- plot_pcr_log / plot_raw_log +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot_log



# Save the combined plot
output_path <- "Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/methods_compare_calanoida_offshore_onshore.png"

# Use ggsave-like saving for grid.arrange
# Convert mm to inches (169mm = 6.65 inches)
ggsave(
  filename = output_path,
  plot = combined_plot,
  width = 6.65,  # Width in inches
  height = 6.65,  # Height scaled to retain aspect ratio (2:1)
  units = "in",
  dpi = 300
)






# Stacked Bar -------------------------------------------------------------

# Prepare the data
pcr_raw_zoo_18s_long <- pcr_raw_zoo_18s_calanoida %>%
  pivot_longer(
    cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr),
    names_to = "Method",
    values_to = "relative_abundance"
  ) %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(
    diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
    diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance)
  ) %>%
  mutate(is_closer = ifelse(
    Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
    ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA)
  )) %>%
  mutate(
    closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05 * relative_abundance, "**", NA),
    worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05 * relative_abundance, "xx", NA)
  )

labels_for_map <- pcr_raw_zoo_18s_calanoida %>%
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange((PC1))

# Create the plot
grouped_bar_all_18s <- pcr_raw_zoo_18s_long %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y",
             labeller = labeller(size_fraction = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"))) +
  theme_minimal() +
  labs(
    title = "Methods Differences 18S",
    x = expression(paste("Offshore ", PC1, " Onshore")),
    y = "Proportion Reads or Biomass",
    fill = "Method"
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(
    values = c("#70BF41", "#4F86F7", "#F78D4F"),
    labels = c("Zooscan Biomass", "PCR-corrected", "Raw Reads")
  ) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short) +
  scale_shape_manual(values = c(16, 17, 18)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"))

# Display the plot
grouped_bar_all_18s




# COI ---------------------------------------------------------------------

#Proportions
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s1_phy_all_and_subpools.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s2_phy_all_and_subpools.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s3_phy_all_and_subpools.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))

#Merge
final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 

#Taxa file from pre-processed fido families for coi
coi_taxa=read.csv(here("PCR_bias_correction/data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv")) %>%
  mutate(Genus = ifelse(Genus == "Genus", paste("unidentified",Family), Genus)) %>%
  filter(Genus %in% rownames(fido_coi_merged_raw)) %>% 
  column_to_rownames("Genus") %>% 
  select(-X) 

#Make final dataframe
phy_taxa_pcr= final_data_all_sizes %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  left_join(.,env_metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr=phy_taxa_pcr %>% mutate(Genus=taxa) %>%
  left_join(.,coi_taxa %>% rownames_to_column("Genus"), by="Genus") %>%
  filter(Order=="Calanoida")


#PCR-RA df ready to join with RRA
pcr_join_prop=phy_taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,max_size,PC1,cycle) %>%
  summarise(n_reads_pcr=sum(n_reads),
            p.2.5=min(p2.5),
            p.97.5=max(p97.5))%>%
  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = max_size) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-2", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )










# =========== Raw Reads Relative Abundance Data using taxa that went into fido model

#All merged
phy_pcr= final_data_all_sizes %>%
  filter(cycle_num==0) %>% 
  group_by(Sample_ID,size) %>%
  summarise(n_reads=sum(n_reads)) %>%
  left_join(.,env_metadata, by="Sample_ID")

#By species
phy_taxa_pcr= final_data_all_sizes %>%
  # filter(str_detect(coord, "Calanoida")) %>%
  filter(cycle_num==0) %>% 
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) %>%
  mutate(taxa=coord)



# Raw Reads ---------------------------------------------------------------

#Predicted proportions
fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Genus, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Genus) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)


merge(fido_s1_raw, fido_s2_raw, by = "Genus", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Genus", all = TRUE)%>%
  column_to_rownames("Genus") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_coi_merged_raw

#Taxa
coi_taxa=read.csv(here("PCR_bias_correction/data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv")) %>%
  mutate(Genus = ifelse(Genus == "Genus", paste("unidentified",Family), Genus)) %>%
  filter(Genus %in% rownames(fido_coi_merged_raw)) %>% 
  column_to_rownames("Genus") %>% 
  select(-X) 


#Make phyloseq objects

#coi
OTU = otu_table(as.matrix(fido_coi_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(env_metadata_phy)
Phy_raw_coi=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.)

phy_coi=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Genus=asv_code) %>%
  select(-asv_code)


#Join with PCR Proportions
pcr_and_raw_coi=phy_coi %>% 
  group_by(Sample_ID,size_fraction, Genus,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Genus, and PC1
  group_by(size_fraction, Genus, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}")) %>% 
  rename(Sample_ID=Sample_ID_short)

calanoida_dna=phy_coi %>% 
  filter(Order == "Calanoida", !Genus %in% c("Eucalanus")) %>% 
  filter(Genus=="Calanus") %>% 
  group_by(Sample_ID,size_fraction, Order,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Genus, and PC1
  group_by(size_fraction, Order, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

eucalanidae_dna=phy_coi %>% 
  filter(Genus %in% c("Eucalanus")) %>%
  group_by(Sample_ID,size_fraction, Genus,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Genus, and PC1
  group_by(size_fraction, Genus, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

oithona_dna=phy_coi %>% 
  filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Genus,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Genus, and PC1
  group_by(size_fraction, Genus, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)

euphausiid_dna=phy_coi %>% 
  filter(Order=="Euphausiacea") %>%
  # filter(Genus=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction,Order,Sample_ID_short) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)%>%
  ungroup() %>% 
  select(-Sample_ID) %>% 
  # Group and summarize across all columns by size_fraction, Genus, and PC1
  group_by(size_fraction, Order, Sample_ID_short) %>%
  summarise(across(everything(), ~ if(is.numeric(.)) mean(., na.rm = TRUE) else unique(.), .names = "{.col}"))%>% 
  rename(Sample_ID=Sample_ID_short)


#Join in counts
phy_coi_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction, Genus) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  rename(size_fraction=size_fraction.y)->pcr_and_raw_coi_counts

#Add PCR-bias corrected counts
counts_raw_all=phy_coi_raw_counts %>% 
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw_sum=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>% 
  select(Sample_ID,size_fraction,n_reads_raw_sum)


# pcr_and_raw_coi_counts=pcr_and_raw_coi_counts %>% 
#   left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction")) %>% 
#   mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



# Zooscan -----------------------------------------------------------------

#Zooscan biomass proportion: Filter to offshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID_short=Sample_ID) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1,Sample_ID) %>% unique()


zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  select(-X)  %>%
  mutate(
    Sample_ID = sample_id, 
    size_fraction = case_when(
      size_fraction=="1-5" ~ "1-2",  # Map values 1-5 to "1-2"
      TRUE ~ as.character(size_fraction)  # Retain other values as is
    )
  ) %>%
  left_join(., clusters, by = "Sample_ID_short")



zooscan_calanoida=zooscan_taxa %>% 
  filter(object_annotation_category=="Calanoida")


zooscan_eucalanus=zooscan_taxa %>% 
  filter(object_annotation_category=="Eucalanidae")


zooscan_euphausiid=zooscan_taxa %>%  
  filter(object_annotation_category %in% c("Euphausiacea", "Eumalacostraca"))

zooscan_oithona=zooscan_taxa %>% 
  filter(object_annotation_category=="Oithonidae")

# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category==taxa_sel)  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id)


#========== COMPARE: Make combined dataframe for comparing all 3 methods

#Propotions

#Calanoida
pcr_raw_zoo_coi_calanoida=zooscan_calanoida %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(calanoida_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Euphausiids
pcr_raw_zoo_coi_euphausiid=zooscan_euphausiid %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(euphausiid_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Eucalanidae
pcr_raw_zoo_coi_eucalanidae=zooscan_eucalanus %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(eucalanidae_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Oithonidae
pcr_raw_zoo_coi_oithonidae=zooscan_oithona %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(oithona_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Counts
pcr_raw_zoo_coi_counts=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">5") %>%
  left_join(pcr_and_raw_coi_counts, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Abundance
pcr_raw_zoo_coi_relab=zooscan_relative %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_coi, by=c("PC1","size_fraction"))




# Plots --------------------------------------------------------------------
#Calanoida first

# Function to calculate R² and p-value for each group
assess_correlation_by_group <- function(data, x_var, y_var, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      r_squared = {
        model <- lm(as.formula(paste(y_var, "~", x_var)), data = cur_data())
        summary(model)$r.squared
      },
      p_value = {
        model <- lm(as.formula(paste(y_var, "~", x_var)), data = cur_data())
        coef(summary(model))["biomass_prop_taxa", "Pr(>|t|)"]
      },
      .groups = "drop"  # Avoids grouped dataframes in the result
    )
}

# Add R² and p-value annotations for each facet
add_r2_p_annotations <- function(data, x_var, y_var, group_var) {
  correlation_results <- assess_correlation_by_group(data, x_var, y_var, group_var)
  
  # Add annotations to the data for plotting
  data <- data %>%
    left_join(correlation_results, by = c(offshore_onshore = group_var)) %>%
    mutate(
      annotation = paste0("R²: ", round(r_squared, 2), "\nP: ", signif(p_value, 3))
    )
  return(data)
}

# Prepare data with annotations for PCR
data_pcr <- add_r2_p_annotations(pcr_raw_zoo_coi_calanoida, "biomass_prop_taxa", "n_reads_pcr", "offshore_onshore")

# Prepare data with annotations for Raw Reads
data_raw <- add_r2_p_annotations(pcr_raw_zoo_coi_calanoida, "biomass_prop_taxa", "n_reads_raw", "offshore_onshore")

# Plotting function
create_facet_plot <- function(data, x_var, y_var, title, y_limit, x_limit) {
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(size = 4, aes(shape = cycle.x, fill = cycle.x)) +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24, "T1" = 23, "T2" = 25)) +
    scale_fill_manual(values = my_palette) +
    scale_color_manual(name = "Cycle", values = my_palette) +  # Change legend label to "Cycle"
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
    facet_wrap(~offshore_onshore, labeller = labeller(offshore_onshore = toupper)) +  # Capitalize subplot titles
    labs(
      title = title,
      x = "Zooscan Proportion Biomass",  # Updated x-axis label
      y = ifelse(y_var == "n_reads_pcr", "PCR Bias Corrected RRA", "Raw Reads")
    ) +
    geom_text(data = distinct(data, offshore_onshore, annotation), 
              aes(x = Inf, y = Inf, label = annotation), 
              hjust = 2.1, vjust = 1.1, size = 4, inherit.aes = FALSE) +
    # coord_cartesian(ylim = y_limit, xlim = x_limit) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Center main title
      strip.text = element_text(size = 10, face = "bold", hjust = 0.5),  # Center and bold facet titles
      axis.text = element_text(size = 10),  # Tick label size set to 10
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.spacing = unit(1, "lines"))
}

# Create individual plots
plot_pcr <- create_facet_plot(data_pcr, "biomass_prop_taxa", "n_reads_pcr", "Calanoida - PCR Bias Corrected", y_limits, x_limits)
plot_raw <- create_facet_plot(data_raw, "biomass_prop_taxa", "n_reads_raw", "Calanoida - Raw Reads", y_limits, x_limits)

# Combine plots into a grid
combined_plot <- plot_pcr / plot_raw +
  plot_annotation(title = "") +  # No title for combined plots
  plot_layout(guides = "collect")  # Collect legends, if any

combined_plot


# Save the combined plot
output_path <- "Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/methods_compare_calanoida_offshore_onshore.png"

# Use ggsave-like saving for grid.arrange
# Convert mm to inches (169mm = 6.65 inches)
ggsave(
  filename = output_path,
  plot = combined_plot,
  width = 6.65,  # Width in inches
  height = 6.65,  # Height scaled to retain aspect ratio (2:1)
  units = "in",
  dpi = 300
)






# Stacked Bar -------------------------------------------------------------

# Prepare the data
pcr_raw_zoo_coi_long <- pcr_raw_zoo_coi_calanoida %>%
  pivot_longer(
    cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr),
    names_to = "Method",
    values_to = "relative_abundance"
  ) %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(
    diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
    diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance)
  ) %>%
  mutate(is_closer = ifelse(
    Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
    ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA)
  )) %>%
  mutate(
    closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05 * relative_abundance, "**", NA),
    worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05 * relative_abundance, "xx", NA)
  )

labels_for_map <- pcr_raw_zoo_coi_calanoida %>%
  ungroup() %>%
  select(Sample_ID_short, PC1) %>%
  unique() %>%
  arrange((PC1))

# Create the plot
grouped_bar_all_coi <- pcr_raw_zoo_coi_long %>%
  ggplot(aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y",
             labeller = labeller(size_fraction = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"))) +
  theme_minimal() +
  labs(
    title = "Methods Differences coi",
    x = expression(paste("Offshore ", PC1, " Onshore")),
    y = "Proportion Reads or Biomass",
    fill = "Method"
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  ) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(
    values = c("#70BF41", "#4F86F7", "#F78D4F"),
    labels = c("Zooscan Biomass", "PCR-corrected", "Raw Reads")
  ) +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short) +
  scale_shape_manual(values = c(16, 17, 18)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F"))

# Display the plot
grouped_bar_all_coi




# Scrap -------------------------------------------------------------------

# 
# 
# 
# 
# # Define the plot function
# create_plot <- function(data, title) {
#   # Calculate limits for the axes
#   # x_max <- max(data$dryweight_C_mg_m2_taxa, na.rm = TRUE) * 1.2
#   # y_max_pcr <- max(data$n_reads_pcr, na.rm = TRUE) * 1.2
#   # y_max_raw <- max(data$n_reads_raw, na.rm = TRUE) * 1.2
#   
#   # Create the left plot (original)
#   plot_left <- data %>%
#     ggplot(aes(x = biomass_prop_taxa, y = n_reads_pcr)) +
#     geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
#     geom_smooth(method = "lm", se = FALSE) +
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#     labs(
#       title = title,
#       x = "Zooscan Dry Weight (C mg/m²)",
#       y = "PCR Bias Corrected RRA",
#       color = "Cycle"
#     ) +
#     # xlim(0, x_max) +  # Set x-axis limits
#     # ylim(0, y_max_pcr) +  # Set y-axis limits
#     theme_minimal() +
#     theme(
#       strip.text = element_text(size = 10),
#       axis.text = element_text(size = 10),
#       axis.title = element_text(size = 12),
#       panel.spacing = unit(1, "lines")
#     )
#   
#   # Create the right plot (with n_reads_raw on y-axis)
#   plot_right <- data %>%
#     ggplot(aes(x = biomass_prop_taxa, y = n_reads_raw)) +
#     geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
#     geom_smooth(method = "lm", se = FALSE) +
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#     labs(
#       title = paste(title, "- Raw Reads"),
#       x = "Zooscan Dry Weight (C mg/m²)",
#       y = "Raw Reads",
#       color = "Cycle"
#     ) +
#     # xlim(0, x_max) +  # Set x-axis limits
#     # ylim(0, y_max_raw) +  # Set y-axis limits
#     theme_minimal() +
#     theme(
#       strip.text = element_text(size = 10),
#       axis.text = element_text(size = 10),
#       axis.title = element_text(size = 12),
#       panel.spacing = unit(1, "lines")
#     )
#   
#   # Combine the two plots side by side using patchwork
#   combined_plot <- plot_left + plot_right
#   
#   return(combined_plot)
# }
# 
# 
# # 1) Zooscan biomass vs PCR-RRA
# # Create individual plots
# plot_calanoida <- create_plot(pcr_raw_zoo_18s_calanoida, "Calanoida")
# plot_calanoida
# plot_euphausiid <- create_plot(pcr_raw_zoo_18s_euphausiid, "Euphausiids")
# plot_euphausiid
# plot_eucalanidae <- create_plot(pcr_raw_zoo_18s_eucalanidae, "Eucalanidae")
# plot_eucalanidae
# plot_oithonidae <- create_plot(pcr_raw_zoo_18s_oithonidae, "Oithonidae")
# plot_oithonidae
# 
# # Combine plots into a 2x2 grid using gridExtra
# 
# 
# # Arrange the plots in a 2x2 grid
# grid.arrange(
#   plot_calanoida, plot_euphausiid,
#   plot_eucalanidae, plot_oithonidae,
#   ncol = 2, nrow = 2
# )
# 
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_proportions_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_pcr_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_proportions.pdf",
#     width = 7, height = 9)
# print(zoobio_pcr_plot)
# dev.off()
# 
# # 2) Zooscan biomass vs Raw RRA
# zoobio_raw_plot <- pcr_raw_zoo_18s %>%
#   ggplot(aes(x = dryweight_C_mg_m2_taxa, y = n_reads_raw, color = cycle.x)) +
#   geom_point(size = 3) +  # Scatter plot
#   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#   # facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
#   labs(
#     x = "Calanoida  Zooscan Dry Weight (C mg/m²)",
#     y = "Calanoida Raw RRA",
#     color = "Cycle",
#     title = "PCR Reads vs Biomass across Cycles and Size Fractions"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   )
# 
# zoobio_raw_plot
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_proportions_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_raw_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_proportions.pdf",
#     width = 7, height = 9)
# print(zoobio_raw_plot)
# dev.off()
# 
# 
# 
# # 3) PCR Counts
# zoobio_pcr_plot <- pcr_raw_zoo_18s_counts %>%
#   ggplot(aes(x = dryweight_C_mg_m2_taxa, y = n_reads_pcr_counts, color = cycle.x)) +
#   geom_point(size = 3) +  # Scatter plot
#   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#   facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
#   labs(
#     x = "Calanoida  Zooscan Dry Weight (C mg/m²)",
#     y = "Calanoida PCR Bias Corrected Reads",
#     color = "Cycle"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   )
# 
# zoobio_pcr_plot
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_reads_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_pcr_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_reads.pdf",
#     width = 7, height = 9)
# print(zoobio_pcr_plot)
# dev.off()
# 
# 
# 
# 
# 
# #4) Raw Counts
# zoobio_raw_plot <- pcr_raw_zoo_18s_counts %>%
#   ggplot(aes(x = dryweight_C_mg_m2_taxa, y = n_reads_raw, color = cycle.x)) +
#   geom_point(size = 3) +  # Scatter plot
#   # geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#   facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
#   labs(
#     x = "Calanoida  Zooscan Dry Weight (C mg/m²)",
#     y = "Calanoida Raw Reads",
#     color = "Cycle",
#     title = "PCR Reads vs Biomass across Cycles and Size Fractions"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   )
# 
# zoobio_raw_plot
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_reads_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_raw_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_reads.pdf",
#     width = 7, height = 9)
# print(zoobio_raw_plot)
# dev.off()
# 
# #5) Multiply times filter biomass PCR
# zoobio_pcr_plot <- pcr_raw_zoo_18s %>%
#   ggplot(aes(x = dryweight_C_mg_m2_sample, y = n_reads_pcr*biomass_mg_m2, color = cycle.x)) +
#   geom_point(size = 3) +  # Scatter plot
#   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#   # facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
#   labs(
#     x = "Calanoida  Zooscan Dry Weight (C mg/m²)",
#     y = "Calanoida PCR Bias Corrected RA* Filter \n Dry Weight (C mg/m²)",
#     color = "Cycle"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   )
# 
# zoobio_pcr_plot
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_biomass_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_pcr_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_pcr_biomass.pdf",
#     width = 7, height = 9)
# print(zoobio_pcr_plot)
# dev.off()
# 
# 
# 
# #6) Multiply times filter biomass PCR
# zoobio_raw_plot <- pcr_raw_zoo_18s %>%
#   ggplot(aes(x = dryweight_C_mg_m2_taxa, y = n_reads_raw*biomass_mg_m2, color = cycle.x)) +
#   geom_point(size = 3) +  # Scatter plot
#   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
#   # facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
#   labs(
#     x = "Calanoida  Zooscan Dry Weight (C mg/m²)",
#     y = "Calanoida Raw RRA* Filter \n Dry Weight (C mg/m²)",
#     color = "Cycle",
#     title = "PCR Reads vs Biomass across Cycles and Size Fractions"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12),
#     panel.spacing = unit(1, "lines")
#   )
# 
# zoobio_raw_plot
# 
# # Save as PNG
# png("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_biomass_18s.png",
#     width = 7, height = 9, units = "in", res = 300)
# print(zoobio_raw_plot)
# dev.off()
# 
# # Save as PDF
# pdf("Q:/Dante/ZoopMetaB/zooscan_DNA_compare_figures/zooscan_vs_raw_biomass.pdf",
#     width = 7, height = 9)
# print(zoobio_raw_plot)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # #DNA counts
# # pcr_raw_zoo_18s_counts %>% 
# #   ggplot(aes(x = dryweight_C_mg_m2_taxa, y = n_reads_pcr, color=cycle.x)) +
# #   geom_point(size = 3) +  # Scatter plot
# #   # geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.8) +  # Linear regression
# #   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
# #   facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
# #   labs(
# #     x = "Log(Dry Weight (C mg/m²)) per Taxa",
# #     y = "Mean PCR Reads",
# #     color = "Cycle",
# #     title = "PCR Reads vs Biomass across Cycles and Size Fractions"
# #   ) +
# #   theme_minimal() +
# #   theme(
# #     strip.text = element_text(size = 10),
# #     axis.text = element_text(size = 10),
# #     axis.title = element_text(size = 12),
# #     panel.spacing = unit(1, "lines")
# #   )
# # 
# # 
# # 
# # ggplot(pcr_raw_zoo_18s, aes(x = dryweight_C_mg_m2_taxa, y = n_reads_raw, color = size_fraction)) +
# #   geom_point(size = 3) +  # Scatter plot with point size adjustment
# #   geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.8) +  # Regression line
# #   facet_wrap(~ size_fraction, ncol = 1) +  # Facet by size_fraction
# #   # geom_text(
# #   #   data = r2_data,  # Use the R² data for annotation
# #   #   aes(x = Inf, y = Inf, label = paste0("R² = ", round(r_squared, 2))),
# #   #   inherit.aes = FALSE,
# #   #   hjust = 1.2, vjust = 1.2, size = 4, color = "black"
# #   # ) +
# #   labs(
# #     x = "Dry Weight (C mg/m²) per Taxa",
# #     y = "PCR Reads x Biomass (mg/m²)",
# #     color = "Cycle"
# #   ) +
# #   theme_minimal() +
# #   theme(
# #     strip.text = element_text(size = 10),
# #     axis.text = element_text(size = 10),
# #     axis.title = element_text(size = 12),
# #     panel.spacing = unit(1, "lines")
# #   )
# 
# 
# # # Calculate R² using residuals and model.frame
# # r2_data <- pcr_raw_zoo_18s_relab %>%
# #   # group_by(cycle.x) %>%
# #   summarise(
# #     model = list(lm(n_reads_raw ~ count, data = cur_data())),  # Fit the model
# #     .groups = "drop"
# #   ) %>%
# #   mutate(
# #     r_squared = map_dbl(model, ~ {
# #       lm_model <- .x
# #       y_actual <- model.frame(lm_model)$n_reads_raw
# #       y_pred <- fitted(lm_model)
# #       ss_total <- sum((y_actual - mean(y_actual))^2)
# #       ss_residual <- sum((y_actual - y_pred)^2)
# #       1 - (ss_residual / ss_total)  # R² formula
# #     })
# #   )
# # 
# # 
# # # Create the plot
# # pcr_raw_zoo_18s_relab %>% 
# #   ggplot(aes(x = count, y = n_reads_pcr*biomass_mg_m2)) +
# #   geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
# #   geom_smooth(method = "lm", se = TRUE, linetype = "dashed", size = 0.8) +  # Linear regression
# #   geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
# #   geom_text(
# #     data = r2_data,
# #     aes(
# #       x = Inf, y = -Inf,  # Place R² at the bottom-right corner
# #       label = paste0("R² = ", round(r_squared, 2))
# #     ),
# #     inherit.aes = FALSE,
# #     hjust = 1.1, vjust = -1.5, size = 4
# #   ) +
# #   labs(
# #     x = "Log(Dry Weight (C mg/m²)) per Taxa",
# #     y = "Mean PCR Reads",
# #     color = "Cycle",
# #     title = "PCR Reads vs Biomass across Cycles and Size Fractions"
# #   ) +
# #   theme_minimal() +
# #   theme(
# #     strip.text = element_text(size = 10),
# #     axis.text = element_text(size = 10),
# #     axis.title = element_text(size = 12),
# #     panel.spacing = unit(1, "lines")
# #   )

# Scrap -------------------------------------------------------------------

