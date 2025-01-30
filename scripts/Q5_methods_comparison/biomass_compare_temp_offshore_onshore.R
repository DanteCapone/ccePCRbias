library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(patchwork)
library(here)
library(gridExtra)

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


# write.csv(here("Zoop_Patterns/data/zoop_other/biomass_processed.csv"))

  #Add size fraction that will match with Zooscan
  rename(size_fraction_numeric = size_fraction) %>% # Rename the existing column
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 2 ~ "1-2", # Assign 1-5
      size_fraction_numeric > 2 ~ ">2", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )



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
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s1_phy_all_and_subpools_offshore.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s2_phy_all_and_subpools_offshore.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s3_phy_all_and_subpools_offshore.csv")) %>%
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
  # filter(Order=="Calanoida") 
  filter(Family=="Oithonidae")


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
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
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
  select(-asv_code) %>% 
  filter(offshore_onshore=="offshore")


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
phy_18s %>% 
  # filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction, Family) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s

calanoida_dna=phy_18s %>% 
  filter(Order=="Calanoida") %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

eucalanidae_dna=phy_18s %>% 
  filter(Family %in% c("Eucalanidae","Rhincalanidae")) %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

oithona_dna=phy_18s %>% 
  filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Family) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

euphausiid_dna=phy_18s %>% 
  filter(Order=="Euphausiacea") %>%
  # filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)
  

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


pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
  left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



# Zooscan -----------------------------------------------------------------


#Add biomass sum, calanoid biomass and proportion of calanoid biomass
taxa_sel="Calanoida"
taxa_sel="Oithonidae"

#Zooscan biomass proportion: Filter to offshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID_short=Sample_ID) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1,Sample_ID) %>% unique()


zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) %>% 
  left_join(.,clusters, by="Sample_ID_short") %>% 
  filter(offshore_onshore=="offshore")




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
pcr_raw_zoo_18s_calanoida_offshore=zooscan_calanoida %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(calanoida_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Euphausiids
pcr_raw_zoo_18s_euphausiid_offshore=zooscan_euphausiid %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(euphausiid_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Eucalanidae
pcr_raw_zoo_18s_eucalanidae_offshore=zooscan_eucalanus %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(eucalanidae_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Oithonidae
pcr_raw_zoo_18s_oithonidae_offshore=zooscan_oithona %>%
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


# Onshore -----------------------------------------------------------------

# Methods Comparison. ------------------------------------------------------
# Use either proportions 

#Proportions
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s1_phy_all_and_subpools_onshore.csv")) %>%
  select(-X) %>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s2_phy_all_and_subpools_onshore.csv")) %>%
  select(-X)%>% 
  mutate(coord=str_remove(coord, "^clr_"))
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_18s_s3_phy_all_and_subpools_onshore.csv")) %>%
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
  # filter(Order=="Calanoida") 
  filter(Family=="Oithonidae")


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
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
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
#Filter to onshore onshore
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
  select(-asv_code) %>% 
  filter(offshore_onshore=="onshore")


#Join with PCR Proportions
phy_18s %>% 
  # filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction, Family) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s

calanoida_dna=phy_18s %>% 
  filter(Order=="Calanoida") %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

eucalanidae_dna=phy_18s %>% 
  filter(Family %in% c("Eucalanidae","Rhincalanidae")) %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

oithona_dna=phy_18s %>% 
  filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Family) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)

euphausiid_dna=phy_18s %>% 
  filter(Order=="Euphausiacea") %>%
  # filter(Family=="Oithonidae") %>%
  group_by(Sample_ID,size_fraction, Order) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)


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


pcr_and_raw_18s_counts=pcr_and_raw_18s_counts %>% 
  left_join(.,counts_raw_all,by=c("Sample_ID","size_fraction")) %>% 
  mutate(n_reads_pcr_counts=n_reads_pcr*n_reads_raw_sum)



# Zooscan -----------------------------------------------------------------

#Zooscan biomass proportion: Filter to onshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID_short=Sample_ID) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1,Sample_ID) %>% unique()


zooscan_taxa=read.csv(here("PCR_bias_correction/data/Zooscan/zooscan_by_sample_biomass.csv"))%>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) %>% 
  left_join(.,clusters, by="Sample_ID_short")  %>% 
  filter(offshore_onshore=="onshore")


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
pcr_raw_zoo_18s_calanoida_onshore=zooscan_calanoida %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(calanoida_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Euphausiids
pcr_raw_zoo_18s_euphausiid_onshore=zooscan_euphausiid %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(euphausiid_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Eucalanidae
pcr_raw_zoo_18s_eucalanidae_onshore=zooscan_eucalanus %>%
  #remove XL size class
  filter(size_fraction != ">2") %>%
  left_join(eucalanidae_dna, by=c("PC1","size_fraction")) %>% 
  unique(.)


#Oithonidae
pcr_raw_zoo_18s_oithonidae_onshore=zooscan_oithona %>%
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
# Define the plot function
create_plot <- function(data, title) {
  # Calculate limits for the axes
  # x_max <- max(data$dryweight_C_mg_m2_taxa, na.rm = TRUE) * 1.2
  # y_max_pcr <- max(data$n_reads_pcr, na.rm = TRUE) * 1.2
  # y_max_raw <- max(data$n_reads_raw, na.rm = TRUE) * 1.2
  
  # Create the left plot (original)
  plot_left <- data %>%
    ggplot(aes(x = biomass_prop_taxa, y = n_reads_pcr)) +
    geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
    labs(
      title = title,
      x = "Zooscan Dry Weight (C mg/m²)",
      y = "PCR Bias Corrected RRA",
      color = "Cycle"
    ) +
    # xlim(0, x_max) +  # Set x-axis limits
    # ylim(0, y_max_pcr) +  # Set y-axis limits
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.spacing = unit(1, "lines")
    )
  
  # Create the right plot (with n_reads_raw on y-axis)
  plot_right <- data %>%
    ggplot(aes(x = biomass_prop_taxa, y = n_reads_raw)) +
    geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
    labs(
      title = paste(title, "- Raw Reads"),
      x = "Zooscan Dry Weight (C mg/m²)",
      y = "Raw Reads",
      color = "Cycle"
    ) +
    # xlim(0, x_max) +  # Set x-axis limits
    # ylim(0, y_max_raw) +  # Set y-axis limits
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.spacing = unit(1, "lines")
    )
  
  # Combine the two plots side by side using patchwork
  combined_plot <- plot_left + plot_right
  
  return(combined_plot)
}


# 1) Zooscan biomass vs PCR-RRA
# Create individual plots

#Offshore
plot_calanoida_offshore <- create_plot(pcr_raw_zoo_18s_calanoida_offshore, "Calanoida")
plot_calanoida_offshore
plot_euphausiid_offshore <- create_plot(pcr_raw_zoo_18s_euphausiid_offshore, "Euphausiids")
plot_euphausiid_offshore
plot_eucalanidae_offshore <- create_plot(pcr_raw_zoo_18s_eucalanidae_offshore, "Eucalanidae")
plot_eucalanidae_offshore
plot_oithonidae_offshore <- create_plot(pcr_raw_zoo_18s_oithonidae_offshore, "Oithonidae")
plot_oithonidae_offshore

grid.arrange(plot_calanoida_offshore,plot_euphausiid_offshore,plot_eucalanidae_offshore,
             plot_oithonidae_offshore)

#Onshore
plot_calanoida_onshore <- create_plot(pcr_raw_zoo_18s_calanoida_onshore, "Calanoida")
plot_calanoida_onshore
plot_euphausiid_onshore <- create_plot(pcr_raw_zoo_18s_euphausiid_onshore, "Euphausiids")
plot_euphausiid_onshore
plot_eucalanidae_onshore <- create_plot(pcr_raw_zoo_18s_eucalanidae_onshore, "Eucalanidae")
plot_eucalanidae_onshore
plot_oithonidae_onshore <- create_plot(pcr_raw_zoo_18s_oithonidae_onshore, "Oithonidae")
plot_oithonidae_onshore




#Now assess correlations


# Function to calculate R^2 and p-value
assess_correlation <- function(data, x_var, y_var) {
  # Ensure the data contains no NA values for the variables
  data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  
  # Fit linear model
  model <- lm(as.formula(paste(y_var, "~", x_var)), data = data)
  
  # Get model summary
  model_summary <- summary(model)
  
  # Extract R^2 and p-value
  r_squared <- model_summary$r.squared
  p_value <- coef(model_summary)[2,4] # p-value for the slope (2nd row in coefficients table)
  
  return(list(r_squared = r_squared, p_value = p_value))
}

# Compare correlations for each dataset
datasets <- list(
  Calanoida = pcr_raw_zoo_18s_calanoida_onshore,
  Euphausiids = pcr_raw_zoo_18s_euphausiid_onshore,
  Eucalanidae = pcr_raw_zoo_18s_eucalanidae_onshore,
  Oithonidae = pcr_raw_zoo_18s_oithonidae_onshore
)


results <- lapply(names(datasets), function(name) {
  data <- datasets[[name]]
  
  # Assess correlation for n_reads_pcr
  pcr_results <- assess_correlation(data, "biomass_prop_taxa", "n_reads_pcr")
  
  # Assess correlation for n_reads_raw
  raw_results <- assess_correlation(data, "biomass_prop_taxa", "n_reads_raw")
  
  # Combine results for this dataset
  return(data.frame(
    Dataset = name,
    R2_PCR = pcr_results$r_squared,
    P_PCR = pcr_results$p_value,
    R2_Raw = raw_results$r_squared,
    P_Raw = raw_results$p_value
  ))
})

# Combine all results into one dataframe
correlation_results <- do.call(rbind, results)

# Display results
print(correlation_results)


## Combine all plots with R2
# Function to create a single plot and add R^2 and p-value as text
create_single_plot <- function(data, x_var, y_var, title) {
  # Assess correlation
  correlation <- assess_correlation(data, x_var, y_var)
  r_squared <- round(correlation$r_squared, 2)
  p_value <- signif(correlation$p_value, 3)
  
  # Create the plot
  data %>%
    ggplot(aes_string(x = x_var, y = y_var)) +
    geom_point(aes(color = cycle.x), size = 3) +  # Scatter plot
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +  # 1:1 reference line
    labs(
      title = title,
      x = "Zooscan Proportion Biomass",
      y = ifelse(y_var == "n_reads_pcr", "PCR Bias Corrected \nRelative Read Abundance", "Raw Relative Read Abundance"),
      color = "Cycle"
    ) +
    annotate("text", x = Inf, y = Inf, label = paste0("R²: ", r_squared, "\nP: ", p_value),
             hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.spacing = unit(1, "lines")
    )
}

# Create plots for offshore datasets
plot_calanoida_offshore_pcr <- create_single_plot(pcr_raw_zoo_18s_calanoida_offshore, "biomass_prop_taxa", "n_reads_pcr", "Calanoida Offshore PCR")
plot_calanoida_offshore_raw <- create_single_plot(pcr_raw_zoo_18s_calanoida_offshore, "biomass_prop_taxa", "n_reads_raw", "Calanoida Offshore Raw")

# plot_euphausiid_offshore_pcr <- create_single_plot(pcr_raw_zoo_18s_euphausiid_offshore, "biomass_prop_taxa", "n_reads_pcr", "Euphausiids Offshore PCR")
# plot_euphausiid_offshore_raw <- create_single_plot(pcr_raw_zoo_18s_euphausiid_offshore, "biomass_prop_taxa", "n_reads_raw", "Euphausiids Offshore Raw")
# 
# plot_eucalanidae_offshore_pcr <- create_single_plot(pcr_raw_zoo_18s_eucalanidae_offshore, "biomass_prop_taxa", "n_reads_pcr", "Eucalanidae Offshore PCR")
# plot_eucalanidae_offshore_raw <- create_single_plot(pcr_raw_zoo_18s_eucalanidae_offshore, "biomass_prop_taxa", "n_reads_raw", "Eucalanidae Offshore Raw")
# 
# plot_oithonidae_offshore_pcr <- create_single_plot(pcr_raw_zoo_18s_oithonidae_offshore, "biomass_prop_taxa", "n_reads_pcr", "Oithonidae Offshore PCR")
# plot_oithonidae_offshore_raw <- create_single_plot(pcr_raw_zoo_18s_oithonidae_offshore, "biomass_prop_taxa", "n_reads_raw", "Oithonidae Offshore Raw")

# Create plots for onshore datasets
plot_calanoida_onshore_pcr <- create_single_plot(pcr_raw_zoo_18s_calanoida_onshore, "biomass_prop_taxa", "n_reads_pcr", "Calanoida Onshore PCR")
plot_calanoida_onshore_raw <- create_single_plot(pcr_raw_zoo_18s_calanoida_onshore, "biomass_prop_taxa", "n_reads_raw", "Calanoida Onshore Raw")

# plot_euphausiid_onshore_pcr <- create_single_plot(pcr_raw_zoo_18s_euphausiid_onshore, "biomass_prop_taxa", "n_reads_pcr", "Euphausiids Onshore PCR")
# plot_euphausiid_onshore_raw <- create_single_plot(pcr_raw_zoo_18s_euphausiid_onshore, "biomass_prop_taxa", "n_reads_raw", "Euphausiids Onshore Raw")
# 
# plot_eucalanidae_onshore_pcr <- create_single_plot(pcr_raw_zoo_18s_eucalanidae_onshore, "biomass_prop_taxa", "n_reads_pcr", "Eucalanidae Onshore PCR")
# plot_eucalanidae_onshore_raw <- create_single_plot(pcr_raw_zoo_18s_eucalanidae_onshore, "biomass_prop_taxa", "n_reads_raw", "Eucalanidae Onshore Raw")
# 
# plot_oithonidae_onshore_pcr <- create_single_plot(pcr_raw_zoo_18s_oithonidae_onshore, "biomass_prop_taxa", "n_reads_pcr", "Oithonidae Onshore PCR")
# plot_oithonidae_onshore_raw <- create_single_plot(pcr_raw_zoo_18s_oithonidae_onshore, "biomass_prop_taxa", "n_reads_raw", "Oithonidae Onshore Raw")

# Arrange the plots in a 2x8 grid
grid.arrange(
  plot_calanoida_offshore_pcr,
  plot_calanoida_onshore_pcr, 
  plot_calanoida_offshore_raw,
  plot_calanoida_onshore_raw
)
