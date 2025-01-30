#This script compares PCR-RA, RRA, and Zoo-PB
#First performs correlations between methods
# Analysis is conducted for both 18S and COI primers


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

metadata=metadata %>% 
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x) %>% 
  column_to_rownames("Sample_ID")




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
  left_join(.,metadata, by="Sample_ID")%>%
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
  left_join(.,metadata, by="Sample_ID")%>%
  mutate(taxa = coord)

#Filter to calanoida
taxa_pcr=phy_taxa_pcr %>% mutate(Family=taxa) %>%
  left_join(.,zhan_taxa %>% rownames_to_column("Family"), by="Family") %>%
  filter(Order=="Calanoida") 


#PCR-RA df ready to join with RRA
pcr_join_prop=taxa_pcr %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  group_by(Sample_ID,size_fraction,PC1,cycle) %>%
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
env_metadata_phy=metadata %>%
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
  select(-asv_code) %>% 
  mutate(Sample_ID=Sample_ID.x) 


#Raw in Proportions
phy_18s=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  select(-asv_code)%>%
  #Add size fraction that will match with Zooscan
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  ) %>% 
  mutate(Sample_ID=Sample_ID.x) 


#Raw in CLR
phy_18s_clr=phyloseq(OTU, TAX, meta) %>% 
  tax_glom(taxrank="Order", NArm=TRUE) %>%
  transform_sample_counts(., clr_convert)%>%
  phyloseq_transform_to_long(.)%>%
  #Add size fraction that will match with Zooscan
  mutate(
    size_fraction = case_when(
      size_fraction_numeric %in% c(0.2, 0.5) ~ "0.2-1", # Combine 0.2 and 0.5
      size_fraction_numeric >= 1 & size_fraction_numeric <= 5 ~ "1-5", # Assign 1-5
      TRUE ~ NA_character_ # Handle any unexpected values
    )
  )

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
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_18s

#Join in counts
phy_18s_raw_counts %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction) %>%
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
  left_join(.,metadata, by=c("PC1","size_fraction"))  

# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category=="Calanoida")  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) 


#Zooscan relative abundance CLR
zooscan_relative_clr=zooscan_relative %>% 
  group_by(sample_id,size_fraction) %>% 
  mutate(abundance_clr=clr_convert(count))%>%
  ungroup() %>% 
  group_by(sample_id) %>% 
  left_join(.,metadata, by=c("PC1","size_fraction"))


#========== COMPARE: Make combined dataframe for comparing all 3 methods

#Propotions
pcr_raw_zoo_18s=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">5") %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction")) %>% 
  unique(.)

#Counts
pcr_raw_zoo_18s_counts=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != ">5") %>%
  left_join(pcr_and_raw_18s_counts, by=c("PC1","size_fraction")) %>% 
  unique(.)


pcr_raw_zoo_18s_clr=zooscan_taxa_clr %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_18s_clr, by=c("PC1","size_fraction"))

pcr_raw_zoo_18s_relab=zooscan_relative %>%
  #remove XL size class
  filter(size_fraction != ">5") %>%
  left_join(pcr_and_raw_18s, by=c("PC1","size_fraction")) %>% 
  unique(.)

pcr_raw_zoo_18s_relab_clr=zooscan_relative_clr %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_18s_clr, by=c("PC1","size_fraction"))





# First: Proportions Analysis ----------------------------------------------------


#=============Grouped Barplot site abundances
#Format Long and add difference metrics
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s, 
                                     cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_18s %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_18s_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
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

saving=1
if (saving==1) {
    grouped_bar_all_18s=grouped_bar_all_18s+
      labs(title = "Methods Differences Biomass Proportion 18S",y = "Proportion Reads or Biomass")
    
    ggsave(
      filename = here("PCR_bias_correction/figures/grouped_bar_relative_abundance_diff_18s_biomass.pdf"),
      # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
      plot = grouped_bar_all_18s,
      width = 8,  # Width in inches
      height = 6  # Height in inches
    )
    
    ggsave(
      filename = here("PCR_bias_correction/figures/grouped_bar_relative_abundance_diff_18s_biomass.png"),
      # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
      plot = grouped_bar_all_18s,
      width = 8,  # Width in inches
      height = 6  # Height in inches
    )
  }
  
  
#2. Relative abundance by counts
#Format Long and add difference metrics
pcr_raw_zoo_18s_long <- pivot_longer(pcr_raw_zoo_18s_relab, 
                                     cols = c(relative_abundance, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "relative_abundance" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "relative_abundance" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "relative_abundance" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "relative_abundance" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_18s_relab %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_18s_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences 18S",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Abundance",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#4F86F7", "#F78D4F","#70BF41"),
                    labels = c("PCR-corrected", "Raw Reads","Zooscan Abundances")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+# Define linetypes for each Method
  scale_color_manual(values = c("#4F86F7", "#F78D4F","#70BF41")) -> grouped_bar_all_18s  # Define colors for each Method

grouped_bar_all_18s

saving=1
if (saving==1) {
  grouped_bar_all_18s=grouped_bar_all_18s+
    labs(title = "Methods Differences Biomass Proportion 18S",y = "Proportion Reads or Abundance")
  
  ggsave(
    filename = here("PCR_bias_correction/figures/grouped_bar_relative_abundance_diff_18s_counts.pdf"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_18s,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("PCR_bias_correction/figures/grouped_bar_relative_abundance_diff_18s_counts.png"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
    plot = grouped_bar_all_18s,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
}



# Error metrics

# Mean squared error ---------------------------------------------------------------------
pcr_raw_zoo_18s_mse=pcr_raw_zoo_18s%>% 
  # filter(!is.na(n_reads_raw)) %>% 
  group_by(size_fraction) %>% 
  mutate(se_zoo_pcr = (log(biomass_prop_taxa) - log(n_reads_pcr))^2,
         se_zoo_raw = (log(biomass_prop_taxa) - log(n_reads_raw))^2,
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
  labs(title = "MedSE Differences between Methods",
       x = "PC1",
       y = "Difference in Square-Error (RRA-PCR-RA)",
       fill="") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none")  +
  geom_text(data = sum_mse, aes(x = 0.5, y = -0.2, label = sprintf("MedSE PCR-RA = %.3f\nMedSE RRA = %.3f", average_mse_by_size_pcr, average_mse_by_size_raw)),
            hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE)->mse_plot
  
mse_plot
saving=1
if (saving==1) {
    
    ggsave(
      filename = here("plots/methods_comparison/error_compare/mse_plot_18s.pdf"),
      plot = mse_plot,
      width = 9,  # Width in inches
      height = 6  # Height in inches
    )
    
    ggsave(
      filename = here("plots/methods_comparison/error_compare/mse_plot_18s.png"),
      plot = mse_plot,
      width = 9,  # Width in inches
      height = 6  # Height in inches
    )
  }
  
  

# Using Relative Abundances from Zooscan
#Format Long and add difference metrics
pcr_raw_zoo_18s_long_relab <- pivot_longer(pcr_raw_zoo_18s_relab, 
                                     cols = c(relative_abundance, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "relative_abundance" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "relative_abundance" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "relative_abundance" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "relative_abundance" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_18s_relab %>% 
  ungroup()%>%
  select(Sample_ID_short.x,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_18s_long_relab %>%
  # filter(Method=="biomass_prop_taxa") %>% 
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Counts",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#4F86F7", "#F78D4F","#70BF41"),
                    labels = c( "PCR-corrected", "Raw Reads","Zooscan Relative \nAbundance")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+# Define linetypes for each Method
  scale_color_manual(values = c( "#4F86F7", "#F78D4F","#70BF41")) -> grouped_bar_all_18s_relab  # Define colors for each Method

grouped_bar_all_18s_relab

saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_relab.pdf"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_18s_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_relab.png"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_18s_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
}




# Correlations: Use CLR ---------------------------------------------------


# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test_18s_pcr <- shapiro.test(pcr_raw_zoo_18s_clr$n_reads_pcr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test_18s_pcr)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = n_reads_pcr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = n_reads_pcr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_clr$n_reads_raw)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = n_reads_raw)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = n_reads_raw)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_clr$biomass_taxa_clr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = biomass_taxa_clr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = biomass_taxa_clr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()
pcr_raw_zoo_18s_clr=pcr_raw_zoo_18s_clr %>%
  left_join(.,clusters,by="PC1")



## Plotting
# Zooscan vs. PCR-RA
pcr_raw_zoo_18s_clr %>%
  filter(!is.na(cycle.x.x))%>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr_clr
zoo_vs_pcr_clr


if (saving==1) {
ggsave(
  filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clr.pdf"),
  plot = zoo_vs_pcr_clr,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### =====  Zooscan vs RRA
pcr_raw_zoo_18s_clr %>%
  filter(!is.na(cycle.x.x))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((n_reads_pcr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  facet_wrap(~size_fraction, nrow=3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 0.8, label.y = 0.8)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)


#. Relative Abundance (CLR) Correlations --------------------------------


# Check distribution & test for normaility
#Zoo-RA 
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_relab_clr$abundance_clr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_relab_clr %>%
  ggplot(aes(x = abundance_clr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_relab_clr %>%
  ggplot(aes(sample = abundance_clr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_18s_relab_clr=pcr_raw_zoo_18s_relab_clr %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_18s_relab %>%
  filter(!is.na(cycle.y))%>%
  ggplot(.,aes(x=((clr_convert(count))), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr_clr
zoo_vs_pcr_clr


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

# Zooscan Relative Abundance vs RRA
pcr_raw_zoo_18s_relab_clr %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((abundance_clr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((n_reads_pcr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  # facet_wrap(~size_fraction, nrow=3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 0.8, label.y = 0.8)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/pcr_vs_raw_correlation.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)



# ## Correlations: Using Proportion data ------------------------------------

# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_18s$n_reads_pcr)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(n_reads_pcr)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_18s$n_reads_raw)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(n_reads_raw)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test_zoo <- shapiro.test(pcr_raw_zoo_18s$biomass_prop_taxa)

# Print the result of the Shapiro-Wilk test
print(shapiro_test_zoo)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(biomass_prop_taxa)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_18s=pcr_raw_zoo_18s %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  ggplot(.,aes(x=asin(sqrt(biomass_prop_taxa)), y=asin(sqrt(n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
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


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_pcr_correlation_18s.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### 18S =====  
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.x.x))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_taxa)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(n_reads_pcr)), y=asin(sqrt(n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance\n (arcsine square-root)", y = "Raw Reads Relative \n (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  facet_wrap(~size_fraction, nrow=3) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 1.2, label.y = 1.2)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s.pdf"),
  plot = pcr_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s.png"),
  plot = pcr_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)


# PART II: COI ---------------------------------------------------------------------


#Predicted proportions
fido_s1=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s1_phy_all_and_subpools.csv")) %>%
  select(-X)
fido_s2=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s2_phy_all_and_subpools.csv")) %>%
  select(-X)
fido_s3=read.csv(here("PCR_bias_correction/data/predicted_og/predicted_og_coi_05_29_2024_s3_phy_all_and_subpools.csv")) %>%
  select(-X)
#Merge
final_data_all_sizes=rbind(fido_s1,fido_s2,fido_s3) %>%
  mutate(Sample_ID = str_extract(replicate, "(?<=predicted )\\S+")) 


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
  left_join(.,env_metadata, by="Sample_ID") %>%
  mutate(taxa=coord)



# Raw Reads ---------------------------------------------------------------


## === RAW READS: Repeat Analysis with raw/normalized reads === ##

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

#Metadata
metadata=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  select(-X, -Sizefractionmm,max_size) %>%
  mutate(Sample_ID=Sample_ID_dot) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1)

env_metadata_phy=metadata %>% 
  left_join(.,dryweights, by = c("Sample_ID_short","max_size"))%>% 
  left_join(.,depths, by="Sample_ID_short") %>%
  mutate(biomass_dry = replace(biomass_dry, which(biomass_dry<0), NA)) %>%
  mutate(biomass_mg_m2=biomass_dry/Volume_Filtered_m3*210) %>%
  select(-Sample_ID.y) %>%
  mutate(Sample_ID=Sample_ID.x) %>% 
  filter(Sample_ID %in% colnames(fido_coi_merged_raw)) %>% 
  column_to_rownames("Sample_ID")


#Make phyloseq objects

#coi
OTU = otu_table(as.matrix(fido_coi_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(coi_taxa))
meta=sample_data(env_metadata_phy)
Phy_raw_coi=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.)

Phy_props_coi=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Genus=asv_code) %>%
  select(-asv_code)



# PCR vs Raw --------------------------------------------------------------

#Raw in Proportions
phy_coi=transform_sample_counts(phyloseq(OTU, TAX, meta), function(x) x / sum(x))%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Family=asv_code) %>%
  mutate(Sample_ID=file_code) %>% 
  select(-asv_code, -file_code) 


#Raw in CLR
phy_coi_clr=phyloseq(OTU, TAX, meta) %>% 
  tax_glom(taxrank="Order", NArm=TRUE) %>%
  transform_sample_counts(., clr_convert)%>%
  phyloseq_transform_to_long(.)%>%
  mutate(Family=asv_code) %>%
  mutate(Sample_ID=file_code) %>% 
  select(-asv_code, -file_code) 

#filter to calanoid copepods
taxa_sel="Calanoida"


#Join with PCR 
phy_coi_clr %>% 
  filter(Order==taxa_sel) %>%
  left_join(pcr_join_clr, by=c("Sample_ID","PC1"), keep = FALSE) %>%
  select(-size_fraction.x) %>%
  mutate(n_reads_raw=n_reads) %>% 
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_coi_clr

#Join with Raw
phy_coi %>% 
  filter(Order==taxa_sel) %>%
  group_by(Sample_ID,size_fraction) %>%
  summarise(n_reads_raw=sum(n_reads), biomass_mg_m2=mean(biomass_mg_m2)) %>%  
  left_join(pcr_join_prop, by="Sample_ID") %>%
  select(-size_fraction.x) %>%
  mutate(size_fraction=size_fraction.y)->pcr_and_raw_coi





# Zooscan -----------------------------------------------------------------


#===== Zooscan comparison

#Need to modify string category for joining
size_mapping <- c("0.2-0.5" = 0.2, "0.5-1" = 0.5, "1-2" = 1, ">2" = 5)


#Read in processed Zooscan data and look at biomass
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


# Zooscan relative abundances
zooscan_relative=read.csv(here("PCR_bias_correction/data/Zooscan/zoop_relative_abundance.csv"))%>%
  filter(object_annotation_category=="Calanoida")  %>%
  select(-X) %>% 
  mutate(Sample_ID=sample_id) 




#========== COMPARE: Make combined dataframe for comparing all 3 methods

#Propotions
pcr_raw_zoo_coi=zooscan_taxa %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_coi, by=c("PC1","size_fraction"))

pcr_raw_zoo_coi_relab=zooscan_relative %>%
  #remove XL size class
  filter(size_fraction != 5) %>%
  left_join(pcr_and_raw_coi, by=c("PC1","size_fraction"))



# First: Proportions Analysis ----------------------------------------------------


#=============Grouped Barplot site abundances
#Format Long and add difference metrics
pcr_raw_zoo_coi_long <- pivot_longer(pcr_raw_zoo_coi, 
                                     cols = c(biomass_prop_taxa, n_reads_raw, n_reads_pcr), 
                                     names_to = "Method", 
                                     values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "biomass_prop_taxa" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "biomass_prop_taxa" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "biomass_prop_taxa" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "biomass_prop_taxa" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_coi %>% 
  ungroup()%>%
  select(Sample_ID_short,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_coi_long %>%
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
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
  scale_color_manual(values = c("#70BF41", "#4F86F7", "#F78D4F")) -> grouped_bar_all_coi  # Define colors for each Method

grouped_bar_all_coi

saving=1
if (saving==1) {
  grouped_bar_all_coi=grouped_bar_all_coi+
    labs(title = "Methods Differences Biomass Proportion COI",y = "Proportion Reads or Biomass")
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi_biomass.pdf"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_coi,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi_biomass.png"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.png"),
    plot = grouped_bar_all_coi,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
}





# Error metrics

# Mean squared error ---------------------------------------------------------------------
pcr_raw_zoo_coi_mse=pcr_raw_zoo_coi%>% 
  # filter(!is.na(n_reads_raw)) %>% 
  group_by(size_fraction) %>% 
  mutate(se_zoo_pcr = (log(biomass_prop_taxa) - log(n_reads_pcr))^2,
         se_zoo_raw = (log(biomass_prop_taxa) - log(n_reads_raw))^2,
         average_se = (se_zoo_pcr + se_zoo_raw) / 2)%>%
  mutate(difference = se_zoo_raw - se_zoo_pcr,
         colorr = ifelse(difference >= 0, "#f5776e", "#8d9af2")) %>% 
  ungroup()

hist(pcr_raw_zoo_coi_mse$se_zoo_pcr)
hist(pcr_raw_zoo_coi_mse$se_zoo_raw)

# Calculate average MSE per size
sum_mse <- pcr_raw_zoo_coi_mse %>%
  group_by(size_fraction) %>%
  summarise(
    average_mse_by_size_pcr = median(se_zoo_pcr, na.rm = TRUE),
    average_mse_by_size_raw = median(se_zoo_raw, na.rm = TRUE)
  )

sum_mse

# Error barplot
ggplot(pcr_raw_zoo_coi_mse, aes(x = as.factor(PC1), y = difference, fill = colorr)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("#f5776e" = "#f5776e", "#8d9af2" = "#8d9af2")) +
  theme_minimal() +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y", 
             labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  labs(title = "MedSE Differences between Methods COI",
       x = "PC1",
       y = "Difference in Square-Error (RRA-PCR-RA)",
       fill="") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none")  +
  geom_text(data = sum_mse, aes(x = 0.5, y = -0.2, label = sprintf("MedSE PCR-RA = %.3f\nMedSE RRA = %.3f", average_mse_by_size_pcr, average_mse_by_size_raw)),
            hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE)->mse_plot

mse_plot
saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi.pdf"),
    plot = mse_plot,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi.png"),
    plot = mse_plot,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
}



# Using Relative Abundances from Zooscan

# Relative Aundance Analysis (Zoo-RA) -------------------------------------



# 18S relative abundances -------------------------------------------------
# Using Relative Abundances from Zooscan
#Format Long and add difference metrics
pcr_raw_zoo_18s_long_relab <- pivot_longer(pcr_raw_zoo_18s_relab, 
                                           cols = c(relative_abundance, n_reads_raw, n_reads_pcr), 
                                           names_to = "Method", 
                                           values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "relative_abundance" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "relative_abundance" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "relative_abundance" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "relative_abundance" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_18s_relab %>% 
  ungroup()%>%
  select(Sample_ID_short.x,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_18s_long_relab %>%
  # filter(Method=="biomass_prop_taxa") %>% 
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Counts",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#4F86F7", "#F78D4F","#70BF41"),
                    labels = c( "PCR-corrected", "Raw Reads","Zooscan Relative \nAbundance")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+# Define linetypes for each Method
  scale_color_manual(values = c( "#4F86F7", "#F78D4F","#70BF41")) -> grouped_bar_all_18s_relab  # Define colors for each Method

grouped_bar_all_18s_relab

saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_relab.pdf"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_18s_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_relab.png"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_18s_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
}

# Median squared error ---------------------------------------------------------------------
pcr_raw_zoo_18s_relab_mse=pcr_raw_zoo_18s_relab%>% 
  # filter(!is.na(n_reads_raw)) %>% 
  group_by(size_fraction) %>% 
  mutate(se_zoo_pcr = (log(relative_abundance) - log(n_reads_pcr))^2,
         se_zoo_raw = (log(relative_abundance) - log(n_reads_raw))^2,
         average_se = (se_zoo_pcr + se_zoo_raw) / 2)%>%
  mutate(difference = se_zoo_raw - se_zoo_pcr,
         colorr = ifelse(difference >= 0, "#f5776e", "#8d9af2")) %>% 
  ungroup()

hist(pcr_raw_zoo_18s_mse$se_zoo_pcr)
hist(pcr_raw_zoo_18s_mse$se_zoo_raw)

# Calculate average MSE per size
sum_mse <- pcr_raw_zoo_18s_relab_mse %>%
  group_by(size_fraction) %>%
  summarise(
    average_mse_by_size_pcr = median(se_zoo_pcr, na.rm = TRUE),
    average_mse_by_size_raw = median(se_zoo_raw, na.rm = TRUE)
  )

sum_mse

# Error barplot
ggplot(pcr_raw_zoo_18s_relab_mse, aes(x = as.factor(PC1), y = difference, fill = colorr)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("#f5776e" = "#f5776e", "#8d9af2" = "#8d9af2")) +
  theme_minimal() +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y", 
             labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  labs(title = "MedSE Differences between Methods",
       x = "PC1",
       y = "Difference in Square-Error (RRA-PCR-RA)",
       fill="") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none")  +
  geom_text(data = sum_mse, aes(x = 0.5, y = -0.2, label = sprintf("MedSE PCR-RA = %.3f\nMedSE RRA = %.3f", average_mse_by_size_pcr, average_mse_by_size_raw)),
            hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE)->mse_plot_18s_relab

mse_plot_18s_relab
saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_18s_relab.pdf"),
    plot = mse_plot_18s_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_18s_relab.png"),
    plot = mse_plot_18s_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
}


# COI relative abundances -------------------------------------------------


#Format Long and add difference metrics
pcr_raw_zoo_coi_long_relab <- pivot_longer(pcr_raw_zoo_coi_relab, 
                                           cols = c(relative_abundance, n_reads_raw, n_reads_pcr), 
                                           names_to = "Method", 
                                           values_to = "relative_abundance") %>%
  filter(!is.na(Sample_ID.y)) %>%
  group_by(Sample_ID.x, size_fraction) %>%
  mutate(diff_pcr = abs(relative_abundance[Method == "n_reads_pcr"] - relative_abundance),
         diff_raw = abs(relative_abundance[Method == "n_reads_raw"] - relative_abundance))%>%
  mutate(is_closer = ifelse(Method == "relative_abundance" & diff_pcr < diff_raw, "*", 
                            ifelse(Method == "relative_abundance" & diff_pcr > diff_raw, "x", NA))) %>%
  mutate(closest = ifelse(Method == "relative_abundance" & diff_pcr < 0.05*relative_abundance, "**",NA)) %>%
  mutate(worse = ifelse(Method == "relative_abundance" & diff_raw < 0.05*relative_abundance, "xx",NA))





#Plot grouped bar plot for RRA, PCR-RA, Zoo-PB proporitions
labels_for_map=pcr_raw_zoo_coi_relab %>% 
  ungroup()%>%
  select(Sample_ID_short.x,PC1) %>%
  unique(.) %>%
  arrange((PC1))

pcr_raw_zoo_coi_long_relab %>%
  # filter(Method=="biomass_prop_taxa") %>% 
  ggplot(., aes(x = as.factor(PC1), y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~size_fraction, nrow=3, scale="free_y", labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  theme_minimal() +
  labs(title = "Methods Differences",
       x = expression(paste("Offshore ", PC1, " Onshore")),
       y = "Proportion Reads or Counts",
       fill = "Method") +
  coord_cartesian(ylim=c(0,1.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_text(aes(label = is_closer), position = position_dodge(width = 0.9), vjust = -1, size = 5) +
  scale_fill_manual(values = c("#4F86F7", "#F78D4F","#70BF41"),
                    labels = c( "PCR-corrected", "Raw Reads","Zooscan Relative \nAbundance")) +
  # geom_errorbar(aes(ymin = p.2.5, ymax = p.97.5), position = position_dodge(width = 0.9), width = 0.25) +
  scale_shape_manual(values = c(16, 17, 18)) +  # Define shapes for each Method
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+# Define linetypes for each Method
  scale_color_manual(values = c( "#4F86F7", "#F78D4F","#70BF41")) -> grouped_bar_all_coi_relab  # Define colors for each Method

grouped_bar_all_coi_relab

saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi_relab.pdf"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_coi_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_coi_relab.png"),
    # filename = here("plots/methods_comparison/line_relative_abundance_diff.pdf"),
    plot = grouped_bar_all_coi_relab,
    width = 8,  # Width in inches
    height = 6  # Height in inches
  )
}

mse_plot_coi_18s_relab_bar=grid.arrange(grouped_bar_all_18s_relab+
                                      labs(title = "18S") +
                                        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)), grouped_bar_all_coi_relab+
                                      labs(title = "COI") +
                                      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24)), ncol=2, widths = c(1, 1))

saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_coi_relab.png"),
    plot = mse_plot_coi_18s_relab_bar,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/grouped_bar_relative_abundance_diff_18s_coi_relab.pdf"),
    plot = mse_plot_coi_18s_relab_bar,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
}

# Median squared error ---------------------------------------------------------------------
pcr_raw_zoo_coi_relab_mse=pcr_raw_zoo_coi_relab%>% 
  # filter(!is.na(n_reads_raw)) %>% 
  group_by(size_fraction) %>% 
  mutate(se_zoo_pcr = (log(relative_abundance) - log(n_reads_pcr))^2,
         se_zoo_raw = (log(relative_abundance) - log(n_reads_raw))^2,
         average_se = (se_zoo_pcr + se_zoo_raw) / 2)%>%
  mutate(difference = se_zoo_raw - se_zoo_pcr,
         colorr = ifelse(difference >= 0, "#f5776e", "#8d9af2")) %>% 
  ungroup()

hist(pcr_raw_zoo_coi_mse$se_zoo_pcr)
hist(pcr_raw_zoo_coi_mse$se_zoo_raw)

# Calculate average MSE per size
sum_mse <- pcr_raw_zoo_coi_relab_mse %>%
  group_by(size_fraction) %>%
  summarise(
    average_mse_by_size_pcr = median(se_zoo_pcr, na.rm = TRUE),
    average_mse_by_size_raw = median(se_zoo_raw, na.rm = TRUE)
  )

sum_mse

# Error barplot
ggplot(pcr_raw_zoo_coi_relab_mse, aes(x = as.factor(PC1), y = difference, fill = colorr)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("#f5776e" = "#f5776e", "#8d9af2" = "#8d9af2")) +
  theme_minimal() +
  facet_wrap(~size_fraction, nrow = 3, scales = "free_y", 
             labeller = label_bquote(rows = .(c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")))) +
  labs(title = "MedSE Differences between Methods",
       x = "PC1",
       y = "Difference in Square-Error (RRA-PCR-RA)",
       fill="") +
  scale_x_discrete(labels = labels_for_map$Sample_ID_short.x)+
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position = "none")  +
  geom_text(data = sum_mse, aes(x = 0.5, y = -0.2, label = sprintf("MedSE PCR-RA = %.3f\nMedSE RRA = %.3f", average_mse_by_size_pcr, average_mse_by_size_raw)),
            hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE)->mse_plot_coi_relab

mse_plot_coi_relab
saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi_relab.pdf"),
    plot = mse_plot_coi_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi_relab.png"),
    plot = mse_plot_coi_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
}


mse_plot_coi_18s_relab=grid.arrange(mse_plot_18s_relab+
               labs(title = "18S") +
               theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24)), mse_plot_coi_relab+
               labs(title = "COI") +
               theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24)), ncol=2)

saving=1
if (saving==1) {
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi_and_18s_relab.pdf"),
    plot = mse_plot_coi_18s_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/error_compare/mse_plot_coi_and_18s_relab.png"),
    plot = mse_plot_coi_18s_relab,
    width = 9,  # Width in inches
    height = 6  # Height in inches
  )
}


# COI Correlations ------------------------------------------------------------

#COI
# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_coi_clr$n_reads_pcr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(x = n_reads_pcr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(sample = n_reads_pcr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_coi_clr$n_reads_raw)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(x = n_reads_raw)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(sample = n_reads_raw)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_coi_clr$biomass_taxa_clr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(x = biomass_taxa_clr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi_clr %>%
  ggplot(aes(sample = biomass_taxa_clr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_coi_clr=pcr_raw_zoo_coi_clr %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_coi_clr %>%
  filter(!is.na(cycle.x.x))%>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr_clr
zoo_vs_pcr_clr


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_coi_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_coi_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### =====  
pcr_raw_zoo_coi_clr %>%
  filter(!is.na(cycle.x.x))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_coi_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_coi_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=asin(sqrt(n_reads_pcr)), y=asin(sqrt(n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance\n (arcsine square-root)", y = "Raw Reads Relative Abundance\n (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance COI")+
  facet_wrap(~size_fraction, nrow=3) +
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 1.2, label.y = 1.25)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_coi.pdf"),
  plot = pcr_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_coi.png"),
  plot = pcr_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)




# ## Correlations: Use Proportion data ------------------------------------

# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_coi$n_reads_pcr)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_coi %>%
  ggplot(aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi %>%
  ggplot(aes(sample = asin(sqrt(n_reads_raw)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_coi$n_reads_raw)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_coi %>%
  ggplot(aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi %>%
  ggplot(aes(sample = asin(sqrt(n_reads_raw)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test_zoo <- shapiro.test(pcr_raw_zoo_coi$biomass_prop_taxa)

# Print the result of the Shapiro-Wilk test
print(shapiro_test_zoo)

# Plot the histogram
pcr_raw_zoo_coi %>%
  ggplot(aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_coi %>%
  ggplot(aes(sample = asin(sqrt(biomass_prop_taxa)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_coi=pcr_raw_zoo_coi %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  ggplot(.,aes(x=((biomass_prop_taxa)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
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


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/correlations/zooscan_vs_pcr_correlation_coi.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_coi.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### =====  
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_prop_taxa)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_coi_non_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_coi_non_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_coi %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((n_reads_pcr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  facet_wrap(~size_fraction, nrow=3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 0.8, label.y = 0.8)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_coi_non_transformed.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_coi_non_transformed.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)




# Correlations 18S --------------------------------------------------------

#18s
# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_clr$n_reads_pcr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = n_reads_pcr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = n_reads_pcr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_clr$n_reads_raw)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = n_reads_raw)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = n_reads_raw)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(pcr_raw_zoo_18s_clr$biomass_taxa_clr)

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(x = biomass_taxa_clr)) +   # Set the data and the variable to plot
  geom_histogram(binwidth = 1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s_clr %>%
  ggplot(aes(sample = biomass_taxa_clr)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_18s_clr=pcr_raw_zoo_18s_clr %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_18s_clr %>%
  filter(!is.na(cycle.x.x))%>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  facet_wrap(~size_fraction, nrow=3) +
  labs(x = "Zooscan Biomass (CLR)", y = "PCR Bias-Mitigated Read Abundance (CLR)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("Pearson Correlation between Zooscan Biomass and PCR Bias-Mitigated Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.2, label.y =.5)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_pcr_clr
zoo_vs_pcr_clr


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### =====  
pcr_raw_zoo_18s_clr %>%
  filter(!is.na(cycle.x.x))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_taxa_clr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.x.x, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_18s_clr.pdf"),
    plot = zoo_vs_pcr_clr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((n_reads_pcr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  facet_wrap(~size_fraction, nrow=3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 0.8, label.y = 0.8)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s_clr.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s_clr.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)




# ## Correlations: Use Proportion data ------------------------------------

# Check distribution & test for normaility

#1) PCR RA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_18s$n_reads_pcr)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(n_reads_pcr)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(n_reads_raw)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")

#1) RRA
# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(asin(sqrt(pcr_raw_zoo_18s$n_reads_raw)))

# Print the result of the Shapiro-Wilk test
print(shapiro_test)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(n_reads_raw)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(n_reads_raw)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")


#3) Zoo-PB  
# Perform the Shapiro-Wilk test for normality
shapiro_test_zoo <- shapiro.test(pcr_raw_zoo_18s$biomass_prop_taxa)

# Print the result of the Shapiro-Wilk test
print(shapiro_test_zoo)

# Plot the histogram
pcr_raw_zoo_18s %>%
  ggplot(aes(x = asin(sqrt(biomass_prop_taxa)))) +   # Set the data and the variable to plot
  geom_histogram(binwidth = .1, color = "black", fill = "#8d9af2", alpha = 0.6) +  # Create the histogram layer
  labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels

# Plot the Q-Q plot
pcr_raw_zoo_18s %>%
  ggplot(aes(sample = asin(sqrt(biomass_prop_taxa)))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot of Random Normal Values", x = "Theoretical Quantiles", y = "Sample Quantiles")



#Add clusters to df for plotting groups
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  select(-Sample_ID_dot) %>% unique()


pcr_raw_zoo_18s=pcr_raw_zoo_18s %>%
  left_join(.,clusters,by="PC1")


pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  ggplot(.,aes(x=((biomass_prop_taxa)), y=((n_reads_pcr))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
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


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/correlations/zooscan_vs_pcr_correlation_18s_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_pcr_correlation_18s_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}

#### Zooscan vs. RRA
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((biomass_prop_taxa)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "Zooscan Biomass Proportion (arcsine square-root)", y = "Raw Relative Abundance (arcsine square-root)", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and Raw Relative Abundance")+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.5)+
  guides(size = FALSE, fill=FALSE) +
  # facet_wrap(~size_fraction, nrow=3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->zoo_vs_raw
zoo_vs_raw


if (saving==1) {
  ggsave(
    filename = here("plots/methods_comparison/correlations/zooscan_vs_raw_correlation_18s_non_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
  
  ggsave(
    filename = here("plots/methods_comparison/zooscan_vs_raw_correlation_18s_non_transformed.pdf"),
    plot = zoo_vs_pcr,
    width = 12,  # Width in inches
    height = 6  # Height in inches
  )
}


### PCR vs Raw
pcr_raw_zoo_18s %>%
  filter(!is.na(cycle.y))%>%
  # filter(cycle.y=="1") %>%
  ggplot(.,aes(x=((n_reads_pcr)), y=((n_reads_raw))))+
  geom_point(aes(shape=cycle.y, size=8,color=as.factor(size_fraction),fill=as.factor(size_fraction)))+
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3"=24, "T1"=23, "T2"=25)) +
  scale_fill_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  scale_color_manual(values=c("#5BA3D5", "#66CC66", "#FF4C38"), labels=c("0.2-0.5 mm","0.5-1 mm","1-2 mm")) +
  # geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +  # Add linear regression line
  labs(x = "PCR Bias-Mitigated Relative Abundance", y = "Raw Reads Relative Abundance", shape = "Cycle", color = "Size Fraction") +  # Add axis labels
  ggtitle("pearson Correlation between Zooscan Biomass Proportion and PCR Bias-Mitigated Relative Abundance")+
  facet_wrap(~size_fraction, nrow=3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))+
  stat_cor(method = "pearson", label.x = 0.8, label.y = 0.8)+
  guides(size = FALSE, fill=FALSE) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))->pcr_vs_raw
pcr_vs_raw

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s_non_transformed.pdf"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)

ggsave(
  filename = here("plots/methods_comparison/correlations/pcr_vs_raw_correlation_18s_non_transformed.png"),
  plot = zoo_vs_raw,
  width = 12,  # Width in inches
  height = 6  # Height in inches
)


