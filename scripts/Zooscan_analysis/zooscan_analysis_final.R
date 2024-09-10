#Zooscan analysis final
librarian::shelf(tidyverse, googledrive, stringr,here,vegan,ggpubr)
here()
source(here("scripts/Zooscan_analysis/zooscan_functions.R"))
source("scripts/helpful_functions/general_helper_functions.R")


#Load in the data
##Lat and lon are funky so add back in from metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.2.2023_for_map.csv"))%>%
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
  # Standardize key columns (some samples are incorerctly notated to compare with Zooscan)
  mutate(sample_id = tolower(str_replace_all(Sample_ID_short, "-", "_"))) %>%
  select(-X,-Sample_ID_dot, -Sizefractionmm, -max_size) %>%
  distinct(.) %>%
  #Make PC1 values opposite for plotting
  mutate(PC1=PC1*-1) 


# List all .tsv files in the folder
tsv_files <- list.files(here("data/Zooscan/"), pattern = "\\.tsv$", full.names = TRUE)

# Find the most recently added file
latest_ecotaxa <- tsv_files[which.max(file.info(tsv_files)$ctime)]
latest_ecotaxa

# Read the .tsv file into a data frame
zooscan_exp <- read.table(latest_ecotaxa, header=TRUE, sep="\t", encoding="latin1")


#Basic plots to show the relative abundance and biomass of zooscan results from each station
zooscan_processed=readEcotaxa(zooscan_exp)%>%
  #Cleaning up fomratting issues
  mutate(sample_id = str_replace_all(sample_id, "-", "_"),
         sample_id = ifelse(sample_id == "c2_t1_h36", "ct2_t1_h36", sample_id),
         sample_id = ifelse(sample_id == "ct1_t8_h10", "c1_t8_h10", sample_id),
         sample_id = ifelse(sample_id == "ct2_t9_h19", "c2_t9_h19", sample_id),
         sample_id = ifelse(sample_id == "c3_bt6_h25", "c3_t6_h25", sample_id)) %>%
  #Need to fix hyperiids
  filter(!(object_annotation_category %in% c("Hyperiidea","part<Crustacea","darksphere", "multiple organisms","head<Chaetognatha",
                                             "egg<Actinopterygii"))) %>%
  #Fix missing volume filtered 
  mutate(sample_tot_vol = case_when(
    sample_id == "ct1_t1_h28" ~ 279,
    sample_id == "ct1_t2_h29" ~ 278,
    TRUE ~ sample_tot_vol  # This line keeps the original values for all other rows
  ))





# Analysis ----------------------------------------------------------------


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
  mutate(size_fraction=factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2", ">2"))) %>%
  left_join(metadata, by="sample_id") %>% 
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
  facet_wrap(~size_fraction, nrow=4)+
  labs(x = "Offfshore \u2190 PC1 \u2192 Onshore", y = expression(italic("H'")), title="COI") +
  scale_color_discrete(name = "Cycle") +  # Adjust color legend label
  stat_cor(method="pearson", label.x = 4, label.y = 3.5)+
  theme_classic()



#Calculate C-biomass
 zooscan_processed %>%
   transform_by_taxa_group(.,"esd") %>%
   mutate(dryweight_C_mg=dryweight_C_ug/1000) %>%
   #Add log biomass, and biomass/m2
   mutate(log10_dryweight_C_mg=log10(dryweight_C_mg))->zooscan_biomass 

#Biomass histogram
 zooscan_biomass %>% 
   filter(object_annotation_category %in% c("Calanoida","Copepoda<Maxillopoda","Oithonidae","Harpacticoida", "Poecilostomatoida"))%>%
   filter(object_annotation_category=="Calanoida") %>%
 ggplot(., aes(x = log10((dryweight_C_ug)))) +   # Set the data and the variable to plot
   geom_histogram(binwidth = 0.2, color = "black", fill = "lightblue", alpha = 0.6) +  # Create the histogram layer
   labs(title = "Histogram of Random Normal Values", x = "Values", y = "Frequency")  # Add titles and labels


#Compute relative abudances 
relative_abundances=zooscan_processed %>%
  group_by(sample_id,size_fraction,object_annotation_category) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         relative_abundance = count / total)





# Merge dataframes: add metadata to relative abudance data
relative_abundances_map <- relative_abundances %>%
  left_join(metadata, by="sample_id") %>%
  distinct(.)

biomass_map = zooscan_biomass %>%
  left_join(.,metadata, by="sample_id") %>% 
  #Make size class a factor
  mutate(size_fraction=factor(size_fraction, levels = c("0.2-0.5", "0.5-1", "1-2", ">2"))) %>%
    distinct(.) 

#Save main biomass df
write.csv(biomass_map,here("data/Zooscan/zooscan_biomass_all.csv"))


#Filter to calalnoids, merge by sample and size and then compute biomass/m2 and log biomass/m2
#General Zooscan dataframe for plotting
zooscan_by_sample = biomass_map %>%
  group_by(size_fraction, sample_id) %>%
  mutate(dryweight_C_mg_sum_sample = sum(dryweight_C_mg, na.rm = TRUE),
            dryweight_C_ug_sum_sample = sum(dryweight_C_ug, na.rm = TRUE)) %>%
  group_by(size_fraction, sample_id,object_annotation_category) %>%
  summarise(dryweight_C_mg_sum_taxa = sum(dryweight_C_mg, na.rm = TRUE),
            dryweight_C_ug_sum_taxa = sum(dryweight_C_ug, na.rm = TRUE),
            dryweight_C_mg_sum_sample=mean(dryweight_C_mg_sum_sample),
            dryweight_C_ug_sum_sample=mean(dryweight_C_ug_sum_sample),
            sample_conc=mean(sample_conc),
            sample_tot_vol=mean(sample_tot_vol),
            acq_sub_part=mean(acq_sub_part))%>%
  left_join(metadata, by="sample_id") %>%
  mutate(log10_dryweight_C_ug_m2_taxa=log(dryweight_C_ug_sum_taxa),
         dryweight_C_mg_m2_taxa=dryweight_C_mg_sum_taxa*sample_conc,
         log10_dryweight_C_ug_m2_sample=log(dryweight_C_ug_sum_sample),
         dryweight_C_mg_m2_sample=dryweight_C_mg_sum_sample*sample_conc) %>%
  mutate(biomass_prop_taxa=dryweight_C_mg_m2_taxa/dryweight_C_mg_m2_sample,
         biomass_taxa_clr=clr_convert(dryweight_C_mg_m2_taxa)) %>% 
  distinct(biomass_taxa_clr, .keep_all = TRUE) 



# Calanoids ---------------------------------------------------------------

#CALANOID DATA FRAME FOR PLOTTING
zoop_calanoid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Calanoida")

#Calanoid dataframe for relative abundances
zoop_calanoid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Calanoida")

#Save Biomass and relative abundances
write.csv(zoop_calanoid_by_sample,here("data/Zooscan/zoop_calanoid_by_sample_biomass.csv"))
write.csv(relative_abundances_map,here("data/Zooscan/zoop_relative_abundance.csv"))
write.csv(zooscan_by_sample,here("data/Zooscan/zooscan_by_sample_biomass.csv"))





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
write.csv(zoop_euphausiid_by_sample,here("data/Zooscan/zoop_euphausiid_by_sample_biomass.csv"))
write.csv(zoop_euphausiid_by_sample_relative_abundance,here("data/Zooscan/zoop_euphausiid_by_sample_relative_abundance.csv"))




# Oithonidae --------------------------------------------------------------

#Euphausiid DATA FRAME FOR PLOTTING
zoop_euphausiid_by_sample = zooscan_by_sample %>%
  filter(object_annotation_category=="Euphausiacea")

#euphausiid dataframe for relative abundances
zoop_euphausiid_by_sample_relative_abundance = relative_abundances_map %>%
  filter(object_annotation_category=="Euphausiacea")

#Save Biomass and relative abundances
write.csv(zoop_euphausiid_by_sample,here("data/Zooscan/zoop_euphausiid_by_sample_biomass.csv"))
write.csv(zoop_euphausiid_by_sample_relative_abundance,here("data/Zooscan/zoop_euphausiid_by_sample_relative_abundance.csv"))





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
