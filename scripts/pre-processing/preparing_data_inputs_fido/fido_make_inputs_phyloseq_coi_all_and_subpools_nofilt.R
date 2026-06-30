# Script to Make inputs from raw OTU data for fido and PCR bias correction
# Barcode: COI
# 


#Load packages
library (tidyverse)
library (here)
library(ggpubr)
library(fido)
library(phyloseq)
library(RColorBrewer)
here()



#Read in the OTU data
#Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("PCR_bias_correction/data/fido/ASV_table_coi_run1.csv")) %>%
  select(-X) 
#Run2, includes pools
asvcoi_run2=read.csv(here("PCR_bias_correction/data/fido/ASV_table_coi_run2.csv")) %>%
  select(-X)



#Phyloseq taxa tables
taxa_coi=read.csv(here("PCR_bias_correction/data/taxa_files/blast_metazoo_coi.csv")) %>% 
  select(-X) %>% 
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>% 
  #Replace Orders that are empty with 'other'
  mutate_all(~replace_na(., "other")) %>% 
  distinct() %>% 
  column_to_rownames("Hash") %>% 
  select(-Subphylum,-Subclass,-Superorder)%>%
  #Deal with weird assignments
  mutate(Genus = ifelse(Genus=="Genus", paste("unidentified ", Family), Genus)) %>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Family = ifelse(Genus == "unidentified Calanoida" & Order == "Calanoida", "unidentified Calanoida", Family))


#Format Long and join as one dataframe
run1_long=asvcoi_run1 %>%
  pivot_longer(cols = 2:ncol(asvcoi_run1), #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )

run2_long=asvcoi_run2%>%
  pivot_longer(cols = 2:ncol(.),  #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )


all_runs=bind_rows(run1_long,run2_long) %>%
  pivot_wider(names_from = Sample_ID, values_from = Nreads)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  column_to_rownames("Hash")

#Replace X in column names
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))

# Calculate sequencing depth per sample
sequencing_depths <- colSums(all_runs)
view(sequencing_depths)

# Compute summary statistics
mean_depth <- mean(sequencing_depths)
sd_depth <- sd(sequencing_depths)
median_depth <- median(sequencing_depths)

# Print summary
cat("Sequencing Depth Summary:\n")
cat(sprintf("Mean ± SD: %.0f ± %.0f reads\n", mean_depth, sd_depth))
cat(sprintf("Median: %.0f reads\n", median_depth))
unique_asvs <- sum(rowSums(all_runs) > 0)
cat(sprintf("Unique ASVs: %d\n", unique_asvs))

#Identify OTUs that aren't in the taxa file

unidentified=all_runs %>% 
  filter(rownames(all_runs) %in% setdiff(rownames(all_runs), rownames(taxa_coi))) %>% 
  mutate(total=rowSums(.))

missing_counts=colSums(unidentified) %>% as.data.frame()


#Separate out by size
#S1
fido_coi_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_coi_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S2"))) %>% 
  filter(rowSums(.) != 0) 
fido_coi_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S3"))) %>% 
  filter(rowSums(.) != 0) 

# Phyloseq filtering: Use phyloseq for filtering and agglomerating
# Phyloseq OTU Tables:
fido_coi_s1_otu=fido_coi_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s2_otu=fido_coi_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s3_otu=fido_coi_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#Phyloseq Taxa table-make an unidentifed Calanoida category. 
# Also handling for fillin in NA and missing data
taxcoi_s1 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s1))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Family = ifelse(Genus == "unidentified Calanoida" & Order == "Calanoida", "unidentified Calanoida", Family))

#Convert to matrix and tax table
taxcoi_s1=  tax_table(as.matrix(taxcoi_s1))

#S2
taxcoi_s2 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s2_otu))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Family = ifelse(Genus == "unidentified Calanoida" & Order == "Calanoida", "unidentified Calanoida", Family))

taxcoi_s2=  tax_table(as.matrix(taxcoi_s2))
#S3
taxcoi_s3 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s3_otu))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Family = ifelse(Genus == "unidentified Calanoida" & Order == "Calanoida", "unidentified Calanoida", Family))

taxcoi_s3=  tax_table(as.matrix(taxcoi_s3))






#Metadata
metacoi=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)



# Agglomerate at the Genus Level -----------------------------------------

# S1 ----------------------------------------------------------------------

fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1, metadata)
fido_coi_s1_genus=tax_glom(fido_coi_s1_phy, taxrank = "Genus")


#Check to make sure we are retaining all OTU counts across agglomeration
# Calculate column sums before tax glomming
colsums_before <- colSums(otu_table(fido_coi_s1_otu))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_coi_s1_genus))

# Find the difference and add the 'other' that was lost to agglomeration
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")

#Check to see which samples had OTUs in lost in agglomeration that need to be added to 'other'
difference==0

#Make inputs for filtering
fido_coi_s1_genus_otu=otu_table(fido_coi_s1_genus) %>% as.data.frame() 

fido_coi_s1_genus_taxa=tax_table(fido_coi_s1_genus) %>% as.data.frame() %>% 
  rownames_to_column("Hash") %>%
  mutate(Hash = ifelse(Genus == "unidentified Calanoida", "unidentified Calanoida", Hash)) %>% 
  distinct() %>% 
  column_to_rownames("Hash")

# Add the new difference row to the dataframe
fido_coi_s1_genus_otu <- bind_rows(fido_coi_s1_genus_otu, difference)

#Check
colSums(fido_coi_s1_otu)==colSums(fido_coi_s1_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s1_genus_taxa, function(x) "other")) %>% 
  column_to_rownames("row_name") %>%
  rbind(.,fido_coi_s1_genus_taxa) -> fido_coi_s1_genus_taxa


#Repeat for other sizes
# S2 ----------------------------------------------------------------------

fido_coi_s2_phy=phyloseq(fido_coi_s2_otu,taxcoi_s2, metadata)
fido_coi_s2_genus=tax_glom(fido_coi_s2_phy, taxrank = "Genus")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums(otu_table(fido_coi_s2_otu))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_coi_s2_genus))

# Find the difference and add the 'other' that was lost to agglomeration
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")
difference==0

#Make inputs for filtering
fido_coi_s2_genus_otu=otu_table(fido_coi_s2_genus) %>% as.data.frame() 

fido_coi_s2_genus_taxa=tax_table(fido_coi_s2_genus) %>% as.data.frame() %>% 
  rownames_to_column("Hash") %>%
  mutate(Hash = ifelse(Genus == "unidentified Calanoida", "unidentified Calanoida", Hash)) %>% 
  distinct() %>% 
  column_to_rownames("Hash")

# Add the new difference row to the dataframe
fido_coi_s2_genus_otu <- bind_rows(fido_coi_s2_genus_otu, difference)

#Check
colSums(fido_coi_s2_otu)==colSums(fido_coi_s2_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s2_genus_taxa, function(x) "other")) %>% 
  column_to_rownames("row_name") %>%
  rbind(.,fido_coi_s2_genus_taxa) -> fido_coi_s2_genus_taxa


# S3 ----------------------------------------------------------------------

fido_coi_s3_phy=phyloseq(fido_coi_s3_otu,taxcoi_s3, metadata)
fido_coi_s3_genus=tax_glom(fido_coi_s3_phy, taxrank = "Genus")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums(otu_table(fido_coi_s3_otu))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_coi_s3_genus))

# Find the difference and add the 'other' that was lost to agglomeration
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")
difference==0

#Make inputs for filtering
fido_coi_s3_genus_otu=otu_table(fido_coi_s3_genus) %>% as.data.frame() 

fido_coi_s3_genus_taxa=tax_table(fido_coi_s3_genus) %>% as.data.frame() %>% 
  rownames_to_column("Hash") %>%
  mutate(Hash = ifelse(Genus == "unidentified Calanoida", "unidentified Calanoida", Hash)) %>% 
  distinct() %>% 
  column_to_rownames("Hash")

# Add the new difference row to the dataframe
fido_coi_s3_genus_otu <- bind_rows(fido_coi_s3_genus_otu, difference)

#Check
colSums(fido_coi_s3_otu)==colSums(fido_coi_s3_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s3_genus_taxa, function(x) "other")) %>% 
  column_to_rownames("row_name") %>%
  rbind(.,fido_coi_s3_genus_taxa) -> fido_coi_s3_genus_taxa




# Merge all sizes  --------------------------------------------------------

#Save aglomerated genus taxa file, replace all columns with 'other' where genus is 'other'
taxcoi_genus_merge=rbind(fido_coi_s1_genus_taxa,fido_coi_s2_genus_taxa,fido_coi_s3_genus_taxa) %>%
  rownames_to_column("Hash") %>% 
  select(-Species) %>% 
  mutate(
    Phylum = if_else(Genus == 'other', 'other', Phylum),
    Class = if_else(Genus == 'other', 'other', Class),

    Order = if_else(Genus == 'other', 'other', Order),
    Family = if_else(Genus == 'other', 'other', Family)
  ) %>% 
  unique() %>% 
  #Filter out unwanted taxa 
  filter(!(Genus == "Strongylocentrotus" & Order =="Camarodonta")) %>%
  filter(!(Genus == "Phascolosoma" & Class =="Polychaeta")) %>%
  filter(!(Genus == "unidentified Calanoida" & Family %in% c("Scolecitrichidae"))) %>%
  filter(!(Genus == "unidentified Calanoida" & Family %in% c("Paracalanidae"))) %>%
  filter(!(Genus == "unidentified Calanoida" & Family %in% c("Euchaetidae"))) %>%
  filter(!(Genus == "unidentified Calanoida" & Family %in% c("other"))) %>%
  filter(!(Genus == "unidentified Calanoida" & Family %in% c(""))) %>% 
  mutate(Hash = ifelse(Genus == "unidentified Calanoida", "unidentified Calanoida", Hash)) 
  
taxcoi_genus=taxcoi_genus_merge %>% select(-Hash) %>% 
  unique()
write.csv(taxcoi_genus,here("PCR_bias_correction/data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv"))



#--------------------------------------------------------------------------------------#


# PART 2: Make OTU and Taxa Inputs for PCR Bias Correction ----------------
# This part will only include the taxa that meet the conditions for the PCR bias correction analysis which include:
# 1. Presence in at least one replicate per treatment per experiment
# 2. The sum of reads in the 'other' category doesn't exceed some threshold, for COI I had to be much looser on this to retain data,
# here I used 55% 

## ==== S1 ====

#Now including subpools, filter so that each final taxa is present in at least one 
#sample in the final df
# Step 1: Define experiments, treatments, and replicates
experiments <- c("Pooled.A1.A3")
treatments <- c("24C", "28C")
replicates <- c("1", "2", "3")

# Generate valid column names
valid_cols <- function(treatment, experiment, replicate) {
  paste0(treatment, "_", experiment, ".", replicate)
}

all_cols <- lapply(experiments, function(exp) {
  sapply(treatments, function(treat) {
    sapply(replicates, function(rep) {
      valid_cols(treat, exp, rep)
    })
  })
}) %>% unlist()

# Retain only columns that exist in the dataset
all_cols <- all_cols[all_cols %in% names(fido_coi_s1_genus_otu)]

# Step 2: Define filtering conditions

# Condition 1: At least one non-zero entry per experiment-treatment combination
filter_condition <- function(df, cols) {
  all(sapply(split(cols, gsub("\\..*$", "", cols)), function(c) {
    any(rowSums(df[c] > 0, na.rm = TRUE) > 0)
  }))
}

# Condition 2: Taxon must be at least x% of the total counts in its experiment-specific columns

filter_percent_experiment <- function(df, all_cols, experiments) {
  percent_set=0.001
  # Create a logical vector to track taxa that pass at least one experiment threshold
  taxa_passes <- rep(FALSE, nrow(df))
  
  for (exp in experiments) {
    # Get columns related to the current experiment
    exp_cols <- all_cols[grep(exp, all_cols)]
    
    if (length(exp_cols) > 0) {
      # Calculate total experiment counts
      total_exp_counts <- colSums(df[exp_cols], na.rm = TRUE)
      min_threshold <- sum(total_exp_counts) * percent_set  # x% of total counts in this experiment
      
      # Identify taxa that meet this threshold in at least one experiment
      taxa_passes <- taxa_passes | rowSums(df[exp_cols], na.rm = TRUE) >= min_threshold
    }
  }
  
  return(taxa_passes)
}

# Step 3: Apply both filters to the dataframe
fido_taxa_filt <- fido_coi_s1_genus_otu %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  filter(
    filter_condition(., all_cols) & 
      filter_percent_experiment(., all_cols, experiments)
  ) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Step 4: Identify "other" taxa (not meeting the criteria)

fido_coi_s1_genus_otu_rename=fido_coi_s1_genus_otu  %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")



other <- fido_coi_s1_genus_otu_rename %>%
  anti_join(fido_taxa_filt) 

# === 1. Join taxonomy to 'other'
other_taxa <- other %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order), by = "Genus")

# === 2. Make "other Calanoida" row
other_calanoida <- other_taxa %>%
  filter(Order == "Calanoida") %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other Calanoida") %>%
  select(Genus, everything())

# === 3. Make "other" row from the rest
other_noncalanoida <- other_taxa %>%
  filter(Order != "Calanoida" | is.na(Order)) %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other") %>%
  select(Genus, everything())

# === 4. Combine both into a final 'other' table
other_final <- bind_rows(other_calanoida, other_noncalanoida)

# === 5. Combine with main filtered taxa
fido_coi_s1_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")


# Recreate `other_final` if needed
# other_final <- bind_rows(other_calanoida, other_noncalanoida)

# Step 5: Visualization of composition within 'other'
# # Join Class info where possible, but label unmatched rows with fallback
# other_bar <- other_taxa %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Class) %>% distinct(), by = "Genus") %>%
#   mutate(Class = case_when(
#     Genus == "other Calanoida" ~ "Calanoida",
#     Genus == "other" ~ "Other",
#     TRUE ~ Class
#   )) %>%
#   pivot_longer(cols = -c(Genus, Class), names_to = "Category", values_to = "Value") %>%
#   group_by(Class, Category) %>%
#   summarise(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum)
# 
# # Stacked bar plot
# ggplot(other_bar, aes(x = Category, y = prop, fill = Class)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of 'Other' by Class",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Class") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Pie chart setup
# n_groups <- length(unique(other_bar$Class))
# base_palette <- brewer.pal(12, "Set3")
# contrasting_palette <- colorRampPalette(base_palette)(n_groups)
# 
# 
# # View for the filtered taxa by Order
# fido_taxa_filt %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Hash, Genus, Order), by = "Genus") %>%
#   select(-Hash, -Genus) %>%
#   pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>%
#   group_by(Category, Order) %>%
#   summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
#   ggplot(aes(x = Category, y = prop, fill = Order)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot by Order (Filtered Taxa)",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Order") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# === Step 6: Identify samples where total 'other' exceeds 80%
# Combine the two other rows before thresholding
other_combined <- other_final %>%
  filter(Genus %in% c("other Calanoida", "other")) %>%
  select(-Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
other_combined[1:5]

# Total reads from original data
sums_fido <- colSums(fido_coi_s1_genus_otu)
sums_fido[1:5]


# Identify samples to drop
thresholds <- sums_fido * 0.8
columns_to_remove <- names(which(colSums(other_combined) > thresholds))

# Prune data
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other_final <- other_final %>% select(c(Genus, setdiff(names(other_final), columns_to_remove)))

# Step 7: Final join
fido_coi_s1_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")

# Convert wide to long format for plotting
fido_coi_s1_final %>%
  rownames_to_column("Genus") %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(total = sum(Reads), prop = Reads / total) %>%
  ungroup() %>%
  ggplot(aes(x = Sample, y = prop, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Genus",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 1))

# Next visualize order
fido_coi_s1_final_order <- fido_coi_s1_final %>%
  rownames_to_column("Genus") %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order) %>% distinct(), by = "Genus") %>%
  mutate(Order = case_when(
    Genus == "other Calanoida" ~ "Calanoida",
    Genus == "other" ~ "Other",
    is.na(Order) ~ "Unclassified",
    TRUE ~ Order
  )) %>%  
  distinct(.)

# Step 2: Aggregate by Order and Sample
fido_order_long <- fido_coi_s1_final_order %>%
  select(-Genus) %>%
  pivot_longer(cols = -Order, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample, Order) %>%
  summarise(Order_sum = sum(Reads), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(Order_sum),
         prop = Order_sum / sample_total)

# Step 3: Plot
ggplot(fido_order_long, aes(x = Sample, y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Order",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))


# Step 8: Save
write.csv(fido_coi_s1_save_genus_phy <- fido_coi_s1_final,
          here("PCR_bias_correction/data/fido/phy/fido_coi_s1_genus_phy_all_subpools_nofilt.csv"))

#Refresh before next size
rm(other)


## ==== S2 ====
#Now including subpools, filter so that each final taxa is present in at least one 
#sample in the final df
# Step 1: Define experiments, treatments, and replicates
experiments <- c("Pooled.A1.A3")
treatments <- c("24C", "28C")
replicates <- c("1", "2", "3")

# Generate valid column names
valid_cols <- function(treatment, experiment, replicate) {
  paste0(treatment, "_", experiment, ".", replicate)
}

all_cols <- lapply(experiments, function(exp) {
  sapply(treatments, function(treat) {
    sapply(replicates, function(rep) {
      valid_cols(treat, exp, rep)
    })
  })
}) %>% unlist()

# Retain only columns that exist in the dataset
all_cols <- all_cols[all_cols %in% names(fido_coi_s2_genus_otu)]

# Step 2: Define filtering conditions

# Condition 1: At least one non-zero entry per experiment-treatment combination
filter_condition <- function(df, cols) {
  all(sapply(split(cols, gsub("\\..*$", "", cols)), function(c) {
    any(rowSums(df[c] > 0, na.rm = TRUE) > 0)
  }))
}

# Condition 2: Taxon must be at least x% of the total counts in its experiment-specific columns

filter_percent_experiment <- function(df, all_cols, experiments) {
  percent_set=0.001
  # Create a logical vector to track taxa that pass at least one experiment threshold
  taxa_passes <- rep(FALSE, nrow(df))
  
  for (exp in experiments) {
    # Get columns related to the current experiment
    exp_cols <- all_cols[grep(exp, all_cols)]
    
    if (length(exp_cols) > 0) {
      # Calculate total experiment counts
      total_exp_counts <- colSums(df[exp_cols], na.rm = TRUE)
      min_threshold <- sum(total_exp_counts) * percent_set  # x% of total counts in this experiment
      
      # Identify taxa that meet this threshold in at least one experiment
      taxa_passes <- taxa_passes | rowSums(df[exp_cols], na.rm = TRUE) >= min_threshold
    }
  }
  
  return(taxa_passes)
}

# Step 3: Apply both filters to the dataframe
fido_taxa_filt <- fido_coi_s2_genus_otu %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  filter(
    filter_condition(., all_cols) & 
      filter_percent_experiment(., all_cols, experiments)
  ) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Step 4: Identify "other" taxa (not meeting the criteria)

fido_coi_s2_genus_otu_rename=fido_coi_s2_genus_otu  %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")



other <- fido_coi_s2_genus_otu_rename %>%
  anti_join(fido_taxa_filt) 

# === 1. Join taxonomy to 'other'
other_taxa <- other %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order), by = "Genus")

# === 2. Make "other Calanoida" row
other_calanoida <- other_taxa %>%
  filter(Order == "Calanoida") %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other Calanoida") %>%
  select(Genus, everything())

# === 3. Make "other" row from the rest
other_noncalanoida <- other_taxa %>%
  filter(Order != "Calanoida" | is.na(Order)) %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other") %>%
  select(Genus, everything())

# === 4. Combine both into a final 'other' table
other_final <- bind_rows(other_calanoida, other_noncalanoida)

# === 5. Combine with main filtered taxa
fido_coi_s2_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")


# Recreate `other_final` if needed
# other_final <- bind_rows(other_calanoida, other_noncalanoida)

# Step 5: Visualization of composition within 'other'
# # Join Class info where possible, but label unmatched rows with fallback
# other_bar <- other_taxa %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Class) %>% distinct(), by = "Genus") %>%
#   mutate(Class = case_when(
#     Genus == "other Calanoida" ~ "Calanoida",
#     Genus == "other" ~ "Other",
#     TRUE ~ Class
#   )) %>%
#   pivot_longer(cols = -c(Genus, Class), names_to = "Category", values_to = "Value") %>%
#   group_by(Class, Category) %>%
#   summarise(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum)
# 
# # Stacked bar plot
# ggplot(other_bar, aes(x = Category, y = prop, fill = Class)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of 'Other' by Class",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Class") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Pie chart setup
# n_groups <- length(unique(other_bar$Class))
# base_palette <- brewer.pal(12, "Set3")
# contrasting_palette <- colorRampPalette(base_palette)(n_groups)
# 
# 
# # View for the filtered taxa by Order
# fido_taxa_filt %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Hash, Genus, Order), by = "Genus") %>%
#   select(-Hash, -Genus) %>%
#   pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>%
#   group_by(Category, Order) %>%
#   summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
#   ggplot(aes(x = Category, y = prop, fill = Order)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot by Order (Filtered Taxa)",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Order") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# === Step 6: Identify samples where total 'other' exceeds 80%
# Combine the two other rows before thresholding
other_combined <- other_final %>%
  filter(Genus %in% c("other Calanoida", "other")) %>%
  select(-Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
other_combined[1:5]

# Total reads from original data
sums_fido <- colSums(fido_coi_s2_genus_otu)
sums_fido[1:5]


# Identify samples to drop
thresholds <- sums_fido * 0.8
columns_to_remove <- names(which(colSums(other_combined) > thresholds))

# Prune data
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other_final <- other_final %>% select(c(Genus, setdiff(names(other_final), columns_to_remove)))

# Step 7: Final join
fido_coi_s2_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")

# Convert wide to long format for plotting
fido_coi_s2_final %>%
  rownames_to_column("Genus") %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(total = sum(Reads), prop = Reads / total) %>%
  ungroup() %>%
  ggplot(aes(x = Sample, y = prop, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Genus",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 1))


# Step 8: Save
write.csv(fido_coi_s2_save_genus_phy <- fido_coi_s2_final,
          here("PCR_bias_correction/data/fido/phy/fido_coi_s2_genus_phy_all_subpools_nofilt.csv"))

# Visualize Order
fido_coi_s2_final_order <- fido_coi_s2_final %>%
  rownames_to_column("Genus") %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order) %>% distinct(), by = "Genus") %>%
  mutate(Order = case_when(
    Genus == "other Calanoida" ~ "Calanoida",
    Genus == "other" ~ "Other",
    is.na(Order) ~ "Unclassified",
    TRUE ~ Order
  ))%>%  
  distinct(.)

# Step 2: Aggregate by Order and Sample
fido_order_long <- fido_coi_s2_final_order %>%
  select(-Genus) %>%
  pivot_longer(cols = -Order, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample, Order) %>%
  summarise(Order_sum = sum(Reads), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(Order_sum),
         prop = Order_sum / sample_total)

# Step 3: Plot
ggplot(fido_order_long, aes(x = Sample, y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Order",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))

#Refresh before next size
rm(other)


## ==== S3 ====
#Now including subpools, filter so that each final taxa is present in at least one 
#sample in the final df
# Step 1: Define experiments, treatments, and replicates
experiments <- c("Pooled.A1.A3")
treatments <- c("24C", "28C")
replicates <- c("1", "2", "3")

# Generate valid column names
valid_cols <- function(treatment, experiment, replicate) {
  paste0(treatment, "_", experiment, ".", replicate)
}

all_cols <- lapply(experiments, function(exp) {
  sapply(treatments, function(treat) {
    sapply(replicates, function(rep) {
      valid_cols(treat, exp, rep)
    })
  })
}) %>% unlist()

# Retain only columns that exist in the dataset
all_cols <- all_cols[all_cols %in% names(fido_coi_s3_genus_otu)]

# Step 2: Define filtering conditions

# Condition 1: At least one non-zero entry per experiment-treatment combination
filter_condition <- function(df, cols) {
  all(sapply(split(cols, gsub("\\..*$", "", cols)), function(c) {
    any(rowSums(df[c] > 0, na.rm = TRUE) > 0)
  }))
}

# Condition 2: Taxon must be at least x% of the total counts in its experiment-specific columns

filter_percent_experiment <- function(df, all_cols, experiments) {
  percent_set=0.001
  # Create a logical vector to track taxa that pass at least one experiment threshold
  taxa_passes <- rep(FALSE, nrow(df))
  
  for (exp in experiments) {
    # Get columns related to the current experiment
    exp_cols <- all_cols[grep(exp, all_cols)]
    
    if (length(exp_cols) > 0) {
      # Calculate total experiment counts
      total_exp_counts <- colSums(df[exp_cols], na.rm = TRUE)
      min_threshold <- sum(total_exp_counts) * percent_set  # x% of total counts in this experiment
      
      # Identify taxa that meet this threshold in at least one experiment
      taxa_passes <- taxa_passes | rowSums(df[exp_cols], na.rm = TRUE) >= min_threshold
    }
  }
  
  return(taxa_passes)
}

# Step 3: Apply both filters to the dataframe
fido_taxa_filt <- fido_coi_s3_genus_otu %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  filter(
    filter_condition(., all_cols) & 
      filter_percent_experiment(., all_cols, experiments)
  ) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Step 4: Identify "other" taxa (not meeting the criteria)

fido_coi_s3_genus_otu_rename=fido_coi_s3_genus_otu  %>%
  rownames_to_column("Hash") %>%
  left_join(
    taxa_coi %>% 
      rownames_to_column("Hash") %>% 
      select(Hash, Genus), 
    by = "Hash"
  ) %>%
  select(-Hash) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")



other <- fido_coi_s3_genus_otu_rename %>%
  anti_join(fido_taxa_filt) 

# === 1. Join taxonomy to 'other'
other_taxa <- other %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order), by = "Genus")

# === 2. Make "other Calanoida" row
other_calanoida <- other_taxa %>%
  filter(Order == "Calanoida") %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other Calanoida") %>%
  select(Genus, everything())

# === 3. Make "other" row from the rest
other_noncalanoida <- other_taxa %>%
  filter(Order != "Calanoida" | is.na(Order)) %>%
  select(-Genus, -Order) %>%
  distinct(.) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  mutate(Genus = "other") %>%
  select(Genus, everything())

# === 4. Combine both into a final 'other' table
other_final <- bind_rows(other_calanoida, other_noncalanoida)

# === 5. Combine with main filtered taxa
fido_coi_s3_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")


# Recreate `other_final` if needed
# other_final <- bind_rows(other_calanoida, other_noncalanoida)

# Step 5: Visualization of composition within 'other'
# # Join Class info where possible, but label unmatched rows with fallback
# other_bar <- other_taxa %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Class) %>% distinct(), by = "Genus") %>%
#   mutate(Class = case_when(
#     Genus == "other Calanoida" ~ "Calanoida",
#     Genus == "other" ~ "Other",
#     TRUE ~ Class
#   )) %>%
#   pivot_longer(cols = -c(Genus, Class), names_to = "Category", values_to = "Value") %>%
#   group_by(Class, Category) %>%
#   summarise(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum)
# 
# # Stacked bar plot
# ggplot(other_bar, aes(x = Category, y = prop, fill = Class)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of 'Other' by Class",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Class") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Pie chart setup
# n_groups <- length(unique(other_bar$Class))
# base_palette <- brewer.pal(12, "Set3")
# contrasting_palette <- colorRampPalette(base_palette)(n_groups)
# 
# 
# # View for the filtered taxa by Order
# fido_taxa_filt %>%
#   left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Hash, Genus, Order), by = "Genus") %>%
#   select(-Hash, -Genus) %>%
#   pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>%
#   group_by(Category, Order) %>%
#   summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
#   ungroup() %>%
#   group_by(Category) %>%
#   mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
#   ggplot(aes(x = Category, y = prop, fill = Order)) +
#   geom_bar(stat = "identity") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot by Order (Filtered Taxa)",
#        x = "Sample",
#        y = "Proportion",
#        fill = "Order") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# === Step 6: Identify samples where total 'other' exceeds 80%
# Combine the two other rows before thresholding
other_combined <- other_final %>%
  filter(Genus %in% c("other Calanoida", "other")) %>%
  select(-Genus) %>%
  summarise(across(everything(), sum, na.rm = TRUE))
other_combined[1:5]

# Total reads from original data
sums_fido <- colSums(fido_coi_s3_genus_otu)
sums_fido[1:5]


# Identify samples to drop
thresholds <- sums_fido * 0.8
columns_to_remove <- names(which(colSums(other_combined) > thresholds))

# Prune data
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other_final <- other_final %>% select(c(Genus, setdiff(names(other_final), columns_to_remove)))

# Step 7: Final join
fido_coi_s3_final <- bind_rows(fido_taxa_filt, other_final) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Genus")

# Convert wide to long format for plotting
fido_coi_s3_final %>%
  rownames_to_column("Genus") %>%
  pivot_longer(cols = -Genus, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample) %>%
  mutate(total = sum(Reads), prop = Reads / total) %>%
  ungroup() %>%
  ggplot(aes(x = Sample, y = prop, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Genus",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 1))



# Next visualize order
fido_coi_s3_final_order <- fido_coi_s3_final %>%
  rownames_to_column("Genus") %>%
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% select(Genus, Order) %>% distinct(), by = "Genus") %>%
  mutate(Order = case_when(
    Genus == "other Calanoida" ~ "Calanoida",
    Genus == "other" ~ "Other",
    is.na(Order) ~ "Unclassified",
    TRUE ~ Order
  )) %>%  
  distinct(.)

# Step 2: Aggregate by Order and Sample
fido_order_long <- fido_coi_s3_final_order %>%
  select(-Genus) %>%
  pivot_longer(cols = -Order, names_to = "Sample", values_to = "Reads") %>%
  group_by(Sample, Order) %>%
  summarise(Order_sum = sum(Reads), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(Order_sum),
         prop = Order_sum / sample_total)

# Step 3: Plot
ggplot(fido_order_long, aes(x = Sample, y = prop, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Final Taxa Composition by Order",
       x = "Sample",
       y = "Proportion of Total Reads",
       fill = "Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))


# Step 8: Save
write.csv(fido_coi_s3_save_genus_phy <- fido_coi_s3_final,
          here("PCR_bias_correction/data/fido/phy/fido_coi_s3_genus_phy_all_subpools_nofilt.csv"))

#Refresh before next size
rm(other)
