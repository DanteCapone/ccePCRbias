#Make inputs that inlcude all sub-pools and All pool for mutliple observartions of amp_eff
library (tidyverse)
library (here)
library(ggpubr)
library(fido)
library(phyloseq)
here()



###18S
#Read in the OTU data
#Run 1 (Non pooled data)
asv18s_run1=read.csv(here("PCR_bias_correction/data/raw_reads/","ASV_table_18s_run1.csv")) %>%
  select(-X) 
#Run2
asv18s_run2=read.csv(here("PCR_bias_correction/data/raw_reads/","ASV_table_18s_run2.csv")) %>%
  select(-X)



#Taxa Tables using new combined taxa file between BLAST and Metazoogene
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
  column_to_rownames("Hash")


#Format Long
run1_long=asv18s_run1 %>%
  pivot_longer(cols = 2:ncol(asv18s_run1), #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )

run2_long=asv18s_run2%>%
  pivot_longer(cols = 2:ncol(.),  #Specify the columns to pivot
               names_to = "Sample_ID", #Name of the new variable column
               values_to = "Nreads" #Name of the new value column
  )



all_runs=bind_rows(run1_long,run2_long) %>%
  pivot_wider(names_from = Sample_ID, values_from = Nreads)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  column_to_rownames("Hash")

#Replace X
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))


#Identify OTUs that aren't in the taxa file

unidentified=all_runs %>% 
  filter(rownames(all_runs) %in% setdiff(rownames(all_runs), rownames(taxa_18s))) %>% 
  mutate(total=rowSums(.))

missing_counts=colSums(unidentified) %>% as.data.frame()



#Separate out by size
#S1
fido_18s_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_18s_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S2"))) %>% 
  filter(rowSums(.) != 0)
fido_18s_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("A1"),
                  contains("B3"),contains("C5"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)


###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#Make taxa table and add a category for unidentified calanoida
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Calanoida", "unidentified Calanoida", Family))%>% 
  mutate(Family = ifelse(Family == "other" & Order == "Collodaria", "unidentified Collodaria", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Siphonophorae", "unidentified Siphonophorae", Family)) 
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Collodaria", "unidentified Collodaria", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Siphonophorae", "unidentified Siphonophorae", Family)) 
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Calanoida", "unidentified Calanoida", Family))%>% 
  mutate(Family = ifelse(Family == "other" & Order == "Collodaria", "unidentified Collodaria", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Siphonophorae", "unidentified Siphonophorae", Family)) 
tax18s_s3=  tax_table(as.matrix(tax18s_s3))






#Metadata
meta18s=read.csv(here("PCR_bias_correction/data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)



# Agglomerate at the Family Level -----------------------------------------

# S1 ----------------------------------------------------------------------

fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1, metadata)
fido_18s_s1_family=tax_glom(fido_18s_s1_phy, taxrank = "Family")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums(otu_table(fido_18s_s1_otu))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s1_family))

# Find the difference and add the 'other' that was lost to agglomeration
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")
difference==0

#Make inputs for filtering
fido_18s_s1_family_otu=otu_table(fido_18s_s1_family) %>% as.data.frame() 

fido_18s_s1_family_taxa=tax_table(fido_18s_s1_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s1_family_otu <- bind_rows(fido_18s_s1_family_otu, difference)

#Check
colSums(fido_18s_s1_otu)==colSums(fido_18s_s1_family_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s1_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s1_family_taxa) -> fido_18s_s1_family_taxa



# S2 ----------------------------------------------------------------------
fido_18s_s2_phy=phyloseq(fido_18s_s2_otu,tax18s_s2, metadata)
fido_18s_s2_family=tax_glom(fido_18s_s2_phy, taxrank = "Family")




#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s2))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s2_family))

# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s2_family_otu=otu_table(fido_18s_s2_family) %>% as.data.frame()
fido_18s_s2_family_taxa=tax_table(fido_18s_s2_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s2_family_otu <- bind_rows(fido_18s_s2_family_otu, difference)


colSums(fido_18s_s2_otu)[1:5]
colSums(fido_18s_s2_family_otu)[1:5]


#Make inputs for filtering
fido_18s_s2_family_otu=otu_table(fido_18s_s2_family) %>% as.data.frame()
fido_18s_s2_family_taxa=tax_table(fido_18s_s2_family) %>% as.data.frame()

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s2_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s2_family_taxa) -> fido_18s_s2_family_taxa



# S3 ----------------------------------------------------------------------
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3)
fido_18s_s3_family=tax_glom(fido_18s_s3_phy, taxrank = "Family")


#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s3))

# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s3_family))

# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s3_family_otu=otu_table(fido_18s_s3_family) %>% as.data.frame() 

fido_18s_s3_family_taxa=tax_table(fido_18s_s3_family) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s3_family_otu <- bind_rows(fido_18s_s3_family_otu, difference)


colSums(fido_18s_s3_otu)[1:5]
colSums(fido_18s_s3_family_otu)[1:5]



#Make inputs for filtering
fido_18s_s3_family_otu=otu_table(fido_18s_s3_family) %>% as.data.frame()
fido_18s_s3_family_taxa=tax_table(fido_18s_s3_family) %>% as.data.frame()

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s3_family_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s3_family_taxa) -> fido_18s_s3_family_taxa

#Save aglomerated family taxa file, replace all columns with 'other' where family is 'other'
tax18s_family=rbind(fido_18s_s1_family_taxa,fido_18s_s2_family_taxa,fido_18s_s3_family_taxa) %>%
  rownames_to_column("Hash") %>% 
  select(-Subphylum,-Subclass,-Superorder,-Species, -Genus,-Hash) %>% 
  mutate(
    Phylum = if_else(Family == 'other', 'other', Phylum),
    # Subphylum = if_else(Family == 'other', 'other', Subphylum),
    Class = if_else(Family == 'other', 'other', Class),
    # Subclass = if_else(Family == 'other', 'other', Subclass),
    # Superorder = if_else(Family == 'other', 'other', Superorder),
    Order = if_else(Family == 'other', 'other', Order),
  ) %>% 
  unique() %>% 
  #Additionally filter out reclassified Corycaeidae (actually a poecilistomatoida)
  filter(!(Family == "Corycaeidae" & Order =="Cyclopoida")) 
  
write.csv(tax18s_family,here("PCR_bias_correction/data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))




# PART 2: MAKE OTU Tables -------------------------------------------------



## ==== S1 ====

#Now including subpools, filter so that each final taxa is present in at least one 
#sample in the final df
# Step 1: Define experiments, treatments, and replicates
experiments <- c("AllPool", "Pooled.A1.A3", "Pooled.B3.B5", "Pooled.C5.C7")
treatments <- c("20C", "24C", "28C")
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
all_cols <- all_cols[all_cols %in% names(fido_18s_s1_family_otu)]

# Step 2: Define filtering conditions

# Condition 1: At least one non-zero entry per experiment-treatment combination
filter_condition <- function(df, cols) {
  all(sapply(split(cols, gsub("\\..*$", "", cols)), function(c) {
    any(rowSums(df[c] > 0, na.rm = TRUE) > 0)
  }))
}

# Condition 2: Taxon must be at least x% of the total counts in its experiment-specific columns

filter_percent_experiment <- function(df, all_cols, experiments) {
  percent_set=0.01
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
fido_taxa_filt <- fido_18s_s1_family_otu %>%
  rownames_to_column("Hash") %>%
  filter(filter_condition(., all_cols) & filter_percent_experiment(., all_cols, experiments)) 

# Step 4: Identify "other" taxa (not meeting the criteria)
other <- fido_18s_s1_family_otu %>%
  anti_join(fido_taxa_filt)

# Step 5: Visualization of taxa composition before proceeding
other %>%
  rownames_to_column("Hash") %>%
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% select(Hash, Order), by = "Hash") %>%
  select(-Hash) %>%
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>%
  group_by(Order, Category) %>%
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
  ungroup() %>%
  group_by(Category) %>%
  mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
  ggplot(aes(x = Category, y = taxa_sum, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# View for the filtered taxa
fido_taxa_filt %>%
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% select(Hash, Order), by = "Hash") %>%
  select(-Hash) %>%
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>%
  group_by(Order, Category) %>%
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
  ungroup() %>%
  group_by(Category) %>%
  mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
  ggplot(aes(x = Category, y = taxa_sum, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot if needed
# ggsave(here("plots/pre_processing/QC/other_order_composition_s1.png"), width = 10, height = 6, units = "in")
# ggsave(here("plots/pre_processing/QC/other_order_composition_s1.pdf"), width = 18, height = 6, units = "in")

# Remove "other" from `fido_18s_s1_family_otu`
other %>%
  summarise_all(sum) %>%
  mutate(Hash = "other") -> other

# Step 6: Identify columns where 'other' exceeds 20% of total sample counts
sums_other <- colSums(other %>% select(-Hash))
sums_fido <- colSums(fido_18s_s1_family_otu)

thresholds <- sums_fido * 0.2
columns_to_remove <- names(which(sums_other > thresholds))

# Remove these columns
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other <- other %>% select(-all_of(columns_to_remove))

# Step 7: Combine data
fido_18s_s1_final <- rbind(fido_taxa_filt, other) %>%
  group_by(Hash) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Hash")

colSums(fido_18s_s1_family_otu)[1:5]
colSums(fido_18s_s1_final)[1:5]

# Step 8: Join with taxonomy file
fido_18s_s1_save_family_phy <- fido_18s_s1_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_18s_s1_family_taxa %>% rownames_to_column("Hash"), by = "Hash") %>%
  mutate(Family = if_else(
    Family == "other" & Order != "other",
    paste0("unidentified ", Order),
    Family
  )) %>%
  mutate(Family = ifelse(Hash == "other", "other", Family)) %>%
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

  
  
  #Save
  write.csv(fido_18s_s1_save_family_phy,here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv"))

# Remove "unidentified Collodaria"

fido_18s_s1_save_family_phy_nocollo= fido_18s_s1_save_family_phy %>%
  filter(Family != "unidentified Collodaria") %>% 
  pivot_longer(cols = -Family, names_to = "sample", values_to = "read_count") %>%
  # filter(!str_detect(sample, "All")) %>%
  mutate(global_median = median(read_count, na.rm = TRUE)) %>%
  group_by(sample) %>%
  mutate(sample_median = median(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(norm_read_count = round(read_count * (global_median / sample_median))) %>%
  select(Family, sample, norm_read_count) %>%
  pivot_wider(names_from = sample, values_from = norm_read_count)

# Save cleaned dataset
write.csv(fido_18s_s1_save_family_phy_nocollo, 
          here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools_nocollodaria.csv"), 
          row.names = FALSE)


## ==== s2 ====
# Filter out any non-existent columns from your dataframe
all_cols <- all_cols[all_cols %in% names(fido_18s_s2_family_otu)]


# Apply the filter across the dataframe
fido_taxa_filt <- fido_18s_s2_family_otu%>%
  rownames_to_column("Hash") %>%
  filter(filter_condition(., all_cols) & filter_percent_experiment(., all_cols, experiments)) 

other <- fido_18s_s2_family_otu %>%
  anti_join(fido_taxa_filt) 


#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  select(-Hash) %>% 
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>% 
  group_by(Order, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  # filter(Order=="Calanoida") %>%
  ggplot(., aes(x = Category, y = taxa_sum, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


other %>% 
  summarise_all(sum) %>% 
  mutate(Hash = "other")->other

#Find columns where other exceeds threshold of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% select(-Hash))
sums_fido <- colSums(fido_18s_s2_family_otu)

# Step 2: Calculate threshold of the column sums of fido_18s_s3_family_otu
thresholds <- sums_fido * 0.2

# Step 3: Identify columns where the sum of `other` is greater than threshold of the sum of `fido_18s_s2_family_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_18s_s3_family_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))


# Combine data
fido_18s_s2_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Join with taxa file
fido_18s_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s2_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Family = if_else(
    Family == "other" & Order != "other",
    paste0("unidentified ", Order),
    Family
  )) %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s2_save_family_phy

#Save
write.csv(fido_18s_s2_save_family_phy,here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv"))

fido_18s_s2_save_family_phy%>% 
  pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>% 
  group_by(Family, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  ggplot(., aes(x = Category, y = taxa_sum, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Family",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Remove "unidentified Collodaria"

fido_18s_s2_save_family_phy_nocollo <- fido_18s_s2_save_family_phy %>%
  filter(Family != "unidentified Collodaria") 

# Save cleaned dataset
write.csv(fido_18s_s2_save_family_phy_nocollo, 
          here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools_nocollodaria.csv"), 
          row.names = FALSE)



## ==== s3 ====

# Filter out any non-existent columns from your dataframe
all_cols <- all_cols[all_cols %in% names(fido_18s_s3_family_otu)]

# Apply the filter across the dataframe
fido_taxa_filt <- fido_18s_s3_family_otu%>%
  rownames_to_column("Hash") %>%
  filter(filter_condition(., all_cols) & filter_percent_experiment(., all_cols, experiments)) 


other <- fido_18s_s3_family_otu %>%
  anti_join(fido_taxa_filt) 

fido_18s_s3_family_otu %>% 
  summarise_all(sum) %>% 
  mutate(Order = "total")->total

#Check composition before proceeding

other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  select(-Hash) %>% 
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>% 
  group_by(Order, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  # filter(Order=="Calanoida") %>%
  ggplot(., aes(x = Category, y = taxa_sum, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Save plot 
# ggsave(here("plots/pre_processing/QC/other_order_composition_s3.png"), width = 10, height = 6, units = "in")
# ggsave(here("plots/pre_processing/QC/other_order_composition_s3.pdf"), width = 18, height = 6, units = "in")

#Find columns where other exceeds threshold of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other)
sums_fido <- colSums(fido_18s_s3_family_otu)

# Step 2: Calculate threshold of the column sums of fido_18s_s3_family_otu
thresholds <- sums_fido * 0.20

# Step 3: Identify columns where the sum of `other` is greater than threshold of the sum of `fido_18s_s3_family_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_18s_s3_family_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))

other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  select(-Hash) %>% 
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>% 
  group_by(Order, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  # filter(Order=="Calanoida") %>%
  ggplot(., aes(x = Category, y = taxa_sum, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

other %>% 
    summarise_all(sum) %>% 
    mutate(Hash = "other")->other
  
  
# Combine data
fido_18s_s3_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Join with taxa file
fido_18s_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s3_family_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Family = if_else(
    Family == "other" & Order != "other",
    paste0("unidentified ", Order),
    Family
  )) %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s3_save_family_phy

#Save
write.csv(fido_18s_s3_save_family_phy,here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv"))

fido_18s_s3_save_family_phy%>% 
  pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>% 
  group_by(Family, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  ggplot(., aes(x = Category, y = taxa_sum, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Family",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Save a version without Collodaria, normalized to the median # of reads
# Remove "unidentified Collodaria"

fido_18s_s3_save_family_phy_nocollo <- fido_18s_s3_save_family_phy %>%
  filter(Family != "unidentified Collodaria") 

# Save cleaned dataset
write.csv(fido_18s_s3_save_family_phy_nocollo, 
          here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools_nocollodaria.csv"), 
          row.names = FALSE)

#Visualize
fido_18s_s3_save_family_phy_nocollo%>% 
  pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>% 
  group_by(Family, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) %>% 
  ggplot(., aes(x = Category, y = taxa_sum, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Family",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --Visualize the proportion of each taxa across samples
library(patchwork)  # For arranging plots

# Define custom colors for taxa
taxa_colors <- c(
  "other" = "grey",
  "Calanidae" = "#77DD77",       # Pastel green
  "Clausocalanidae" = "#72872d", # Bright pastel green
  "Eucalanidae" = "#89CFF0",     # Baby blue
  "Metridinidae" = "#9370DB",    # Pastel purple
  "Rhincalanidae" = "#4682B4",   # Steel blue
  "Paracalanidae" = "#B0E57C",   # Lime pastel green
  "Oithonidae" = "#f0ca62",      # Light blue
  "Euphausiidae" = "#FF6961",    # Pastel red
  "Salpidae" = "#FFB6C1",        # Pastel pink
  "unidentified Calanoida" = "#2d8087",
  "unidentified Collodaria" = "#FFD580",    # Pastel orange
  "unidentified Siphonophorae" = "#CDA4DE"  # Pastel lavender
)

# Function to process each dataset
process_data <- function(df, dataset_name) {
  df %>%
    pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>%
    group_by(Family, Category) %>%
    summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
    ungroup() %>%
    group_by(Category) %>%
    mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
    mutate(Dataset = dataset_name)
}

# Process datasets separately
plot_data_s1 <- process_data(fido_18s_s1_save_family_phy, "S1")
plot_data_s2 <- process_data(fido_18s_s2_save_family_phy, "S2")
plot_data_s3 <- process_data(fido_18s_s3_save_family_phy, "S3")

# Function to create a proportional stacked bar plot for each dataset
create_plot <- function(data, dataset_name) {
  ggplot(data, aes(x = Category, y = prop, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +  
    scale_fill_manual(values = taxa_colors) +  
    theme_minimal() +
    labs(
      title = paste("Proportional Stacked Bar Plot by Family (", dataset_name, ")", sep = ""),
      x = "Sample",
      y = "Proportion",
      fill = "Taxa"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create separate plots
p1 <- create_plot(plot_data_s1, "S1")
p2 <- create_plot(plot_data_s2, "S2")
p3 <- create_plot(plot_data_s3, "S3")

# Arrange plots in a 3-row grid layout without shared x-axes
final_plot1 <- p1 / p2 / p3  # Stitches them vertically with their own x-axes


final_plot1
# Save first plot
save_path1 <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/CCE_Zooplankton_Metabarcoding_Pub/PCR_bias_correction/figures/miscellaneous/calibration_experiment_taxa_proportions_18s_all_samples")
ggsave(filename = paste0(save_path1, ".pdf"), plot = final_plot1, width = 14, height = 12, dpi = 300)
ggsave(filename = paste0(save_path1, ".png"), plot = final_plot1, width = 14, height = 12, dpi = 300)

# ---- SECOND PLOT: Only Calibration Experiment Samples ----

# Function to filter calibration samples
process_filtered_data <- function(df, dataset_name) {
  df %>%
    select(c(contains("All"), contains("A1"), contains("B3"), contains("C5"), contains("Family"))) %>%
    pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>%
    group_by(Family, Category) %>%
    summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
    ungroup() %>%
    group_by(Category) %>%
    mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum) %>%
    mutate(Dataset = dataset_name)
}

# Process and combine filtered datasets
plot_data_filtered <- bind_rows(
  process_filtered_data(fido_18s_s1_save_family_phy, "S1"),
  process_filtered_data(fido_18s_s2_save_family_phy, "S2"),
  process_filtered_data(fido_18s_s3_save_family_phy, "S3")
)

# Create proportional stacked bar plot (Filtered Calibration Samples)
p2 <- ggplot(plot_data_filtered, aes(x = Category, y = prop, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_fill_manual(values = taxa_colors) +  
  theme_minimal() +
  labs(
    title = "Proportional Stacked Bar Plot by Family (Calibration Experiment)",
    x = "Sample",
    y = "Proportion",
    fill = "Taxa"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p2

# Save second plot
save_path2 <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/CCE_Zooplankton_Metabarcoding_Pub/PCR_bias_correction/figures/miscellaneous/calibration_experiment_taxa_proportions_18s")
ggsave(filename = paste0(save_path2, ".pdf"), plot = p2, width = 14, height = 10, dpi = 300)
ggsave(filename = paste0(save_path2, ".png"), plot = p2, width = 14, height = 10, dpi = 300)




# Extra and Scrap ---------------------------------------------------------

#1) Examine how thresholding impacts composition
# Define the percent thresholds to test
percent_values <- c(0.05, 0.01, 0.005, 0.001)

# Function to filter data by percent threshold
filter_percent_experiment <- function(df, all_cols, experiments, percent_set) {
  taxa_passes <- rep(FALSE, nrow(df))
  
  for (exp in experiments) {
    exp_cols <- all_cols[grep(exp, all_cols)]
    
    if (length(exp_cols) > 0) {
      total_exp_counts <- colSums(df[exp_cols], na.rm = TRUE)
      min_threshold <- sum(total_exp_counts) * percent_set
      
      taxa_passes <- taxa_passes | rowSums(df[exp_cols], na.rm = TRUE) >= min_threshold
    }
  }
  
  return(taxa_passes)
}

# Create an empty list to store results
filtered_data_list <- list()

# Loop over different percent_set values
for (percent_set in percent_values) {
  filtered_taxa <- fido_18s_s1_family_otu %>%
    rownames_to_column("Hash") %>%
    filter(filter_condition(., all_cols) & filter_percent_experiment(., all_cols, experiments, percent_set)) %>%
    mutate(percent_set = as.factor(percent_set))  # Store percent value for plotting
  
  # Store in the list
  filtered_data_list[[as.character(percent_set)]] <- filtered_taxa
}

# Combine all filtered data into one dataframe
filtered_data_combined <- bind_rows(filtered_data_list)

# Count number of families retained at each percentage
family_counts <- filtered_data_combined %>%
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% select(Hash, Family), by = "Hash") %>%
  group_by(percent_set) %>%
  summarise(families_retained = n_distinct(Family), .groups = "drop")

# Join with taxa file and prepare for plotting
plot_data <- filtered_data_combined %>%
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% select(Hash, Family), by = "Hash") %>%
  select(-Hash) %>%
  pivot_longer(cols = -c(Family, percent_set), names_to = "Category", values_to = "Value") %>%
  group_by(Family, percent_set, Category) %>%
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
  ungroup() %>%
  group_by(percent_set, Category) %>%
  mutate(sample_sum = sum(taxa_sum), prop = taxa_sum / sample_sum)

# Merge with family counts for annotations
plot_data <- plot_data %>%
  left_join(family_counts, by = "percent_set")

# Create stacked bar plot
ggplot(plot_data, aes(x = percent_set, y = prop, fill = Family)) +
  geom_bar(stat = "identity") +
  
  # Color Palette
  scale_fill_viridis_d(option = "C", direction = -1) +
  
  # Add text annotations for # of families retained
  geom_text(aes(label = families_retained), 
            stat = "identity", vjust = -1.5, size = 6, fontface = "bold") +
  
  theme_minimal() +
  labs(title = "Taxa Retained Across Different Percent Thresholds",
       x = "Percent Threshold",
       y = "Proportion of Counts",
       fill = "Taxa") +
  
  # Improve aesthetics for publication
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        legend.position = "right")


# #Finally, just for Salpidae, absolute abundance
# # Define custom color for Salpidae
# taxa_colors <- c("Salpidae" = "#FFB6C1")  # Pastel pink for Salpidae
# 
# # Function to process each dataset for Salpidae (absolute abundance)
# process_salpidae_data <- function(df, dataset_name) {
#   df %>%
#     pivot_longer(cols = -Family, names_to = "Category", values_to = "Value") %>%
#     filter(Family == "Salpidae") %>%  # Keep only Salpidae
#     group_by(Family, Category) %>%
#     summarize(taxa_sum = sum(Value), .groups = 'drop') %>%
#     mutate(Dataset = dataset_name)
# }
# 
# # Process datasets separately
# plot_data_s1 <- process_salpidae_data(fido_18s_s1_save_family_phy, "S1")
# plot_data_s2 <- process_salpidae_data(fido_18s_s2_save_family_phy, "S2")
# plot_data_s3 <- process_salpidae_data(fido_18s_s3_save_family_phy, "S3")
# 
# # Function to create an absolute abundance stacked bar plot for Salpidae
# create_salpidae_plot <- function(data, dataset_name) {
#   ggplot(data, aes(x = Category, y = taxa_sum, fill = Family)) +
#     geom_bar(stat = "identity") +  # Absolute abundance
#     scale_fill_manual(values = taxa_colors) +  
#     theme_minimal() +
#     labs(
#       title = paste("Absolute Abundance of Salpidae (", dataset_name, ")", sep = ""),
#       x = "Sample",
#       y = "Absolute Abundance",
#       fill = "Taxa"
#     ) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# }
# 
# # Create separate plots for each dataset
# p1 <- create_salpidae_plot(plot_data_s1, "S1")
# p2 <- create_salpidae_plot(plot_data_s2, "S2")
# p3 <- create_salpidae_plot(plot_data_s3, "S3")
#                            
# 
# # Arrange plots in a 3-row grid layout without shared x-axes
# final_plot3 <- p1 / p2 / p3  # Stitches them vertically with their own x-axes
# final_plot3
# 
# # Save first plot
# save_path3 <- here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/CCE_Zooplankton_Metabarcoding_Pub/PCR_bias_correction/figures/miscellaneous/calibration_experiment_taxa_proportions_18s_salpidae")
# ggsave(filename = paste0(save_path3, ".pdf"), plot = final_plot3, width = 14, height = 12, dpi = 300)
# ggsave(filename = paste0(save_path3, ".png"), plot = final_plot3, width = 14, height = 12, dpi = 300)
