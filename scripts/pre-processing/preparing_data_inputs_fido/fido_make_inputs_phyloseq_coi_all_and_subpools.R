# Script to Make inputs from raw OTU data for fido and PCR bias correction
# Barcode: COI
# 


#Load packages
library (tidyverse)
library (here)
library(ggpubr)
library(fido)
library(phyloseq)
here()



#Read in the OTU data
#Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("data/fido/ASV_table_coi_run1.csv")) %>%
  select(-X) 
#Run2, includes pools
asvcoi_run2=read.csv(here("data/fido/ASV_table_coi_run2.csv")) %>%
  select(-X)



#Phyloseq taxa tables
taxa_coi=read.csv(here("data/taxa_files/blast_metazoo_coi.csv")) %>% 
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
  select(-Subphylum,-Subclass,-Superorder)


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


#Separate out by size
#S1
fido_coi_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_coi_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) %>% 
  filter(rowSums(.) != 0)
fido_coi_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) %>% 
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
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

#Make phyloseq objects
fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1)
fido_coi_s2_phy=phyloseq(fido_coi_s2_otu,taxcoi_s2)
fido_coi_s3_phy=phyloseq(fido_coi_s3_otu,taxcoi_s3)



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
write.csv(taxcoi_genus,here("data/phyloseq_bio_data/COI/fido_coi_genus_tax_table.csv"))



#--------------------------------------------------------------------------------------#


# PART 2: Make OTU and Taxa Inputs for PCR Bias Correction ----------------
# This part will only include the taxa that meet the conditions for the PCR bias correction analysis which include:
# 1. Presence in at least one replicate per treatment per experiment
# 2. The sum of reads in the 'other' category doesn't exceed some threshold, for COI I had to be much looser on this to retain data,
# here I used 55% 


## ==== S1 ====
# Separate rows based appearance in the calibration samples
#Pools
experiments <- c("AllPool", "Pooled.A1.A3", "Pooled.B3.B5", "Pooled.C5.C7")

#PCR Cycles
treatments <- c("20C", "24C", "28C")

#TEechnical Replicates
replicates <- c("1", "2", "3")

# Function to create column names based on treatment, experiment, and replicate
valid_cols <- function(treatment, experiment, replicate) {
  paste0(treatment, "_", experiment, ".", replicate)
}

# Generate all possible column names for the conditions
all_cols <- lapply(experiments, function(exp) {
  sapply(treatments, function(treat) {
    sapply(replicates, function(rep) {
      valid_cols(treat, exp, rep)
    })
  })
})

# Flatten the list of column names
all_cols <- unlist(all_cols)

# Filter out any non-existent columns from my dataframe (Relevant for failed 20C)
all_cols <- all_cols[all_cols %in% names(fido_coi_s1_genus_otu)]

# Define a condition to check presence in at least one replicate per treatment per experiment
filter_condition <- function(df, cols) {
  # Applying row-wise check to ensure at least one non-zero entry per experiment-treatment combination
  all(sapply(split(cols, gsub("\\..*$", "", cols)), function(c) {
    any(rowSums(df[c] > 0, na.rm = TRUE) > 0)
  }))
}

# Apply the filter across the dataframe
fido_taxa_filt <- fido_coi_s1_genus_otu%>%
  rownames_to_column("Hash") %>%
  rowwise() %>%
  filter(filter_condition(cur_data(), all_cols)) %>% as.data.frame()%>% 
  left_join(fido_coi_s1_genus_taxa %>% rownames_to_column("Hash") , by = "Hash") %>% 
  mutate(across(c(Genus, Family, Hash), ~ ifelse(is.na(Genus), "unidentified Calanoida", .))) %>% 
  select(-Phylum:-Species)



other <- fido_coi_s1_genus_otu %>%
  anti_join(fido_taxa_filt) 


#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s1_genus_taxa %>% rownames_to_column("Hash") %>% 
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


# Add taxa file and make an other and unidentified calanoida category
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s1_genus_taxa %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  mutate(Category = ifelse(Order == "Calanoida", "unidentified Calanoida", "other")) %>%
  mutate(Category = ifelse(is.na(Category), "other", Category)) %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(Hash=Category) %>% 
  select(-Category)->other



#Find columns where other exceeds 5% of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% filter(Hash!="unidentified Calanoida")%>% select(-Hash) )
sums_fido <- colSums(fido_coi_s1_genus_otu)

# Step 2: Calculate 50% of the column sums of fido_coi_s1_genus_otu
thresh_fraction=0.5
thresholds <- sums_fido * thresh_fraction

# Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_coi_s3_genus_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_coi_s1_genus_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))


# Combine data
fido_coi_s1_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Check to ensure the values are conserved
colSums(fido_coi_s1_genus_otu)[1:5]
colSums(fido_coi_s1_final)[1:5]


#Join with taxa file
fido_coi_s1_save_genus_phy <- fido_coi_s1_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_coi_s1_genus_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
  mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
  select(-Phylum, -Class, -Family, -Order, -Hash) %>%
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


#Save
write.csv(fido_coi_s1_save_genus_phy,here("data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv"))

#Refresh before next size
rm(other)


## ==== S2 ====
# Filter out any non-existent columns from my dataframe
all_cols <- all_cols[all_cols %in% names(fido_coi_s2_genus_otu)]


# Apply the filter across the dataframe
fido_taxa_filt <- fido_coi_s2_genus_otu%>%
  rownames_to_column("Hash") %>%
  rowwise() %>%
  filter(filter_condition(cur_data(), all_cols)) %>% as.data.frame()%>% 
  left_join(fido_coi_s2_genus_taxa %>% rownames_to_column("Hash") , by = "Hash") %>% 
  mutate(across(c(Genus, Family, Hash), ~ ifelse(is.na(Genus), "unidentified Calanoida", .))) %>% 
  select(-Phylum:-Species)



other <- fido_coi_s2_genus_otu %>%
  anti_join(fido_taxa_filt) 


#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s2_genus_taxa %>% rownames_to_column("Hash") %>% 
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


# Add taxa file and make an other and unidentified calanoida category
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s2_genus_taxa %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  mutate(Category = ifelse(Order == "Calanoida", "unidentified Calanoida", "other")) %>%
  mutate(Category = ifelse(is.na(Category), "other", Category)) %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(Hash=Category) %>% 
  select(-Category)->other



#Find columns where other exceeds 5% of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% filter(Hash!="unidentified Calanoida")%>% select(-Hash) )
sums_fido <- colSums(fido_coi_s2_genus_otu)

# Step 2: Calculate 5% of the column sums of fido_coi_s2_genus_otu
thresh_fraction=0.5
thresholds <- sums_fido * thresh_fraction

# Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_coi_s3_genus_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_coi_s3_genus_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))


# Combine data
fido_coi_s2_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

colSums(fido_coi_s2_genus_otu)[1:5]
colSums(fido_coi_s2_final)[1:5]


#Join with taxa file
fido_coi_s2_save_genus_phy <- fido_coi_s2_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_coi_s2_genus_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
  mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
  select(-Phylum, -Class, -Family, -Order, -Hash) %>%
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

#Save
write.csv(fido_coi_s2_save_genus_phy,here("data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv"))



## ==== S3 ====

# Filter out any non-existent columns from my dataframe
all_cols <- all_cols[all_cols %in% names(fido_coi_s3_genus_otu)]
# Apply the filter across the dataframe
fido_taxa_filt <- fido_coi_s3_genus_otu%>%
  rownames_to_column("Hash") %>%
  rowwise() %>%
  filter(filter_condition(cur_data(), all_cols)) %>% as.data.frame()%>% 
  left_join(fido_coi_s3_genus_taxa %>% rownames_to_column("Hash") , by = "Hash") %>% 
  mutate(across(c(Genus, Family, Hash), ~ ifelse(is.na(Genus), "unidentified Calanoida", .))) %>% 
  select(-Phylum:-Species)



other <- fido_coi_s3_genus_otu %>%
  anti_join(fido_taxa_filt) 


#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s3_genus_taxa %>% rownames_to_column("Hash") %>% 
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


# Add taxa file and make an other and unidentified calanoida category
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(fido_coi_s3_genus_taxa %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  mutate(Category = ifelse(Order == "Calanoida", "unidentified Calanoida", "other")) %>%
  mutate(Category = ifelse(is.na(Category), "other", Category)) %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(Hash=Category) %>% 
  select(-Category)->other



#Find columns where other exceeds 5% of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% filter(Hash!="unidentified Calanoida")%>% select(-Hash) )
sums_fido <- colSums(fido_coi_s3_genus_otu)

# Step 2: Calculate 5% of the column sums of fido_coi_s3_genus_otu
thresh_fraction=0.5
thresholds <- sums_fido * thresh_fraction

# Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_coi_s3_genus_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_coi_s3_genus_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))


# Combine data
fido_coi_s3_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

colSums(fido_coi_s3_genus_otu)[1:5]
colSums(fido_coi_s3_final)[1:5]


#Join with taxa file
fido_coi_s3_save_genus_phy <- fido_coi_s3_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_coi_s3_genus_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
  mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
  select(-Phylum, -Class, -Family, -Order, -Hash) %>%
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))



#Save
write.csv(fido_coi_s3_save_genus_phy,here("data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv"))
