#Preparing inputs for PCR bias-mitigation for COI

library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)
here()



#Read in the OTU data
#Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("data/fido/ASV_table_coi_run1.csv")) %>%
  select(-X) 
#Run2
asvcoi_run2=read.csv(here("data/fido/ASV_table_coi_run2.csv")) %>%
  select(-X)



#Taxa tables
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
  column_to_rownames("Hash")


#Format Long
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

#Replace X
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

###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_coi_s1_otu=fido_coi_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s2_otu=fido_coi_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_coi_s3_otu=fido_coi_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
taxcoi_s1 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s1))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus))
taxcoi_s1=  tax_table(as.matrix(taxcoi_s1))

#S2
taxcoi_s2 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s2_otu))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus))
taxcoi_s2=  tax_table(as.matrix(taxcoi_s2))
#S3
taxcoi_s3 = taxa_coi %>% filter(rownames(taxa_coi) %in% rownames(fido_coi_s3_otu))%>%
  mutate(Genus = ifelse(is.na(Genus), Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "", 'other', Genus)) %>% 
  mutate(Genus = ifelse(Genus == "" & Order == "Calanoida", "unidentified Calanoida", Genus)) %>% 
  mutate(Genus = ifelse(Genus == "other" & Order == "Calanoida", "unidentified Calanoida", Genus))
taxcoi_s3=  tax_table(as.matrix(taxcoi_s3))






#Metadata
metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1)
fido_coi_s2_phy=phyloseq(fido_coi_s2_otu,taxcoi_s2)
fido_coi_s3_phy=phyloseq(fido_coi_s3_otu,taxcoi_s3)



#PHYLOSEQ
#Agglomerate at the genus level

#S1
fido_coi_s1_phy=phyloseq(fido_coi_s1_otu,taxcoi_s1, metadata)
fido_coi_s1_genus=tax_glom(fido_coi_s1_phy, taxrank = "Genus")


#Check colsums
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
difference==0

#Make inputs for filtering
fido_coi_s1_genus_otu=otu_table(fido_coi_s1_genus) %>% as.data.frame() 

fido_coi_s1_genus_taxa=tax_table(fido_coi_s1_genus) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_coi_s1_genus_otu <- bind_rows(fido_coi_s1_genus_otu, difference)

#Check
colSums(fido_coi_s1_otu)==colSums(fido_coi_s1_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s1_genus_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_coi_s1_genus_taxa) -> fido_coi_s1_genus_taxa



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

fido_coi_s2_genus_taxa=tax_table(fido_coi_s2_genus) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_coi_s2_genus_otu <- bind_rows(fido_coi_s2_genus_otu, difference)

#Check
colSums(fido_coi_s2_otu)==colSums(fido_coi_s2_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s2_genus_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
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

fido_coi_s3_genus_taxa=tax_table(fido_coi_s3_genus) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_coi_s3_genus_otu <- bind_rows(fido_coi_s3_genus_otu, difference)

#Check
colSums(fido_coi_s3_otu)==colSums(fido_coi_s3_genus_otu)

#Need to add 'other' row to taxa table to accomodate this new category
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_coi_s3_genus_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_coi_s3_genus_taxa) -> fido_coi_s3_genus_taxa



#Save aglomerated genus taxa file, replace all columns with 'other' where genus is 'other'
taxcoi_genus=rbind(fido_coi_s1_genus_taxa,fido_coi_s2_genus_taxa,fido_coi_s3_genus_taxa) %>%
  rownames_to_column("Hash") %>% 
  select(-Subphylum,-Subclass,-Superorder,-Species,-Hash) %>% 
  mutate(
    Phylum = if_else(Genus == 'other', 'other', Phylum),
    # Subphylum = if_else(Genus == 'other', 'other', Subphylum),
    Class = if_else(Genus == 'other', 'other', Class),
    # Subclass = if_else(Genus == 'other', 'other', Subclass),
    # Superorder = if_else(Genus == 'other', 'other', Superorder),
    Order = if_else(Genus == 'other', 'other', Order),
    Family = if_else(Genus == 'other', 'other', Family)
  ) %>% 
  unique() 
write.csv(taxcoi_genus,here("data/phyloseq_bio_data/coi/fido_coi_genus_tax_table.csv"))




#Visualize ECDF
fido_coi_s1_genus_otu[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

## ==== S1 ====
# Separate rows based appearance in the calibration samples
fido_taxa_filt <- fido_coi_s1_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_coi_s1_genus_otu %>%
  anti_join(fido_coi_s1_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 


#Check composition before proceeding

  #The other here would be Orders that weren't identified to that level, as with the NA
  other %>% 
    rownames_to_column("Hash")%>% 
    left_join(taxa_coi %>% rownames_to_column("Hash") %>% 
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
  
  #Find columns where other exceeds 5% of the sample total
  # Step 1: Calculate column sums for both DataFrames
  sums_other <- colSums(other %>% select(-Hash))
  sums_fido <- colSums(fido_coi_s1_genus_otu)
  
  # Step 2: Calculate 5% of the column sums of fido_coi_s3_genus_otu
  thresholds <- sums_fido * 0.9
  
  # Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_coi_s3_genus_otu`
  columns_to_remove <- names(which(sums_other > thresholds))
  
  # Step 4: Remove these columns from the fido_coi_s3_genus_otu DataFrame
  fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
  other=other%>% select(-all_of(columns_to_remove))
  
  
  # Combine data
  fido_coi_s1_final <- rbind(fido_taxa_filt,other)  %>% 
    group_by(Hash) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
    column_to_rownames("Hash")
  
  colSums(fido_coi_s1_genus_otu)[1:5]
  colSums(fido_coi_s1_final)[1:5]
  
  
  #Join with taxa file
  fido_coi_s1_save_genus_phy <- fido_coi_s1_final %>%
    rownames_to_column("Hash") %>%
    left_join(fido_coi_s1_genus_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
    mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
    select(-Phylum, -Class, -Family, -Order, -Species, -Hash) %>%
    group_by(Genus) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  
  #Save
  write.csv(fido_coi_s1_save_genus_phy,here("data/fido/phy/fido_coi_s1_ecdf_genus_phy.csv"))


rm(other)


## ==== s2 ====
fido_taxa_filt <- fido_coi_s2_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")

other <- fido_coi_s2_genus_otu %>%
  anti_join(fido_coi_s2_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 


#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% 
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

#Find columns where other exceeds 5% of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% select(-Hash))
sums_fido <- colSums(fido_coi_s2_genus_otu)

# Step 2: Calculate 5% of the column sums of fido_coi_s3_genus_otu
thresholds <- sums_fido * 0.99

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

#Join with taxa file
fido_coi_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_coi_s2_genus_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
  select(-Phylum, -Class, -Family, -Order, -Species, -Hash) %>%
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_coi_s2_save_genus_phy

#Save
write.csv(fido_coi_s2_save_genus_phy,here("data/fido/phy/fido_coi_s2_ecdf_genus_phy.csv"))





## ==== s3 ====
fido_taxa_filt <- fido_coi_s3_genus_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_coi_s3_genus_otu %>%
  anti_join(fido_coi_s3_genus_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 

fido_coi_s3_genus_otu %>% 
  summarise_all(sum) %>% 
  mutate(Order = "total")->total

#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_coi %>% rownames_to_column("Hash") %>% 
              select(Hash,Order), by = "Hash") %>% 
  select(-Hash) %>% 
  rbind(total) %>% 
  pivot_longer(cols = -Order, names_to = "Category", values_to = "Value") %>% 
  group_by(Order, Category) %>% 
  summarize(taxa_sum = sum(Value), .groups = 'drop') %>% 
  ungroup() %>%  group_by(Category) %>% 
  mutate(sample_sum=sum(taxa_sum),prop=taxa_sum/sample_sum) ->other_for_plot

other_for_plot %>% 
  filter(Order != "total") %>% 
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

#Find columns where other exceeds 5% of the sample total
# Step 1: Calculate column sums for both DataFrames
sums_other <- colSums(other %>% select(-Hash))
sums_fido <- colSums(fido_coi_s3_genus_otu)

# Step 2: Calculate 5% of the column sums of fido_coi_s3_genus_otu
thresholds <- sums_fido * 0.99

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

#Join with taxa file
fido_coi_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_coi_s3_genus_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Genus = ifelse(Hash == "other", "other", Genus)) %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_coi_s3_save_genus_phy

#Save
write.csv(fido_coi_s3_save_genus_phy,here("data/fido/phy/fido_coi_s3_ecdf_genus_phy.csv"))


