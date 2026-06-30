# From these runs, we didn't get data from the 20C or various pools for COI runs, but it worked for 18s
# try first with 18s

library (tidyverse)
library (here)
library(ggpubr)
library(fido)
library(phyloseq)
here()



###18S
#Read in the OTU data
#Run 1 (Non pooled data)
asv18s_run1=read.csv(here("data/past/","ASV_table_18s_run1.csv")) %>%
  select(-X) 
#Run2
asv18s_run2=read.csv(here("data/past/","ASV_table_18s_run2.csv")) %>%
  select(-X)



#Taxa Tables 
taxa_18s=read.csv(here("data/taxa_files/blast_metazoo_18s.csv")) %>% 
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
  dplyr::select(c(contains("All"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 
fido_18s_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) %>% 
  filter(rowSums(.) != 0)
fido_18s_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)


###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa tables
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Calanoida", "unidentified Calanoida", Family))
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Calanoida", "unidentified Calanoida", Family))
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))%>%
  mutate(Family = ifelse(is.na(Family), Order, Family)) %>%
  mutate(Family = ifelse(Family == "", 'other', Family)) %>% 
  mutate(Family = ifelse(Family == "" & Order == "Calanoida", "unidentified Calanoida", Family)) %>% 
  mutate(Family = ifelse(Family == "other" & Order == "Calanoida", "unidentified Calanoida", Family))
tax18s_s3=  tax_table(as.matrix(tax18s_s3))






#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
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
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3, metadata)
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
  unique() 
write.csv(tax18s_family,here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))



## ==== S1 ====
# Separate rows based appearance in the calibration samples
fido_taxa_filt <- fido_18s_s1_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_18s_s1_family_otu %>%
  anti_join(fido_18s_s1_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 


#Check composition before proceeding
filter_other=1

if(filter_other==1){
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
  
  #Find columns where other exceeds 5% of the sample total
  # Step 1: Calculate column sums for both DataFrames
  sums_other <- colSums(other %>% select(-Hash))
  sums_fido <- colSums(fido_18s_s1_family_otu)
  
  # Step 2: Calculate 5% of the column sums of fido_18s_s3_family_otu
  thresholds <- sums_fido * 0.05
  
  # Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_18s_s3_family_otu`
  columns_to_remove <- names(which(sums_other > thresholds))
  
  # Step 4: Remove these columns from the fido_18s_s3_family_otu DataFrame
  fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
  other=other%>% select(-all_of(columns_to_remove))
  
  
  # Combine data
  fido_18s_s1_final <- rbind(fido_taxa_filt,other)  %>% 
    group_by(Hash) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
    column_to_rownames("Hash")
  
  colSums(fido_18s_s1_family_otu)[1:5]
  colSums(fido_18s_s1_final)[1:5]
  
  
  #Join with taxa file
  fido_18s_s1_save_family_phy <- fido_18s_s1_final %>%
    rownames_to_column("Hash") %>%
    left_join(fido_18s_s1_family_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
    mutate(Family = ifelse(Hash == "other", "other", Family)) %>% 
    select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
    group_by(Family) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  
  #Save
  write.csv(fido_18s_s1_save_family_phy,here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"))
}else{
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
  
  # Combine data
  fido_18s_s1_final <- rbind(fido_taxa_filt,other)  %>% 
    group_by(Hash) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
    column_to_rownames("Hash")
  
  colSums(fido_18s_s1_family_otu)[1:5]
  colSums(fido_18s_s1_final)[1:5]
  
  
  #Join with taxa file
  fido_18s_s1_save_family_phy <- fido_18s_s1_final %>%
    rownames_to_column("Hash") %>%
    left_join(fido_18s_s1_family_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
    mutate(Family = ifelse(Hash == "other", "other", Family)) %>% 
    select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
    group_by(Family) %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  
  #Save
  write.csv(fido_18s_s1_save_family_phy,here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"))
  
}

rm(other)


## ==== s2 ====
fido_taxa_filt <- fido_18s_s2_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")

other <- fido_18s_s2_family_otu %>%
  anti_join(fido_18s_s2_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 


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
  filter(Order=="Calanoida") %>%
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
sums_fido <- colSums(fido_18s_s2_family_otu)

# Step 2: Calculate 5% of the column sums of fido_18s_s3_family_otu
thresholds <- sums_fido * 0.1

# Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_18s_s3_family_otu`
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
  mutate(Family = ifelse(Hash == "other", "other", Family)) %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s2_save_family_phy

#Save
write.csv(fido_18s_s2_save_family_phy,here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"))





## ==== s3 ====
fido_taxa_filt <- fido_18s_s3_family_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_18s_s3_family_otu %>%
  anti_join(fido_18s_s3_family_otu %>%
              filter(rowSums(select(., 1:9) == 0) <= 2)) 

fido_18s_s3_family_otu %>% 
  summarise_all(sum) %>% 
  mutate(Order = "total")->total
  
#Check composition before proceeding

#The other here would be Orders that weren't identified to that level, as with the NA
other %>% 
  rownames_to_column("Hash")%>% 
  left_join(taxa_18s %>% rownames_to_column("Hash") %>% 
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
sums_fido <- colSums(fido_18s_s3_family_otu)

# Step 2: Calculate 5% of the column sums of fido_18s_s3_family_otu
thresholds <- sums_fido * 0.1

# Step 3: Identify columns where the sum of `other` is greater than 5% of the sum of `fido_18s_s3_family_otu`
columns_to_remove <- names(which(sums_other > thresholds))

# Step 4: Remove these columns from the fido_18s_s3_family_otu DataFrame
fido_taxa_filt <- fido_taxa_filt %>% select(-all_of(columns_to_remove))
other=other%>% select(-all_of(columns_to_remove))



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
  mutate(Family = ifelse(Hash == "other", "other", Family)) %>% 
  select(-Phylum, -Class, -Genus, -Order, -Species, -Hash) %>%
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s3_save_family_phy

#Save
write.csv(fido_18s_s3_save_family_phy,here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"))

