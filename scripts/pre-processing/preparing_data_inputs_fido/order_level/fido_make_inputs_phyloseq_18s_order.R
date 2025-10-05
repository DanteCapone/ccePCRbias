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
taxa_18s=read.csv(here("data/taxa_files/blast_metazoo_18s.csv"))%>%
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

# MPN: Just to clarify, what is the difference between the ".1" and ".2" samples?
# E.g., C1.T7.H9.1 versus C1.T7.H9.2?


###Phyloseq filtering: Use phyloseq for filtering and agglomerating

fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa tables
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))%>%
  mutate(Order = ifelse(is.na(Order), Class, Order)) %>% 
  mutate(Order = ifelse(Order=="", 'other', Order))
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))%>%
  mutate(Order = ifelse(is.na(Order), Class, Order)) %>% 
  mutate(Order = ifelse(Order=="", 'other', Order))
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))%>%
  mutate(Order = ifelse(is.na(Order), Class, Order)) %>% 
  mutate(Order = ifelse(Order=="", 'other', Order))
tax18s_s3=  tax_table(as.matrix(tax18s_s3))






#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1)
fido_18s_s2_phy=phyloseq(fido_18s_s2_otu,tax18s_s2)
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3)



# Agglomerate at the Order Level -----------------------------------------

# S1 ----------------------------------------------------------------------

fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1, metadata)
fido_18s_s1_order=tax_glom(fido_18s_s1_phy, taxrank = "Order")

#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums(otu_table(fido_18s_s1))
colsums_before[1:5]
# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s1_order))
colsums_after[1:5]
colsums_before == colsums_after
#Soem columns are different
# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s1_order_otu=otu_table(fido_18s_s1_order) %>% as.data.frame() 

fido_18s_s1_order_taxa=tax_table(fido_18s_s1_order) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s1_order_otu <- bind_rows(fido_18s_s1_order_otu, difference)


colSums(fido_18s_s1_otu)[1:5]
colSums(fido_18s_s1_order_otu)[1:5]

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s1_order_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s1_order_taxa) -> fido_18s_s1_order_taxa



# S2 ----------------------------------------------------------------------
fido_18s_s2_phy=phyloseq(fido_18s_s2_otu,tax18s_s2, metadata)
fido_18s_s2_order=tax_glom(fido_18s_s2_phy, taxrank = "Order")

#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s2))
colsums_before[1:5]
# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s2_order))
colsums_after[1:5]
# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s2_order_otu=otu_table(fido_18s_s2_order) %>% as.data.frame() 

fido_18s_s2_order_taxa=tax_table(fido_18s_s2_order) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s2_order_otu <- bind_rows(fido_18s_s2_order_otu, difference)


colSums(fido_18s_s2_otu)[1:5]
colSums(fido_18s_s2_order_otu)[1:5]

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s2_order_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s2_order_taxa) -> fido_18s_s2_order_taxa


# S3 ----------------------------------------------------------------------
fido_18s_s3_phy=phyloseq(fido_18s_s3_otu,tax18s_s3, metadata)
fido_18s_s3_order=tax_glom(fido_18s_s3_phy, taxrank = "Order")

#Check colsums
# Calculate column sums before tax glomming
colsums_before <- colSums((fido_18s_s3))
colsums_before[1:5]
# Calculate column sums after tax glomming
colsums_after <- colSums(otu_table(fido_18s_s3_order))
colsums_after[1:5]
# Find the difference
difference <- colsums_before - colsums_after %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Hash="other") %>% 
  column_to_rownames("Hash")


#Make inputs for filtering
fido_18s_s3_order_otu=otu_table(fido_18s_s3_order) %>% as.data.frame() 

fido_18s_s3_order_taxa=tax_table(fido_18s_s3_order) %>% as.data.frame() 

# Add the new difference row to the dataframe
fido_18s_s3_order_otu <- bind_rows(fido_18s_s3_order_otu, difference)


colSums(fido_18s_s3_otu)[1:5]
colSums(fido_18s_s3_order_otu)[1:5]

#Need to add 'other' row to taxa table
data.frame(
  row_name = "other",
  stringsAsFactors = FALSE,
  lapply(fido_18s_s3_order_taxa, function(x) "other")
) %>% column_to_rownames("row_name") %>%
  rbind(.,fido_18s_s3_order_taxa) -> fido_18s_s3_order_taxa

#Save aglomerated order taxa file, replace all columns with 'other' where order is 'other'
tax18s_order=rbind(fido_18s_s1_order_taxa,fido_18s_s2_order_taxa,fido_18s_s3_order_taxa) %>%
  unique() %>%
  mutate(
    Phylum = if_else(Order == 'other', 'other', Phylum),
    Subphylum = if_else(Order == 'other', 'other', Subphylum),
    Class = if_else(Order == 'other', 'other', Class),
    Subclass = if_else(Order == 'other', 'other', Subclass),
    Superorder = if_else(Order == 'other', 'other', Superorder),
  ) %>% 
  rownames_to_column("Hash") %>% 
  select(-Species, -Genus, -Family,-Hash,-Superorder,-Subclass,-Subphylum) %>% 
  distinct() 
write.csv(tax18s_order,here("data/phyloseq_bio_data/18S/fido_18s_order_tax_table.csv"))




# Filtration based on Calibration Samples ---------------------------------


## ==== S1 ====
# Separate rows based appearance in the calibration samples
fido_taxa_filt <- fido_18s_s1_order_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_18s_s1_order_otu %>%
  anti_join(fido_18s_s1_order_otu %>%
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
  summarize(Total = sum(Value), .groups = 'drop') ->a
  ggplot(., aes(x = Category, y = Total, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> other_check_p1
other_check_p1

other %>% 
  summarise_all(sum) %>% 
  mutate(Hash = "other")->other
  
  
  

# Combine data
fido_18s_s1_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

colSums(fido_18s_s1_order_otu)[1:5]
colSums(fido_18s_s1_final)[1:5]


#Join with taxa file
fido_18s_s1_save_order_phy <- fido_18s_s1_final %>%
  rownames_to_column("Hash") %>%
  left_join(fido_18s_s1_order_taxa %>% rownames_to_column("Hash"), by = "Hash") %>% 
  mutate(Order = ifelse(Hash == "other", "other", Order)) %>% 
  select(-Phylum, -Class, -Genus, -Family, -Species, -Hash) %>%
  group_by(Order) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))


#Save
write.csv(fido_18s_s1_save_order_phy,here("data/fido/phy/fido_18s_s1_ecdf_order_phy.csv"))



## ==== s2 ====
fido_taxa_filt <- fido_18s_s2_order_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")

other <- fido_18s_s2_order_otu %>%
  anti_join(fido_18s_s1_order_otu %>%
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
  summarize(Total = sum(Value), .groups = 'drop') %>% 
  ggplot(., aes(x = Category, y = Total, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> other_check_p2
other_check_p2

other %>% 
  summarise_all(sum) %>% 
  mutate(Hash = "other")->other

# Combine data
fido_18s_s2_final <- rbind(fido_taxa_filt,other)  %>% 
  group_by(Hash) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Hash")

#Join with taxa file
fido_18s_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(fido_18s_s2_order_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Order = ifelse(Hash == "other", "other", Order)) %>% 
  select(-Phylum, -Class, -Genus, -Family, -Species, -Hash) %>%
  group_by(Order) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s2_save_order_phy

#Save
write.csv(fido_18s_s2_save_order_phy,here("data/fido/phy/fido_18s_s2_ecdf_order_phy.csv"))

rm(other)



## ==== s3 ====
fido_taxa_filt <- fido_18s_s3_order_otu %>% filter(rowSums(select(., 1:9) == 0) <= 2) %>%
  rownames_to_column("Hash")


other <- fido_18s_s3_order_otu %>%
  anti_join(fido_18s_s1_order_otu %>%
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
  summarize(Total = sum(Value), .groups = 'drop') %>% 
  ggplot(., aes(x = Category, y = Total, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot by Order",
       x = "Order",
       y = "Sum of Values",
       fill = "Category")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> other_check_p3


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
  left_join(fido_18s_s3_order_taxa %>% rownames_to_column("Hash"), by="Hash")%>%
  mutate(Order = ifelse(Hash == "other", "other", Order)) %>% 
  select(-Phylum, -Class, -Genus, -Family, -Species, -Hash) %>%
  group_by(Order) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE))->fido_18s_s3_save_order_phy

#Save
write.csv(fido_18s_s3_save_order_phy,here("data/fido/phy/fido_18s_s3_ecdf_order_phy.csv"))


# Visualize Other composition

other_check_p1
other_check_p2
other_check_p3
