# From these runs, we didn't get data from the 20C or various pools for COI runs, but it worked for 18s
# try first with 18s

library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
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
taxa_18s=read.csv(here("data/past/metazoopruned18s_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>%
  column_to_rownames("Hash")

# 2) Merging and manipulation (updated 8/24/2023 to create a new 18S input for fido where
# I don't average technical replicates)

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

##MPN: Another common method is to amalgamate to the genus level and add anything that doesn't fit to an "other" category. This is very common in microbiome.

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



###DC: 12/14/2023-Use phyloseq for filtering and agglomerating
##Phyloseq filtering
fido_18s_s1_otu=fido_18s_s1 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu=fido_18s_s2 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu=fido_18s_s3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
tax18s_s1 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1))
tax18s_s1=  tax_table(as.matrix(tax18s_s1))

#S2
tax18s_s2 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu))
tax18s_s2=  tax_table(as.matrix(tax18s_s2))
#S3
tax18s_s3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu))
tax18s_s3=  tax_table(as.matrix(tax18s_s3))

#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1)


#=== 2/12/2024: Make an input for the raw read analysis

#S1
fido_18s_s1 %>%
  rownames_to_column("Hash") %>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) -> fido_18s_s1_merged 

#Save
write.csv(fido_18s_s1_merged,here("data/raw_reads/otus_18_prefilt_s1.csv"))

sum#S2
fido_18s_s2 %>%
  rownames_to_column("Hash") %>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)-> fido_18s_s2_merged 

#Save
write.csv(fido_18s_s2_merged,here("data/raw_reads/otus_18_prefilt_s2.csv"))

#s3
fido_18s_s3 %>%
  rownames_to_column("Hash") %>%
  pivot_longer(cols = -Hash, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Hash) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0) -> fido_18s_s3_merged 

#Save
write.csv(fido_18s_s3_merged,here("data/raw_reads/otus_18_prefilt_s3.csv"))


#All
merge(fido_18s_s1_merged, fido_18s_s2_merged, by = "Hash", all = TRUE) %>%
  merge(.,fido_18s_s3_merged, by = "Hash", all = TRUE)%>%
  column_to_rownames("Hash") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_all

write.csv(fido_18s_merged_all,here("data/raw_reads/otu_18s_prefilt_all.csv"))

#====== End additional code======#


#PHYLOSEQ
#Agglomerate at the genus level
fido_18s_s1_phy=phyloseq(fido_18s_s1_otu,tax18s_s1, metadata)
fido_18s_s1_genus=tax_glom(fido_18s_s1_phy, taxrank = "Family") 
fido_18s_s1_genus_otu=otu_table(fido_18s_s1_genus) %>% as.data.frame()
fido_18s_s1_genus_taxa=tax_table(fido_18s_s1_genus) %>% as.data.frame()





# Function to filter taxa based on ECDF
filter_taxa_ecdf <- function(physeq, quantile_threshold) {
  # Extract abundance data
  abundances <- otu_table(physeq)
  
  # Check if abundances matrix is empty
  if (nrow(abundances) == 0) {
    stop("Abundance matrix is empty. No taxa to filter.")
  }
  
  # Compute ECDF for each taxon
  ecdf_values <- apply(abundances, 1, ecdf)
  
  # Identify taxa exceeding the quantile_threshold
  taxa_to_keep <- rownames(abundances)[apply(ecdf_values, 1, function(x) x(quantile_threshold)) > 0]
  
  # Create a new category for taxa below the threshold
  taxa_to_aggregate <- rownames(abundances)[apply(ecdf_values, 1, function(x) x(quantile_threshold)) == 0]
  
  # Aggregate taxa below threshold into a single "Other" category
  other_abundance <- rowSums(abundances[taxa_to_aggregate, ])
  
  # Combine filtered taxa and "Other" category
  if (length(taxa_to_aggregate) > 0) {
    other_taxa <- tax_table(physeq)[taxa_to_aggregate, ]
    other_taxa$Genus <- "Other"
    other_taxa <- tax_table(other_taxa)
    other_otu <- otu_table(other_abundance)
    physeq_filtered <- merge_phyloseq(physeq, phyloseq(other_otu, other_taxa))
  } else {
    physeq_filtered <- physeq
  }
  
  # Filter taxa above the threshold
  physeq_filtered <- prune_taxa(taxa_to_keep, physeq_filtered)
  
  return(physeq_filtered)
}

# Set the quantile threshold (e.g., 0.95)
quantile_threshold <- 0.9

# Apply the filter function
physeq_filtered <- filter_taxa_ecdf(fido_18s_s1_genus, quantile_threshold)





##MPN: Why not just use phyloseq?
#Set ECDF threshold
thresh_val=0.9


#9/26/2023
#ECDF plot for determining criteria 

#90% threshold
#S1

## MPN: Do you want to agglomerate the taxa at say the genus level? If so, should do before filtering.
##MPN: to use the ecdf method, we want to look at the plot. Is this where you came up with the 90% threshold?
##MPN: also might want to 
fido_18s_s1[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)
###end of added code

fido_18s_s1=fido_18s_s1 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s1$rowsum, thresh_val)

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
above_threshold <- fido_18s_s1 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s1 %>% 
  ##MPN: are you throwing out otus that don't appear in the pooled samples?
  filter(rowsum < threshold) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s1_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
##MPN: are the counts of any of these large? Why are they missing?
rows_not_in_taxa_18s <- rownames(fido_18s_s1_final)[!rownames(fido_18s_s1_final) %in% taxa_18s$Hash]

#Check if in the other classification
taxa_18s2=read.csv(here("data/BLAST_taxa_class/zhang_taxa.csv"))  %>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>%
  column_to_rownames("Hash")

rows_only_in_taxa_18s2 <- taxa_18s2 %>% filter(rownames(taxa_18s2) %in% rows_not_in_taxa_18s)


 #Add to taxa 18s df


for(row in rows_not_in_taxa_18s) {
  if(row != "other") {
    fido_18s_s1_final["other",] <- fido_18s_s1_final["other",] + fido_18s_s1_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_18s_s1_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s1_save_spp_hash

#Hash
  select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s1_save_hash
  
  #Save
  write.csv(fido_18s_s1_save_spp_hash,here("data/fido/fido_18s_s1_ecdf_spp_hash.csv"))
  write.csv(fido_18s_s1_save_hash,here("data/fido/fido_18s_s1_ecdf_hash.csv"))

#Save
# write.csv(fido_18s_s1_final,here("data/EDCF/18s/fido_18s_s1_ecdf_spp_hash.csv"))



#Plot ECDF
fido_18s_s1 %>%  ggplot(., aes(rowSums(.))) +
  stat_ecdf(geom = "point")+
  geom_vline(aes(xintercept = threshold), linetype = "dashed", color = "red") +
  labs(title = "ECDF of Row Sums",
       x = "Row Sum",
       y = "ECDF") +
  theme_minimal()






### ========== S2 ===============
#s2
fido_18s_s2[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)
###end of added code

fido_18s_s2=fido_18s_s2 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s2$rowsum, thresh_val)

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
above_threshold <- fido_18s_s2 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s2 %>% 
  ##MPN: are you throwing out otus that don't appear in the pooled samples?
  filter(rowsum < threshold) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s2_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
##MPN: are the counts of any of these large? Why are they missing?
rows_not_in_taxa_18s <- rownames(fido_18s_s2_final)[!rownames(fido_18s_s2_final) %in% taxa_18s$Hash]

#Check if in the other classification
taxa_18s2=read.csv(here("data/BLAST_taxa_class/zhang_taxa.csv"))  %>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>%
  column_to_rownames("Hash")

rows_only_in_taxa_18s2 <- taxa_18s2 %>% filter(rownames(taxa_18s2) %in% rows_not_in_taxa_18s)


#Add to taxa 18s df


for(row in rows_not_in_taxa_18s) {
  if(row != "other") {
    fido_18s_s2_final["other",] <- fido_18s_s2_final["other",] + fido_18s_s2_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_18s_s2_save <- fido_18s_s2_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s2_save_spp_hash
  
  #Hash
  select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s2_save_hash

#Save
write.csv(fido_18s_s2_save_spp_hash,here("data/fido/fido_18s_s2_ecdf_spp_hash.csv"))
write.csv(fido_18s_s2_save_hash,here("data/fido/fido_18s_s2_ecdf_hash.csv"))




### =============== S3 =========
#s3
fido_18s_s3[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)
###end of added code

fido_18s_s3=fido_18s_s3 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s3$rowsum, thresh_val)

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
above_threshold <- fido_18s_s3 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s3 %>% 
  ##MPN: are you throwing out otus that don't appear in the pooled samples?
  filter(rowsum < threshold) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s3_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
##MPN: are the counts of any of these large? Why are they missing?
rows_not_in_taxa_18s <- rownames(fido_18s_s3_final)[!rownames(fido_18s_s3_final) %in% taxa_18s$Hash]

#Check if in the other classification
taxa_18s2=read.csv(here("data/BLAST_taxa_class/zhang_taxa.csv"))  %>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() %>%
  column_to_rownames("Hash")

rows_only_in_taxa_18s2 <- taxa_18s2 %>% filter(rownames(taxa_18s2) %in% rows_not_in_taxa_18s)


# #Add to taxa 18s df
# 
# 
# for(row in rows_not_in_taxa_18s) {
#   if(row != "other") {
#     fido_18s_s3_final["other",] <- fido_18s_s3_final["other",] + fido_18s_s3_final[row,]
#   }
# }

# Remove the identified rows (excluding 'other')
fido_18s_s3_save <- fido_18s_s3_final %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Genus = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Genus = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # column_to_rownames("spp_hash") %>%
  # select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash)->fido_18s_s3_save_spp_hash

  #Hash
  # column_to_rownames("Hash") %>%
  select(-rowsum,-Phylum,-Class,-Family,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s3_save_hash

#Save
write.csv(fido_18s_s3_save_spp_hash,here("data/fido/fido_18s_s3_ecdf_spp_hash.csv"))
write.csv(fido_18s_s3_save_hash,here("data/fido/fido_18s_s3_ecdf_hash.csv"))
