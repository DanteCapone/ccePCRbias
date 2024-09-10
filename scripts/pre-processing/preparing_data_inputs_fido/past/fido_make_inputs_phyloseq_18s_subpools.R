# Make subpool inputs -----------------------------------------------------


# From these runs, we didn't get data from the 20C or various pools for 18s runs, but it worked for 18s
# try first with 18s




library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)
here()



# Data description:
#   We ran 3 pools of samples at 3 different PCR Cycles (20,24,28 cycles) with 3 replicates per cycle.


# Make strings with Pool samples ------------------------------------------


# A1-A3:
a1_a3=c(
  "C1.T7.H9_S1",
  "C1.T7.H9_S2",
  "C1.T7.H9_S3",
  "C1.T8.H10_S1",
  "C1.T8.H10_S2",
  "C1.T8.H10_S3",
  "C2.T8.H18_S1",
  "C2.T8.H18_S2",
  "C2.T8.H18_S3",
  "C2.T9.H19_S1",
  "C2.T9.H19_S2",
  "C2.T9.H19_S3",
  "C3.T6.H25_S1",
  "C3.T6.H25_S2",
  "C3.T6.H25_S3",
  "C3.T7.H26_S1",
  "C3.T7.H26_S2",
  "A1.A3"
)

# B3-B5:
b3_b5=c("C3.T7.H26_S3",
        "CT1.T1.H28_S1",
        "CT1.T1.H28_S2",
        "CT1.T1.H28_S3",
        "CT1.T2.H29_S1",
        "CT1.T2.H29_S2",
        "CT1.T2.H29_S3",
        "CT1.T3.H30_S1",
        "CT1.T3.H30_S2",
        "CT1.T3.H30_S3",
        "CT1.T4.H31_S1",
        "CT1.T4.H31_S2",
        "CT1.T4.H31_S3",
        "CT1.T5.H32_S1",
        "CT1.T5.H32_S2",
        "CT1.T5.H32_S3",
        "CT1.T6.H33_S1",
        "B3.B5"
)
# 
# 
# C5-C7:
c5_c7=c("CT1.T6.H33_S2",
        "CT1.T6.H33_S3",
        "CT1.T7.H34_S1",
        "CT1.T7.H34_S2",
        "CT1.T7.H34_S3",
        "CT1.T8.H35_S1",
        "CT1.T8.H35_S2",
        "CT1.T8.H35_S3",
        "CT2.T1.H36_S1",
        "CT2.T1.H36_S2",
        "CT2.T1.H36_S3",
        "CT2.T4.H39_S1",
        "CT2.T4.H39_S2",
        "CT2.T4.H39_S3",
        "CT2.T8.H43_S1",
        "CT2.T8.H43_S2",
        "CT2.T8.H43_S3",
        "C5.C7"
)


# ,Read in the data
# ,Run 1 (Non pooled data)
asv18s_run1=read.csv(here("data","fido","ASV_table_18s_run1.csv")) %>% 
  dplyr::select(-X)


# ,Run2
asv18s_run2=read.csv(here("data","fido","ASV_table_18s_run2.csv"))%>% 
  dplyr::select(-X)


# ,Taxa Tables 
taxa_18s=read.csv(here("data/phyloseq_bio_data/18s/metazoopruned18s_tax.csv")) %>%
  column_to_rownames("Hash")



# 2) Merging and manipulation (updated 8/24/2023 to create a new 18s input for fido where
# I don't average technical replicates)
# First need to average technical replicates
# To do this i need to format long
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

##MPN: Another common method is to amalgamate to the family level and add anything that doesn't fit to an "other" category. This is very common in microbiome.

#Replace X
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))



# Pool 1: A1-A3 -----------------------------------------------------------

#Separate out by size and pool
#A1-A3
fido_18s_s1_a1_a3=all_runs%>%
  select(contains(a1_a3)) %>% 
  select(c(contains("A1.A3"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s2_a1_a3=all_runs%>%
  select(contains(a1_a3)) %>% 
  select(c(contains("A1.A3"),contains("S2"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s3_a1_a3=all_runs%>%
  select(contains(a1_a3)) %>% 
  select(c(contains("A1.A3"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)




##Phyloseq filtering
fido_18s_s1_otu_a1_a3=fido_18s_s1_a1_a3 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu_a1_a3=fido_18s_s2_a1_a3 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu_a1_a3=fido_18s_s3_a1_a3 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
tax18s_s1_a1_a3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1_otu_a1_a3))
tax18s_s1_a1_a3 =  tax_table(as.matrix(tax18s_s1_a1_a3))

#S2
tax18s_s2_a1_a3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu_a1_a3))
tax18s_s2_a1_a3 =  tax_table(as.matrix(tax18s_s2_a1_a3))
#S3
tax18s_s3_a1_a3 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu_a1_a3))
tax18s_s3_a1_a3 =  tax_table(as.matrix(tax18s_s3_a1_a3))

#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy_a1_a3=phyloseq(fido_18s_s1_otu_a1_a3,tax18s_s1_a1_a3)
fido_18s_s2_phy_a1_a3=phyloseq(fido_18s_s2_otu_a1_a3,tax18s_s2_a1_a3)
fido_18s_s3_phy_a1_a3=phyloseq(fido_18s_s3_otu_a1_a3,tax18s_s3_a1_a3)



#PHYLOSEQ
#Agglomerate at the family level

#S1
fido_18s_s1_phy_a1_a3=phyloseq(fido_18s_s1_otu_a1_a3,tax18s_s1_a1_a3, metadata)
fido_18s_s1_family_a1_a3=tax_glom(fido_18s_s1_phy_a1_a3, taxrank = "Family")

#Make inputs for filtering
fido_18s_s1_family_otu_a1_a3=otu_table(fido_18s_s1_family_a1_a3) %>% as.data.frame()
fido_18s_s1_family_taxa_a1_a3=tax_table(fido_18s_s1_family_a1_a3) %>% as.data.frame()

#==S2
fido_18s_s2_phy_a1_a3=phyloseq(fido_18s_s2_otu_a1_a3,tax18s_s2_a1_a3, metadata)
fido_18s_s2_family_a1_a3=tax_glom(fido_18s_s2_phy_a1_a3, taxrank = "Family")

#Make inputs for filtering
fido_18s_s2_family_otu_a1_a3=otu_table(fido_18s_s2_family_a1_a3) %>% as.data.frame()
fido_18s_s2_family_taxa_a1_a3=tax_table(fido_18s_s2_family_a1_a3) %>% as.data.frame()

#==s3
fido_18s_s3_phy_a1_a3=phyloseq(fido_18s_s3_otu_a1_a3,tax18s_s3_a1_a3, metadata)
fido_18s_s3_family_a1_a3=tax_glom(fido_18s_s3_phy_a1_a3, taxrank = "Family")

#Make inputs for filtering
fido_18s_s3_family_otu_a1_a3=otu_table(fido_18s_s3_family_a1_a3) %>% as.data.frame()
fido_18s_s3_family_taxa_a1_a3=tax_table(fido_18s_s3_family_a1_a3) %>% as.data.frame()


#Save aglomerated family taxa file
tax18s_family_a1_a3=rbind(fido_18s_s1_family_taxa_a1_a3,fido_18s_s2_family_taxa_a1_a3,fido_18s_s3_family_taxa_a1_a3) %>%
  unique()

# Create a new row with "other" in all columns
other <- data.frame(lapply(tax18s_family_a1_a3, function(x) "other"))

# Bind the new row to the original dataframe
tax18s_family_a1_a3 <- bind_rows(tax18s_family_a1_a3,other)

write.csv(tax18s_family_a1_a3,here("data/phyloseq_bio_data/18s/sub_pools/fido_18s_family_tax_table_a1_a3.csv"))




#Visualize ECDF
fido_18s_s1_family_otu_a1_a3[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s1_family_otu_a1_a3=fido_18s_s1_family_otu_a1_a3 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_a1_a3 <- fido_18s_s1_family_otu_a1_a3 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s1_family_otu_a1_a3 %>%
  anti_join(fido_18s_s1_family_otu_a1_a3 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s1_final_a1_a3 <- bind_rows(fido_taxa_sel_a1_a3, other)


# Remove the identified rows (excluding 'other')
fido_18s_s1_final_a1_a3 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s1_save_taxa_phy_a1_a3

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s1_save_hash_phy

#Save
write.csv(fido_18s_s1_save_taxa_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_a1_a3.csv"))
# write.csv(fido_18s_s1_save_hash_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_hash_phy_a1_a3.csv"))



### ========== S2 ===============


#Visualize ECDF
fido_18s_s2_family_otu_a1_a3[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s2_family_otu_a1_a3=fido_18s_s2_family_otu_a1_a3 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_a1_a3 <- fido_18s_s2_family_otu_a1_a3 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s2_family_otu_a1_a3 %>%
  anti_join(fido_18s_s2_family_otu_a1_a3 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s2_final_a1_a3 <- bind_rows(fido_taxa_sel_a1_a3, other)


# Remove the identified rows (excluding 'other')
fido_18s_s2_final_a1_a3 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s2_save_taxa_phy_a1_a3

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s2_save_hash_phy

#Save
write.csv(fido_18s_s2_save_taxa_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_a1_a3.csv"))
# write.csv(fido_18s_s2_save_hash_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_hash_phy_a1_a3.csv"))





### =============== S3 =========

#Visualize ECDF
fido_18s_s3_family_otu_a1_a3[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s3_family_otu_a1_a3=fido_18s_s3_family_otu_a1_a3 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_a1_a3 <- fido_18s_s3_family_otu_a1_a3 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s3_family_otu_a1_a3 %>%
  anti_join(fido_18s_s3_family_otu_a1_a3 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s3_final_a1_a3 <- bind_rows(fido_taxa_sel_a1_a3, other)


# Remove the identified rows (excluding 'other')
fido_18s_s3_final_a1_a3 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s3_save_taxa_phy_a1_a3

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s3_save_hash_phy

#Save
write.csv(fido_18s_s3_save_taxa_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_a1_a3.csv"))
# write.csv(fido_18s_s3_save_hash_phy_a1_a3,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_hash_phy_a1_a3.csv"))





# Pool 2: B3-B5  ----------------------------------------------------------
#B3-B5
fido_18s_s1_b3_b5=all_runs%>%
  select(contains(b3_b5)) %>% 
  select(c(contains("B3.B5"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s2_b3_b5=all_runs%>%
  select(contains(b3_b5)) %>% 
  select(c(contains("B3.B5"),contains("S2"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s3_b3_b5=all_runs%>%
  select(contains(b3_b5)) %>% 
  select(c(contains("B3.B5"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)


##Phyloseq filtering
fido_18s_s1_otu_b3_b5=fido_18s_s1_b3_b5 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu_b3_b5=fido_18s_s2_b3_b5 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu_b3_b5=fido_18s_s3_b3_b5 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
tax18s_s1_b3_b5 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1_otu_b3_b5))
tax18s_s1_b3_b5 =  tax_table(as.matrix(tax18s_s1_b3_b5))

#S2
tax18s_s2_b3_b5 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu_b3_b5))
tax18s_s2_b3_b5 =  tax_table(as.matrix(tax18s_s2_b3_b5))
#S3
tax18s_s3_b3_b5 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu_b3_b5))
tax18s_s3_b3_b5 =  tax_table(as.matrix(tax18s_s3_b3_b5))

#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy_b3_b5=phyloseq(fido_18s_s1_otu_b3_b5,tax18s_s1_b3_b5)
fido_18s_s2_phy_b3_b5=phyloseq(fido_18s_s2_otu_b3_b5,tax18s_s2_b3_b5)
fido_18s_s3_phy_b3_b5=phyloseq(fido_18s_s3_otu_b3_b5,tax18s_s3_b3_b5)



#PHYLOSEQ
#Agglomerate at the family level
fido_18s_s1_phy_b3_b5=phyloseq(fido_18s_s1_otu_b3_b5,tax18s_s1_b3_b5, metadata)
fido_18s_s1_family_b3_b5=tax_glom(fido_18s_s1_phy_b3_b5, taxrank = "Family")



#Make inputs for filtering
fido_18s_s1_family_otu_b3_b5=otu_table(fido_18s_s1_family_b3_b5) %>% as.data.frame()
fido_18s_s1_family_taxa_b3_b5=tax_table(fido_18s_s1_family_b3_b5) %>% as.data.frame()

#==S2
fido_18s_s2_phy_b3_b5=phyloseq(fido_18s_s2_otu_b3_b5,tax18s_s2_b3_b5, metadata)
fido_18s_s2_family_b3_b5=tax_glom(fido_18s_s2_phy_b3_b5, taxrank = "Family")

#Make inputs for filtering
fido_18s_s2_family_otu_b3_b5=otu_table(fido_18s_s2_family_b3_b5) %>% as.data.frame()
fido_18s_s2_family_taxa_b3_b5=tax_table(fido_18s_s2_family_b3_b5) %>% as.data.frame()

#==s3
fido_18s_s3_phy_b3_b5=phyloseq(fido_18s_s3_otu_b3_b5,tax18s_s3_b3_b5, metadata)
fido_18s_s3_family_b3_b5=tax_glom(fido_18s_s3_phy_b3_b5, taxrank = "Family")

#Make inputs for filtering
fido_18s_s3_family_otu_b3_b5=otu_table(fido_18s_s3_family_b3_b5) %>% as.data.frame()
fido_18s_s3_family_taxa_b3_b5=tax_table(fido_18s_s3_family_b3_b5) %>% as.data.frame()


#Save aglomerated family taxa file
tax18s_family_b3_b5=rbind(fido_18s_s1_family_taxa_b3_b5,fido_18s_s2_family_taxa_b3_b5,fido_18s_s3_family_taxa_b3_b5) %>%
  unique()

# Create a new row with "other" in all columns
other <- data.frame(lapply(tax18s_family_b3_b5, function(x) "other"))

# Bind the new row to the original dataframe
tax18s_family_b3_b5 <- bind_rows(tax18s_family_b3_b5,other)

write.csv(tax18s_family_b3_b5,here("data/phyloseq_bio_data/18s/sub_pools/fido_18s_family_tax_table_b3_b5.csv"))



# S1 ----------------------------------------------------------------------
#Visualize ECDF
fido_18s_s1_family_otu_b3_b5[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)



#Add rowsums
fido_18s_s1_family_otu_b3_b5=fido_18s_s1_family_otu_b3_b5 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_b3_b5 <- fido_18s_s1_family_otu_b3_b5 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s1_family_otu_b3_b5 %>%
  anti_join(fido_18s_s1_family_otu_b3_b5 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s1_final_b3_b5 <- bind_rows(fido_taxa_sel_b3_b5, other)


# Remove the identified rows (excluding 'other')
fido_18s_s1_final_b3_b5 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s1_save_taxa_phy_b3_b5

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s1_save_hash_phy

#Save
write.csv(fido_18s_s1_save_taxa_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_b3_b5.csv"))
# write.csv(fido_18s_s1_save_hash_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_hash_phy_b3_b5.csv"))



### ========== S2 ===============


#Visualize ECDF
fido_18s_s2_family_otu_b3_b5[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s2_family_otu_b3_b5=fido_18s_s2_family_otu_b3_b5 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_b3_b5 <- fido_18s_s2_family_otu_b3_b5 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s2_family_otu_b3_b5 %>%
  anti_join(fido_18s_s2_family_otu_b3_b5 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s2_final_b3_b5 <- bind_rows(fido_taxa_sel_b3_b5, other)


# Remove the identified rows (excluding 'other')
fido_18s_s2_final_b3_b5 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s2_save_taxa_phy_b3_b5

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s2_save_hash_phy

#Save
write.csv(fido_18s_s2_save_taxa_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_b3_b5.csv"))
# write.csv(fido_18s_s2_save_hash_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_hash_phy_b3_b5.csv"))





### =============== S3 =========

#Visualize ECDF
fido_18s_s3_family_otu_b3_b5[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s3_family_otu_b3_b5=fido_18s_s3_family_otu_b3_b5 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_b3_b5 <- fido_18s_s3_family_otu_b3_b5 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s3_family_otu_b3_b5 %>%
  anti_join(fido_18s_s3_family_otu_b3_b5 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s3_final_b3_b5 <- bind_rows(fido_taxa_sel_b3_b5, other)


# Remove the identified rows (excluding 'other')
fido_18s_s3_final_b3_b5 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s3_save_taxa_phy_b3_b5

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s3_save_hash_phy

#Save
write.csv(fido_18s_s3_save_taxa_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_b3_b5.csv"))
# write.csv(fido_18s_s3_save_hash_phy_b3_b5,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_hash_phy_b3_b5.csv"))







# Pool 3: C5-C7 -----------------------------------------------------------

#C5-C7
fido_18s_s1_c5_c7=all_runs%>%
  select(contains(c5_c7)) %>% 
  select(c(contains("C5.C7"),contains("S1"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s2_c5_c7=all_runs%>%
  select(contains(c5_c7)) %>% 
  select(c(contains("C5.C7"),contains("S2"))) %>% 
  filter(rowSums(.) != 0) 

fido_18s_s3_c5_c7=all_runs%>%
  select(contains(c5_c7)) %>% 
  select(c(contains("C5.C7"),contains("S3"))) %>% 
  filter(rowSums(.) != 0)




##Phyloseq filtering
fido_18s_s1_otu_c5_c7=fido_18s_s1_c5_c7 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s2_otu_c5_c7=fido_18s_s2_c5_c7 %>% 
  otu_table(taxa_are_rows = TRUE)

fido_18s_s3_otu_c5_c7=fido_18s_s3_c5_c7 %>% 
  otu_table(taxa_are_rows = TRUE)

#taxa table
tax18s_s1_c5_c7 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s1_otu_c5_c7))
tax18s_s1_c5_c7 =  tax_table(as.matrix(tax18s_s1_c5_c7))

#S2
tax18s_s2_c5_c7 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s2_otu_c5_c7))
tax18s_s2_c5_c7 =  tax_table(as.matrix(tax18s_s2_c5_c7))
#S3
tax18s_s3_c5_c7 = taxa_18s %>% filter(rownames(taxa_18s) %in% rownames(fido_18s_s3_otu_c5_c7))
tax18s_s3_c5_c7 =  tax_table(as.matrix(tax18s_s3_c5_c7))

#Metadata
meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(-c(Sizefractionmm,offshore_onshore,clust_group,PC1,cycle, max_size)) %>%
  sample_data(.)

fido_18s_s1_phy_c5_c7=phyloseq(fido_18s_s1_otu_c5_c7,tax18s_s1_c5_c7)
fido_18s_s2_phy_c5_c7=phyloseq(fido_18s_s2_otu_c5_c7,tax18s_s2_c5_c7)
fido_18s_s3_phy_c5_c7=phyloseq(fido_18s_s3_otu_c5_c7,tax18s_s3_c5_c7)



#PHYLOSEQ
#Agglomerate at the family level

#S1
fido_18s_s1_phy_c5_c7=phyloseq(fido_18s_s1_otu_c5_c7,tax18s_s1_c5_c7, metadata)
fido_18s_s1_family_c5_c7=tax_glom(fido_18s_s1_phy_c5_c7, taxrank = "Family")

#Make inputs for filtering
fido_18s_s1_family_otu_c5_c7=otu_table(fido_18s_s1_family_c5_c7) %>% as.data.frame()
fido_18s_s1_family_taxa_c5_c7=tax_table(fido_18s_s1_family_c5_c7) %>% as.data.frame()

#==S2
fido_18s_s2_phy_c5_c7=phyloseq(fido_18s_s2_otu_c5_c7,tax18s_s2_c5_c7, metadata)
fido_18s_s2_family_c5_c7=tax_glom(fido_18s_s2_phy_c5_c7, taxrank = "Family")

#Make inputs for filtering
fido_18s_s2_family_otu_c5_c7=otu_table(fido_18s_s2_family_c5_c7) %>% as.data.frame()
fido_18s_s2_family_taxa_c5_c7=tax_table(fido_18s_s2_family_c5_c7) %>% as.data.frame()

#==s3
fido_18s_s3_phy_c5_c7=phyloseq(fido_18s_s3_otu_c5_c7,tax18s_s3_c5_c7, metadata)
fido_18s_s3_family_c5_c7=tax_glom(fido_18s_s3_phy_c5_c7, taxrank = "Family")

#Make inputs for filtering
fido_18s_s3_family_otu_c5_c7=otu_table(fido_18s_s3_family_c5_c7) %>% as.data.frame()
fido_18s_s3_family_taxa_c5_c7=tax_table(fido_18s_s3_family_c5_c7) %>% as.data.frame()


#Save aglomerated family taxa file
tax18s_family_c5_c7=rbind(fido_18s_s1_family_taxa_c5_c7,fido_18s_s2_family_taxa_c5_c7,fido_18s_s3_family_taxa_c5_c7) %>%
  unique()

# Create a new row with "other" in all columns
other <- data.frame(lapply(tax18s_family_c5_c7, function(x) "other"))

# Bind the new row to the original dataframe
tax18s_family_c5_c7 <- bind_rows(tax18s_family_c5_c7,other)

write.csv(tax18s_family_c5_c7,here("data/phyloseq_bio_data/18s/sub_pools/fido_18s_family_tax_table_c5_c7.csv"))




#Visualize ECDF
fido_18s_s1_family_otu_c5_c7[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s1_family_otu_c5_c7=fido_18s_s1_family_otu_c5_c7 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_c5_c7 <- fido_18s_s1_family_otu_c5_c7 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s1_family_otu_c5_c7 %>%
  anti_join(fido_18s_s1_family_otu_c5_c7 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s1_final_c5_c7 <- bind_rows(fido_taxa_sel_c5_c7, other)


# Remove the identified rows (excluding 'other')
fido_18s_s1_final_c5_c7 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s1_save_taxa_phy_c5_c7

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s1_save_hash_phy

#Save
write.csv(fido_18s_s1_save_taxa_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_c5_c7.csv"))
# write.csv(fido_18s_s1_save_hash_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_hash_phy_c5_c7.csv"))



### ========== S2 ===============


#Visualize ECDF
fido_18s_s2_family_otu_c5_c7[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s2_family_otu_c5_c7=fido_18s_s2_family_otu_c5_c7 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_c5_c7 <- fido_18s_s2_family_otu_c5_c7 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s2_family_otu_c5_c7 %>%
  anti_join(fido_18s_s2_family_otu_c5_c7 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s2_final_c5_c7 <- bind_rows(fido_taxa_sel_c5_c7, other)


# Remove the identified rows (excluding 'other')
fido_18s_s2_final_c5_c7 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s2_save_taxa_phy_c5_c7

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s2_save_hash_phy

#Save
write.csv(fido_18s_s2_save_taxa_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_c5_c7.csv"))
# write.csv(fido_18s_s2_save_hash_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_hash_phy_c5_c7.csv"))





### =============== S3 =========

#Visualize ECDF
fido_18s_s3_family_otu_c5_c7[,-c(1:10)] %>% rowSums() %>% ecdf() %>% plot() %>% abline(v=1637)

#Add rowsums
fido_18s_s3_family_otu_c5_c7=fido_18s_s3_family_otu_c5_c7 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))

# Separate rows based on threshold
##MPN: Why do you think some of the hashes are appearing quite high in some samples but not in any of the pooled samples?
fido_taxa_sel_c5_c7 <- fido_18s_s3_family_otu_c5_c7 %>% filter(rowSums(select(., 1:9) == 0) <= 2)
other <- fido_18s_s3_family_otu_c5_c7 %>%
  anti_join(fido_18s_s3_family_otu_c5_c7 %>%
              filter(rowSums(select(., 1:9) == 0) <= 2))%>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s3_final_c5_c7 <- bind_rows(fido_taxa_sel_c5_c7, other)


# Remove the identified rows (excluding 'other')
fido_18s_s3_final_c5_c7 %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s %>% rownames_to_column("Hash"), by="Hash")%>%
  #Fill in if spp is missing
  mutate(Order = if_else(is.na(Order), Class, Order)) %>%
  mutate(Order = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Family),Family, Family )) %>%
  mutate(Family = if_else(Family=="",Family, Family )) %>%
  
  mutate(Species = if_else(is.na(Species), Family, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  #Uncomment to save spp_hash
  select(-rowsum,-Phylum,-Class,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash,-spp_hash)->fido_18s_s3_save_taxa_phy_c5_c7

#Hash
# select(-rowsum,-Phylum,-Class,-Genus,-Genus,-Order,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash)->fido_18s_s3_save_hash_phy

#Save
write.csv(fido_18s_s3_save_taxa_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_c5_c7.csv"))
# write.csv(fido_18s_s3_save_hash_phy_c5_c7,here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_hash_phy_c5_c7.csv"))

