# ,From these runs, we didn't get data from the 20C or various pools for COI runs, but it worked for 18s
# ,try first with 18s

library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)
here()



###18S
# ,Read in the data
# ,Run 1 (Non pooled data)
asv18s_run1=read.csv(here("data","ASV_table_18s_run1.csv")) %>%
  select(-X) 
# column_to_rownames("Hash")

# ,Run2
asv18s_run2=read.csv(here("data","ASV_table_18s_run2.csv")) %>%
  select(-X)
# column_to_rownames("Hash")



# ,Taxa Tables 
taxa_18s=read.csv(here("data/metazoopruned18s_tax.csv "))
  # column_to_rownames("Hash")


# 2) Merging and manipulation (updated 8/24/2023 to create a new 18S input for fido where
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
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))%>%
  #Add taxa hash
  left_join(taxa_18s, by="Hash")%>%
  #Fill in if spp is missing
  mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  mutate(Family = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>%
  column_to_rownames("Hash") %>%
  select(-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-spp_hash) 


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


#Set ECDF threshold
thresh_val=0.9


#9/26/2023
#ECDF plot for determining criteria 

#90% threshold
#S1
fido_18s_s1=fido_18s_s1 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s1$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_18s_s1 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s1 %>% 
  filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s1_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_18s <- rownames(fido_18s_s1_final)[!rownames(fido_18s_s1_final) %in% taxa_18s$Hash]

for(row in rows_not_in_taxa_18s) {
  if(row != "other") {
    fido_18s_s1_final["other",] <- fido_18s_s1_final["other",] + fido_18s_s1_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_18s_s1_final <- fido_18s_s1_final[!(rownames(fido_18s_s1_final) %in% rows_not_in_taxa_18s & rownames(fido_18s_s1_final) != "other"),] %>%
  dplyr::select(-rowsum) %>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s, by="Hash")%>%
  #Fill in if spp is missing
  mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  mutate(Family = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>%
  column_to_rownames("spp_hash") %>%
  select(-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash) 


#Save
write.csv(fido_18s_s1_final,"data/fido/fido_18s_s1_ecdf_spp_hash.csv")

#Plot ECDF
fido_18s_s1 %>%  ggplot(., aes(rowSums(.))) +
  stat_ecdf(geom = "point")+
  geom_vline(aes(xintercept = threshold), linetype = "dashed", color = "red") +
  labs(title = "ECDF of Row Sums",
       x = "Row Sum",
       y = "ECDF") +
  theme_minimal()


###S2
#s2
fido_18s_s2=fido_18s_s2 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s2$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_18s_s2 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s2 %>% 
  filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s2_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_18s <- rownames(fido_18s_s2_final)[!rownames(fido_18s_s2_final) %in% taxa_18s$Hash]

for(row in rows_not_in_taxa_18s) {
  if(row != "other") {
    fido_18s_s2_final["other",] <- fido_18s_s2_final["other",] + fido_18s_s2_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_18s_s2_final <- fido_18s_s2_final[!(rownames(fido_18s_s2_final) %in% rows_not_in_taxa_18s & rownames(fido_18s_s2_final) != "other"),] %>%
  dplyr::select(-rowsum)%>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s, by="Hash")%>%
  #Fill in if spp is missing
  mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  mutate(Family = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>%
  column_to_rownames("spp_hash") %>%
  select(-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash) 


#Save
write.csv(fido_18s_s2_final,"data/fido/fido_18s_s2_ecdf_spp_hash.csv")

###S3
#s3
fido_18s_s3=fido_18s_s3 %>% mutate(rowsum = rowSums(.[, 10:ncol(.)]))
threshold <- quantile(fido_18s_s3$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_18s_s3 %>% filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any))
below_threshold_sum <- fido_18s_s3 %>% 
  filter(rowsum > threshold & !apply(.[, 1:9] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_18s_s3_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_18s <- rownames(fido_18s_s3_final)[!rownames(fido_18s_s3_final) %in% taxa_18s$Hash]

for(row in rows_not_in_taxa_18s) {
  if(row != "other") {
    fido_18s_s3_final["other",] <- fido_18s_s3_final["other",] + fido_18s_s3_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_18s_s3_final <- fido_18s_s3_final[!(rownames(fido_18s_s3_final) %in% rows_not_in_taxa_18s & rownames(fido_18s_s3_final) != "other"),] %>%
  dplyr::select(-rowsum)%>%
  rownames_to_column("Hash")%>%
  #Add taxa hash
  left_join(taxa_18s, by="Hash")%>%
  #Fill in if spp is missing
  mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  mutate(Family = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>%
  column_to_rownames("spp_hash") %>%
  select(-Species,-Kingdom,-Subphylum,-Subclass,-Superorder,-Hash) 


#Save
write.csv(fido_18s_s3_final,"data/fido/fido_18s_s3_ecdf_spp_hash.csv")


########### PREVIOUS Threshold approach
#Threshold criteria (using >1 count in 30% of the samples)
fido_18s_s1_filt=fido_18s_s1%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)%>%
  filter(rowSums(select(., starts_with("C")) > 1) > 0,
         !rowSums(select(., starts_with("C")) == 0)/ length(select(., starts_with("C"))) > 0.3)
# ,Select inverse and sum to create an "other" category
summed18s_s1=anti_join(fido_18s_s1, fido_18s_s1_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")
# ,Now merge with the pooled samples
all18s_s1=bind_rows(fido_18s_s1_filt,summed18s_s1)
write.csv(all18s_s1,"data/fido/fido_18s_s1.csv")

#S2
fido_18s_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) 

fido_18s_s2_filt=fido_18s_s2%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)%>%
  filter(rowSums(select(., starts_with("C")) > 1) > 0,
         !rowSums(select(., starts_with("C")) == 0)/ length(select(., starts_with("C"))) > 0.3)
# ,Select inverse and sum to create an "other" category
summed18s_s2=anti_join(fido_18s_s2, fido_18s_s2_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")
# ,Now merge with the pooled samples
all18s_s2=bind_rows(fido_18s_s2_filt,summed18s_s2)
write.csv(all18s_s2,"data/fido/fido_18s_s2.csv")

#S3
fido_18s_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) 

fido_18s_s3_filt=fido_18s_s3%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)%>%
  filter(rowSums(select(., starts_with("C")) > 1) > 0,
         !rowSums(select(., starts_with("C")) == 0)/ length(select(., starts_with("C"))) > 0.3)
# ,Select inverse and sum to create an "other" category
summed18s_s3=anti_join(fido_18s_s3, fido_18s_s3_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")

# ,Now merge with the pooled samples
all18s_s3=bind_rows(fido_18s_s3_filt,summed18s_s3)
write.csv(all18s_s3,"data/fido/fido_18s_s3.csv")









### All sizes (not combined)
#Condense down "other" taxa
# ,,Filter samples that don't have a count >1 in more than 30% of samples
all_runs2_filt= all_runs%>%
  mutate(pooled_column = rowSums(. > 1, na.rm = TRUE)) %>%
  filter(rowSums(. > 1, na.rm = TRUE) / ncol(.) >= 0.3) %>%
  select(-pooled_column)

just_pools=all_runs%>%
  dplyr::select(c(contains("C5"),contains("A1"),contains("B3"),contains("All"),contains("Hash")))


# ,Select inverse and sum to create an "other" category
summed18s=anti_join(all_runs, all_runs2_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")



# ,Now merge with the pooled samples
all18s=bind_rows(all_runs2_filt,summed18s)


write.csv(all18s,"data/fido/fido_18s_unaveraged_all.csv")





##All sizes merged in phyloseq
#Read in phyloseq object of merged samples
phy_18s_allsizes=readRDS("data/fido/phy_merged_sizes_18s_fido.rds")

#Extract OTU table, change column names to match, transpose and then add spp hash
allsizes_18s=otu_table(phy_18s_allsizes) %>% t() %>%
  as.data.frame()%>%
  rename_with(~ gsub("-", ".", .))%>%
  rownames_to_column("Hash") %>%
  #Add taxa hash
  left_join(taxa_18s, by="Hash")%>%
  #Fill in if spp is missing
  mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  mutate(Family = if_else(Order=="", Class, Order)) %>%
  
  
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Family, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste(Species,Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>% %>%
  select(-Hash,-Species,-Kingdom,-Subphylum,-Subclass,-Superorder) 



#Join with the all pools
all_runs_justpooled=all_runs%>%
  dplyr::select(c(contains("All"))) %>%
  rownames_to_column("spp_hash")

allsizes_18s_join=merge(all_runs_justpooled,allsizes_18s, by="spp_hash") %>%
  column_to_rownames("spp_hash")

#Filter

allsizes_filt=allsizes_18s_join %>%
  mutate(pooled_column = rowSums(. > 1, na.rm = TRUE)) %>%
  filter(rowSums(. > 1, na.rm = TRUE) / ncol(.) >= 0.3) %>%
  select(-pooled_column)


# ,Select inverse and sum to create an "other" category
summed18s_allsizes=anti_join(allsizes_18s_join, allsizes_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")

# ,Now merge with the pooled samples
all18s_allsizes=bind_rows(allsizes_filt,summed18s_allsizes)
write.csv(all18s_allsizes,"data/fido/fido_18s_allsizes_merged.csv")



###Metadata
#Extract column names and substrings
fido_input=fido_18s_s1_final
column_names <- colnames(fido_input[1:ncol(fido_input)])
substrings <-  substr(column_names, 1,2) %>% as.numeric()
substrings[is.na(substrings)] <- 30
pool <-  substr(column_names, nchar(column_names), nchar(column_names))
sample <-  substr(column_names, 1, nchar(column_names)-2) 



#Create new dataframe for metadata
meta_18s <- data.frame(
  Sample_name = column_names,
  cycle_num =  as.integer(substrings),
  Pool=pool,
  sample_id=sample
) %>%
  mutate(sample_num = ifelse(startsWith(Sample_name, "C"), sample_id, "Calibration"))


write.csv(meta_18s,"data/fido/meta_18s_unaveraged_s1.csv")






