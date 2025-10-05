library (tidyverse)
library (here)
library (lubridate)
library(matrixStats)
library(ggpubr)
library(fido)
library(phyloseq)




################
###COI###

# ,Read in the data
# ,Run 1 (Non pooled data)
asvcoi_run1=read.csv(here("data","fido","ASV_table_coi_run1.csv")) %>% 
  dplyr::select(-X)


# ,Run2
asvcoi_run2=read.csv(here("data","fido","ASV_table_coi_run2.csv"))%>% 
  dplyr::select(-X)


# ,Taxa Tables 
taxa_coi=read.csv(here("data/metazooprunedcoi_tax.csv ")) 



# 2) Merging and manipulation (updated 8/24/2023 to create a new 18S input for fido where
# I don't average technical replicates)
# First need to average technical replicates
# To do this i need to format long
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
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))%>%
  #Add taxa hash
  left_join(taxa_coi, by="Hash")%>%
  # #Fill in if spp is missing
  # mutate(Family = if_else(is.na(Order), Class, Order)) %>%
  # mutate(Family = if_else(Order=="", Class, Order)) %>%
  # 
  # 
  mutate(Family = if_else(is.na(Family), Order, Family)) %>%
  mutate(Family = if_else(Family=="", Order, Family)) %>%
  # 
  mutate(Family = if_else(is.na(Genus),Family, Genus )) %>%
  mutate(Family = if_else(Genus=="",Family, Genus )) %>%
  # 
  mutate(Species = if_else(is.na(Species), Genus, Species))%>%
  mutate(Species = if_else(Species== "", Genus, Species)) %>%
  
  dplyr::select(-Phylum,-Class,-Family,-Genus,-Order) %>%
  mutate(spp_hash=paste0(Species,".",Hash)) %>%
  # mutate(order_hash=paste(Order,Hash)) %>%
  column_to_rownames("Hash") %>%
  select(-Species,-Kingdom,-Subphylum,-Subclass,-Superorder) 

#Replace X
colnames(all_runs) <- gsub("^X", "", colnames(all_runs))


######## ECDF Methof 9/28

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


#Set ECDF threshold
thresh_val=0.9


#9/26/2023
#ECDF plot for determining criteria 

#90% threshold
#S1
fido_coi_s1=fido_coi_s1 %>% mutate(rowsum = rowSums(.[, 7:ncol(.)]))
threshold <- quantile(fido_coi_s1$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_coi_s1 %>% filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any))
below_threshold_sum <- fido_coi_s1 %>% 
  filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s1_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_coi <- rownames(fido_coi_s1_final)[!rownames(fido_coi_s1_final) %in% taxa_coi$Hash]

for(row in rows_not_in_taxa_coi) {
  if(row != "other") {
    fido_coi_s1_final["other",] <- fido_coi_s1_final["other",] + fido_coi_s1_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_coi_s1_final <- fido_coi_s1_final[!(rownames(fido_coi_s1_final) %in% rows_not_in_taxa_coi & rownames(fido_coi_s1_final) != "other"),] %>%
  dplyr::select(-rowsum)


#Save
write.csv(fido_coi_s1_final,"data/fido/fido_coi_s1_ecdf.csv")

#Plot ECDF
fido_coi_s1 %>%  ggplot(., aes(rowSums(.))) +
  stat_ecdf(geom = "point")+
  geom_vline(aes(xintercept = threshold), linetype = "dashed", color = "red") +
  labs(title = "ECDF of Row Sums",
       x = "Row Sum",
       y = "ECDF") +
  theme_minimal()


###S2
#s2
fido_coi_s2=fido_coi_s2 %>% mutate(rowsum = rowSums(.[, 7:ncol(.)]))
threshold <- quantile(fido_coi_s2$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_coi_s2 %>% filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any))
below_threshold_sum <- fido_coi_s2 %>% 
  filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s2_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_coi <- rownames(fido_coi_s2_final)[!rownames(fido_coi_s2_final) %in% taxa_coi$Hash]

for(row in rows_not_in_taxa_coi) {
  if(row != "other") {
    fido_coi_s2_final["other",] <- fido_coi_s2_final["other",] + fido_coi_s2_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_coi_s2_final <- fido_coi_s2_final[!(rownames(fido_coi_s2_final) %in% rows_not_in_taxa_coi & rownames(fido_coi_s2_final) != "other"),] %>%
  dplyr::select(-rowsum)


#Save
write.csv(fido_coi_s2_final,"data/fido/fido_coi_s2_ecdf.csv")

###S3
#s3
fido_coi_s3=fido_coi_s3 %>% mutate(rowsum = rowSums(.[, 7:ncol(.)]))
threshold <- quantile(fido_coi_s3$rowsum, thresh_val)

# Separate rows based on threshold
above_threshold <- fido_coi_s3 %>% filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any))
below_threshold_sum <- fido_coi_s3 %>% 
  filter(rowsum > threshold & !apply(.[, 1:6] == 0, 1, any)) %>%
  summarise_all(sum) %>% 
  mutate(rowname = "other") %>%
  column_to_rownames("rowname")


# Combine data
fido_coi_s3_final <- bind_rows(above_threshold, below_threshold_sum)

#Remove mismatched rows and merge
rows_not_in_taxa_coi <- rownames(fido_coi_s3_final)[!rownames(fido_coi_s3_final) %in% taxa_coi$Hash]

for(row in rows_not_in_taxa_coi) {
  if(row != "other") {
    fido_coi_s3_final["other",] <- fido_coi_s3_final["other",] + fido_coi_s3_final[row,]
  }
}

# Remove the identified rows (excluding 'other')
fido_coi_s3_final <- fido_coi_s3_final[!(rownames(fido_coi_s3_final) %in% rows_not_in_taxa_coi & rownames(fido_coi_s3_final) != "other"),] %>%
  dplyr::select(-rowsum)


#Save
write.csv(fido_coi_s3_final,"data/fido/fido_coi_s3_ecdf.csv")




































#9/1/2023
#Separate out by size
#S1
fido_coi_s1=all_runs%>%
  dplyr::select(c(contains("All"),contains("S1"))) 

fido_coi_s1_filt=fido_coi_s1%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)
# ,Select inverse and sum to create an "other" category
summedcoi_s1=anti_join(fido_coi_s1, fido_coi_s1_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")
# ,Now merge with the pooled samples
allcoi_s1=bind_rows(fido_coi_s1_filt,summedcoi_s1)
write.csv(allcoi_s1,"data/fido/fido_coi_s1.csv")

#S2
fido_coi_s2=all_runs%>%
  dplyr::select(c(contains("All"),contains("S2"))) 

fido_coi_s2_filt=fido_coi_s2%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)
# ,Select inverse and sum to create an "other" category
summedcoi_s2=anti_join(fido_coi_s2, fido_coi_s2_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")
# ,Now merge with the pooled samples
allcoi_s2=bind_rows(fido_coi_s2_filt,summedcoi_s2)
write.csv(allcoi_s2,"data/fido/fido_coi_s2.csv")

#S3
fido_coi_s3=all_runs%>%
  dplyr::select(c(contains("All"),contains("S3"))) 

fido_coi_s3_filt=fido_coi_s3%>%
  filter(rowSums(select(., starts_with("2")) > 1) > 0,
         !rowSums(select(., starts_with("2")) == 0) > 0)
# ,Select inverse and sum to create an "other" category
summedcoi_s3=anti_join(fido_coi_s3, fido_coi_s3_filt)%>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  summarize(across(everything(), sum)) %>%
  mutate(Hash="Other") %>%
  column_to_rownames("Hash")

# ,Now merge with the pooled samples
allcoi_s3=bind_rows(fido_coi_s3_filt,summedcoi_s3)
write.csv(allcoi_s3,"data/fido/fido_coi_s3.csv")


###Metadata
#Extract column names and substrings
fido_input=allcoi_s3
column_names <- colnames(fido_input[1:ncol(fido_input)])
substrings <-  substr(column_names, 1,2) %>% as.numeric()
substrings[is.na(substrings)] <- 30
pool <-  substr(column_names, nchar(column_names), nchar(column_names))
sample <-  substr(column_names, 1, nchar(column_names)-2) 



#Create new dataframe for metadata
meta_coi <- data.frame(
  Sample_name = column_names,
  cycle_num =  as.integer(substrings),
  Pool=pool,
  sample_id=sample
) %>%
  mutate(sample_num = ifelse(startsWith(Sample_name, "C"), sample_id, "Calibration"))


# write.csv(meta_coi,"data/fido/meta_coi_unaveraged_s1.csv")
# write.csv(meta_coi,"data/fido/meta_coi_unaveraged_s2.csv")
write.csv(meta_coi,"data/fido/meta_coi_unaveraged_s3.csv")

