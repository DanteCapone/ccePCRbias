
#Script Using Fido to Predict the DNA Proportions in each sample at PCR Cycle 0

# Libraries, data, etc ----------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggpubr)
library(fido)
library(stringr)
library(here)
library(gridExtra)
here()

# Back-extrapolation target cycle (see scripts/pcr_correction_pipeline/config.R)
if (!exists("target_cycle")) target_cycle <- 10



###Load in the filtered data for the 18s primer using long format species and hash name so I ca identify taxa
##First Size 1

fido_input_filt=read.csv(here("data/fido/asv_level/fido_18s_s1_asv_nofilt.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>% 
  group_by(Hash) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Hash") %>%
  dplyr::select(c(contains("All"),contains("A1"),contains("C5"),contains("S1")))

#Metadata
meta_18s=read.csv(here("data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt)) %>%
  mutate(Run = case_when(
    grepl("Pool", Sample_name) ~ "2",  # Prioritize "Pool" samples for Run 2
    grepl("\\.1|\\.2", Sample_name) ~ "1",  # Otherwise, assign Run 1 if it ends in .1 or .2
    grepl("\\.3", Sample_name) ~ "2",  # Assign Run 2 if it ends in .3
    TRUE ~ NA_character_  # Default case
  )) %>%
  mutate(Run = factor(Run, levels = c("1", "2")))  # Explicitly define levels
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

#Add in vars for onshore offshore
# read in metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  mutate(sample_id=Sample_ID_dot) %>%
  select(sample_id,day_night_0_1, offshore_onshore) %>% 
  distinct(.) 

meta_18s <- meta_18s %>% 
  left_join(metadata, by="sample_id")%>%
  mutate(
    offshore = ifelse(offshore_onshore == "offshore", 1, 0),
    onshore = ifelse(offshore_onshore == "onshore", 1, 0),
    day = ifelse(day_night_0_1 == 0, 1, 0),
    night = ifelse(day_night_0_1 == 1, 1, 0)
  ) %>%
  replace_na(list(offshore = 0, onshore = 0, day = 0, night = 0)) %>% # Replace NA values with 0
  select(-offshore_onshore, -day_night_0_1)  # Drop the original categorical columns


##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),] 

meta_18s=meta_18s%>%
  filter(!is.na(Sample_name))

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+sample_num+Run+offshore+onshore+day+night  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 



#Fit pibble model 
#Loop thru values for Gamma
gamma <- c(1,2,3,5,8,10,15,20,50,100,200,300,400,500,700,1000)
logML <- rep(NA, length(gamma))
for(i in 1:length(gamma)){
  fit <- pibble(Y_s1, X, Gamma = gamma[i]*diag(nrow(X)), n_samples=5000)
  logML[i] <- fit$logMarginalLikelihood
  print(i)
}

plot(gamma, logML, type = "l")
points(gamma, logML)
#20 seems good based on LML
gamma=20

#Specify the remaining priors with default values
upsilon <- nrow(Y_s1)+3 
Omega <- diag(nrow(Y_s1))
G <- cbind(diag(nrow(Y_s1)-1), -1)
Xi <- (upsilon-nrow(Y_s1))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s1)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = gamma*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)
#summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
##Looks ok, centered at zero
##end of added code

##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s1 <- to_clr(fit)

#Convert to Proportions
fit_prop_1 <- to_proportions(fit_s1)





############


X.tmp.s1 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s1) <- rownames(X)
X.tmp.s1["cycle_num", ] <- target_cycle  # evaluate fitted line at target_cycle


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s1 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num","Run2",
                        "offshore","onshore","day","night"))%>% as.data.frame()->samples_to_loop 

#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s1 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s1[s,] <-1
  
  
  #
  predicted_s1 <- predict(fit_prop_1, newdata=X.tmp.s1, summary=TRUE) %>% 
    mutate(cycle_num = c(target_cycle)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s1 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>% 
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", "")))) %>% 
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("0.2-0.5mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s1,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  ##MPN: NOt exactly sure what you are trying to show with this plot
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  
  
  
  
  final_data_s1 <- bind_rows(final_data_s1, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s1[s,] <-0
  
}

#BEEP to notify when finished lol
beepr::beep(12)

#Save final data
current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s1,here(paste0("data/predicted_og/asv/predicted_og_18s_",current_date,"_s1_phy_all_and_subpools_nofilt_asv.csv")))





# Repeat for other sizes

############First 0.5-1############
#Phyloseq Filtered

fido_input_filt=read.csv(here("data/fido/asv_level/fido_18s_s2_asv_nofilt.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  group_by(Hash) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Hash") %>%
  dplyr::select(c(contains("All"),contains("A1"),contains("C5"),contains("s2")))

#Metadata
meta_18s=read.csv(here("data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt)) %>%
  mutate(Run = case_when(
    grepl("Pool", Sample_name) ~ "2",  # Prioritize "Pool" samples for Run 2
    grepl("\\.1|\\.2", Sample_name) ~ "1",  # Otherwise, assign Run 1 if it ends in .1 or .2
    grepl("\\.3", Sample_name) ~ "2",  # Assign Run 2 if it ends in .3
    grepl("\\.4", Sample_name) ~ "2",  # Assign Run 2 if it ends in .3
    TRUE ~ NA_character_  # Default case
  )) %>%
  mutate(Run = factor(Run, levels = c("1", "2")))  # Explicitly define levels
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

#Add in vars for onshore offshore
# read in metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  mutate(sample_id=Sample_ID_dot) %>%
  select(sample_id,day_night_0_1, offshore_onshore) %>% 
  distinct(.) 

meta_18s <- meta_18s %>% 
  left_join(metadata, by="sample_id")%>%
  mutate(
    offshore = ifelse(offshore_onshore == "offshore", 1, 0),
    onshore = ifelse(offshore_onshore == "onshore", 1, 0),
    day = ifelse(day_night_0_1 == 0, 1, 0),
    night = ifelse(day_night_0_1 == 1, 1, 0)
  ) %>%
  replace_na(list(offshore = 0, onshore = 0, day = 0, night = 0)) %>% # Replace NA values with 0
  select(-offshore_onshore, -day_night_0_1)  # Drop the original categorical columns


##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),] 

meta_18s=meta_18s%>%
  filter(!is.na(Sample_name))

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+sample_num+Run+offshore+onshore+day+night  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 






#Fit pibble model 
gamma=20

#Specify the remaining priors with default values
upsilon <- nrow(Y_s2)+3 
Omega <- diag(nrow(Y_s2))
G <- cbind(diag(nrow(Y_s2)-1), -1)
Xi <- (upsilon-nrow(Y_s2))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s2)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = gamma*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)


##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s2 <- to_clr(fit)

#Convert to Proportions
fit_prop_2 <- to_proportions(fit_s2)





############


X.tmp.s2 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s2) <- rownames(X)
X.tmp.s2["cycle_num", ] <- target_cycle  # evaluate fitted line at target_cycle


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s2 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num","Run2",
                        "offshore","onshore","day","night"))%>% as.data.frame()->samples_to_loop 
#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s2 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s2[s,] <-1
  
  
  #
  predicted_s2 <- predict(fit_prop_2, newdata=X.tmp.s2, summary=TRUE) %>% 
    mutate(cycle_num = c(target_cycle)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s2 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>% 
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", "")))) %>% 
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("0.2-0.5mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s2,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  ##MPN: NOt exactly sure what you are trying to show with this plot
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  
  
  
  
  final_data_s2 <- bind_rows(final_data_s2, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s2[s,] <-0
  
}

#BEEP to notify when finished lol
beepr::beep(10)

#Save final data
current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s2,here(paste0("data/predicted_og/asv/predicted_og_18s_",current_date,"_s2_phy_all_and_subpools_nofilt_asv.csv")))



# S3 ----------------------------------------------------------------------


######### Final size
##### 1-2mm####
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/asv_level/fido_18s_s3_asv_nofilt.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  group_by(Hash) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Hash")%>%
  #Remove b3 pool
  dplyr::select(c(contains("All"),contains("A1"),contains("C5"),contains("S3")))

#Metadata
meta_18s=read.csv(here("data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt)) %>%
  mutate(Run = case_when(
    grepl("Pool", Sample_name) ~ "2",  # Prioritize "Pool" samples for Run 2
    grepl("\\.1|\\.2", Sample_name) ~ "1",  # Otherwise, assign Run 1 if it ends in .1 or .2
    grepl("\\.3", Sample_name) ~ "2",  # Assign Run 2 if it ends in .3
    TRUE ~ NA_character_  # Default case
  )) %>%
  mutate(Run = factor(Run, levels = c("1", "2")))  # Explicitly define levels
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

#Add in vars for onshore offshore
# read in metadata
metadata=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv"))%>%
  mutate(sample_id=Sample_ID_dot) %>%
  select(sample_id,day_night_0_1, offshore_onshore) %>% 
  distinct(.) 

meta_18s <- meta_18s %>% 
  left_join(metadata, by="sample_id")%>%
  mutate(
    offshore = ifelse(offshore_onshore == "offshore", 1, 0),
    onshore = ifelse(offshore_onshore == "onshore", 1, 0),
    day = ifelse(day_night_0_1 == 0, 1, 0),
    night = ifelse(day_night_0_1 == 1, 1, 0)
  ) %>%
  replace_na(list(offshore = 0, onshore = 0, day = 0, night = 0)) %>% # Replace NA values with 0
  select(-offshore_onshore, -day_night_0_1)  # Drop the original categorical columns


##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),] 

meta_18s=meta_18s%>%
  filter(!is.na(Sample_name))

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+sample_num+Run+offshore+onshore+day+night  -1, data = meta_18s))
Y_s3=fido_input_filt%>% as.matrix() 



#Fit pibble model 
gamma=20

#Specify the remaining priors with default values
upsilon <- nrow(Y_s3)+3 
Omega <- diag(nrow(Y_s3))
G <- cbind(diag(nrow(Y_s3)-1), -1)
Xi <- (upsilon-nrow(Y_s3))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s3)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = gamma*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)
#summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
##Looks ok, centered at zero
##end of added code

##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s3 <- to_clr(fit)

#Convert to Proportions
fit_prop_3 <- to_proportions(fit_s3)





############


X.tmp.s3 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s3) <- rownames(X)
X.tmp.s3["cycle_num", ] <- target_cycle  # evaluate fitted line at target_cycle


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s3 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num","Run2",
                        "offshore","onshore","day","night"))%>% as.data.frame()->samples_to_loop 
#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s3 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s3[s,] <-1
  
  
  #
  predicted_s3 <- predict(fit_prop_3, newdata=X.tmp.s3, summary=TRUE) %>% 
    mutate(cycle_num = c(target_cycle)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s3 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>% 
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", "")))) %>% 
    rownames_to_column("coord") %>%
    pivot_longer(cols = c(-coord),
                 names_to = "replicate",
                 values_to = "n_reads") %>%
    mutate(size=rep("0.2-0.5mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s3,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  ##MPN: NOt exactly sure what you are trying to show with this plot
  sample_temp_sel%>% 
    filter(coord %in% taxa_sel) %>% 
    ggplot(.,aes(x=cycle_num,y=n_reads, fill=coord))+
    geom_line(aes(color=coord), size=2)+
    geom_point(aes(color=coord), shape=5)+
    facet_wrap(~coord, nrow=3)+
    theme_classic()
  
  
  
  
  final_data_s3 <- bind_rows(final_data_s3, sample_temp_sel)
  
  #Clear X.tmo
  X.tmp.s3[s,] <-0
  
}

#BEEP to notify when finished lol
beepr::beep(12)

#Save final data
current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s3,here(paste0("data/predicted_og/asv/predicted_og_18s_",current_date,"_s3_phy_all_and_subpools_nofilt_asv.csv")))



# Compile AEs -------------------------------------------------------------



fit_s1_df_all_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="All_S1_18S")
fit_s2_df_all_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="All_S2_18S")
fit_s3_df_all_18s=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="All_S3_18S")


all_fits_18s=rbind(fit_s1_df_all_18s,fit_s2_df_all_18s,fit_s3_df_all_18s) %>% 
  mutate(pool_type = str_extract(pool, "^[^_]+"),
         Size = str_extract(pool, "(S1|S2|S3)"),
         primer = str_extract(pool, "([^_]+$)"),
         size_fraction= case_when(
           Size == "S1" ~ 0.2,
           Size == "S2" ~ 0.5,
           Size == "S3" ~ 1)) %>%
  mutate(rank_pool = case_when(
    pool_type == "All" ~ 1,
    pool_type == "A1" ~ 2,
    pool_type == "B3" ~ 3,
    pool_type == "C5" ~ 4,
    TRUE ~ NA_integer_  # Handle any other cases
  ))


all_fits_18s %>%
  ggplot(., aes(x = Lambda.mean, y = Lambda.coord)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(pool_type)), height = 0.2, size=4) +
  geom_point(aes(color = as.factor(pool_type)),size = 8, alpha=0.8) + 
  labs(title="18S Family Amplification Efficiencies",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Pool") +
  facet_wrap(~size_fraction, nrow=3)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 

#By taxa
amp_effs_all_and_subpools_by_taxa_18s=all_fits_18s %>%
  ggplot(., aes(x = as.factor(size_fraction), y = Lambda.mean)) + 
  geom_hline(aes(yintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_point(aes(color = as.factor(rank_pool),), size = 8, alpha=1, stroke=2) + 
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5,color = as.factor(rank_pool)), size=1, width=0.2) +
  labs(title = "18S Family Amplification Efficiencies",
       x = "Size Fraction",
       y = "Centered Log Ratio",
       color = "Pool") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-0.25, 0.25, by = 0.05)) +
  scale_x_discrete(labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm"))+
  facet_wrap(~Lambda.coord) +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 16),
        axis.title.x =element_text(size = 16),
        axis.title.y =element_text(size = 16),
        legend.title = element_text(size = 14),  # Increase legend title size
        legend.text = element_text(size = 14))+  # Increase legend entries size+
  scale_color_manual(values=c("#FFABAB","#c996d4", "#93d182", "#a8bbe3"),labels=c("All","Pool 1","Pool 2","Pool 3"))
amp_effs_all_and_subpools_by_taxa_18s

#18S
all_fits_18s=rbind(fit_s1_df_all_18s,fit_s2_df_all_18s,fit_s3_df_all_18s) %>% 
  mutate(pool_type = str_extract(pool, "^[^_]+"),
         Size = str_extract(pool, "(S1|S2|S3)"),
         primer = str_extract(pool, "([^_]+$)"),
         size_fraction= case_when(
           Size == "S1" ~ 0.2,
           Size == "S2" ~ 0.5,
           Size == "S3" ~ 1)) %>%
  mutate(rank_pool = case_when(
    pool_type == "All" ~ 1,
    pool_type == "A1" ~ 2,
    pool_type == "B3" ~ 3,
    pool_type == "C5" ~ 4,
    TRUE ~ NA_integer_  # Handle any other cases
  ))

write.csv(all_fits_18s,here("data/amp_effs/all_amp_effs_18s_all_sub.csv"))
