
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



###Load in the filtered data for the 18s primer using long format species and hash name so I ca identify taxa
##First Size 1

fido_input_filt=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("PCR_bias_correction/data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),] 

meta_18s=meta_18s%>%
  filter(!is.na(Sample_name))

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
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
#400 seems good based on LML
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


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s1 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

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
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s1 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
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
write.csv(final_data_s1,here(paste0("PCR_bias_correction/data/predicted_og/predicted_og_18s_",current_date,"_s1_phy_all_and_subpools.csv")))







# Repeat for other sizes

############First 0.5-1############
#Phyloseq Filtered
fido_input_filt=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("PCR_bias_correction/data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))


##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 



#Fit pibble model 
##MPN: Please remind me how did you choose the 20? Was it using the log marginal likelihood? If so, that code should probably be included here. Happy to chat about this more.
##MPN: This is assuming the default priors for Theta, upsilon, and Xi. Probably reasonable here, but, may want to look in prior predictive checks
##Basically, would run this. These first few rows are just setting the defaults (which fido auto does in the line you have)
upsilon <- nrow(Y_s2)+3 
Omega <- diag(nrow(Y_s2))
G <- cbind(diag(nrow(Y_s2)-1), -1)
Xi <- (upsilon-nrow(Y_s2))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s2)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = 20*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)
#summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
##Looks ok, centered at zero
##end of added code

##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s2 <- to_clr(fit)

#Convert to Proportions
fit_prop_1 <- to_proportions(fit_s2)


############




#Select sample to predict on
#Sample select
sample_sel="sample_numC1.T7.H9_s2"



X.tmp.s2 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s2) <- rownames(X)


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s2 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s2 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s2[s,] <-1
  
  
  #
  predicted_s2 <- predict(fit_prop_1, newdata=X.tmp.s2, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("0.5-1mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s2 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
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
beepr::beep(12)

#Save final data
current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s2,here(paste0("PCR_bias_correction/data/predicted_og/predicted_og_18s_",current_date,"_s2_phy_all_and_subpools.csv")))



# S3 ----------------------------------------------------------------------


######### Final size
##### 1-2mm####
#Phyloseq Filtered
fido_input_filt=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Family")


#Metadata
meta_18s=read.csv(file.path("PCR_bias_correction/data/fido/meta_18s_unaveraged_all.csv"), header=TRUE) %>% 
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))


##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s3=fido_input_filt%>% as.matrix() 


upsilon <- nrow(Y_s3)+3 
Omega <- diag(nrow(Y_s3))
G <- cbind(diag(nrow(Y_s3)-1), -1)
Xi <- (upsilon-nrow(Y_s3))*G%*%Omega%*%t(G)
Theta <- matrix(0, nrow(Y_s3)-1, nrow(X))
priors <- pibble(NULL, X, Gamma = 20*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
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
fit_prop_1 <- to_proportions(fit_s3)


############




#Select sample to predict on
#Sample select



X.tmp.s3 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s3) <- rownames(X)


#Samples to loop thru-for each iteration of the loop I will set the one I want to predict on to '1' from '0'
X.tmp.s3 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

#Create a dataframe to fill with Cycle 0 proportions thru the loop
final_data_s3 <- data.frame()


#Here begins the loop
for(s in samples_to_loop$sample){
  #Print sample name as a sanity check
  print(s)
  
  #Set selected sample to 1
  X.tmp.s3[s,] <-1
  
  
  #
  predicted_s3 <- predict(fit_prop_1, newdata=X.tmp.s3, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("1-2mm"))%>%
    mutate(coord = str_replace(coord, "^prop_", "")) %>%
    rename(n_reads = mean) %>%
    mutate(replicate=rep(paste("predicted",c(str_replace(s, "^sample_num", "")))))
  
  #Compare with original count data after 30 cycles
  Y_s3 %>% as.data.frame() %>% 
    #convert to proportions
    mutate(across(everything(), ~ ./sum(.))) %>%
    dplyr::select(starts_with(c(str_replace(s, "^sample_num", ""))))%>%
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
write.csv(final_data_s3,here(paste0("PCR_bias_correction/data/predicted_og/predicted_og_18s_",current_date,"_s3_phy_all_and_subpools.csv")))

