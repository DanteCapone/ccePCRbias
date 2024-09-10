#Size based fido model
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


###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa
##First Size 1
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_order_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Order")

  #Metadata
  meta_18s=read.csv(file.path("data/fido/meta_18s_unaveraged_s1.csv"), header=TRUE) %>%
    select(-c(X)) %>%
    filter(Sample_name %in% colnames(fido_input_filt))
  colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))
  
  ##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
  meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]
  
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
  #20 seems good based on LML
  gamma=20
  upsilon <- nrow(Y_s1)+3 
  Omega <- diag(nrow(Y_s1))
  G <- cbind(diag(nrow(Y_s1)-1), -1)
  Xi <- (upsilon-nrow(Y_s1))*G%*%Omega%*%t(G)
  Theta <- matrix(0, nrow(Y_s1)-1, nrow(X))
  priors <- pibble(NULL, X, Gamma = gamma*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
  print(priors)
  priors <- to_clr(priors)
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  ##Looks ok, centered at zero
  ##end of added code
  
  ##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
  fit <- pibble(Y_s1, X, Gamma = gamma*diag(nrow(X)), n_samples = 10000)
  
  #Convert to centered log ratio coordinates
  fit_s1 <- to_clr(fit)

  #Convert to Proportions
  fit_prop_1 <- to_proportions(fit_s1)
  
  
  ############




#Select sample to predict on




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
write.csv(final_data_s1,here(paste0("data/predicted_og/predicted_og_18s_",current_date,"_s1_phy_order.csv")))



############### Let's repeat for other sizes now ###############

############First 0.5-1############
fido_input_filt=read.csv(file.path("data/fido/phy/fido_18s_s2_ecdf_order_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  column_to_rownames("Order")
#Metadata
meta_18s=read.csv(file.path("data/fido","meta_18s_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##MPN: Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))

Y_s2=fido_input_filt%>% as.matrix() 


fit <- pibble(Y_s2, X, Gamma = gamma*diag(nrow(X)), n_samples = 10000)

# ,Convert to centered log ratio coordinates
fit_s2 <- to_clr(fit)
###Proportions
fit_prop_2 <- to_proportions(fit_s2)


############
#Predict at cycle 0
X.tmp.s2 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s2) <- rownames(X)


#Samples to loop thru
X.tmp.s2 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

final_data_s2 <- data.frame()

for(s in samples_to_loop$sample){
  print(s)
  X.tmp.s2[s,] <-1
  
  predicted_s2 <- predict(fit_prop_2, newdata=X.tmp.s2, summary=TRUE) %>% 
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
    mutate(size=rep("0.5-1mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s2,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  
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


beepr::beep(4)

current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s2,here(paste0("data/predicted_og/predicted_og_18s_",current_date,"_s2_phy_order.csv")))



######### Final size
##### 1-2mm####
fido_input_filt=read.csv(file.path("data/fido/phy/fido_18s_s3_ecdf_order_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  column_to_rownames("Order")


#Metadata
meta_18s=read.csv(file.path("data/fido","meta_18s_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##MPN: Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))

Y_s3=fido_input_filt%>% as.matrix() 

fit <- pibble(Y_s3, X, Gamma = gamma*diag(nrow(X)), n_samples = 10000)

# ,Convert to centered log ratio coordinates
fit_s3 <- to_clr(fit)

##Proportions
fit_prop_3 <- to_proportions(fit_s3)

#
####Predict at cycle 0
#Make X.tmp for loop
X.tmp.s3 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s3) <- rownames(X)


#Samples to loop thru
X.tmp.s3 %>% as.data.frame() %>% rownames_to_column("sample") %>%
  select("sample") %>%
  filter(!sample %in% c("sample_numCalibration","cycle_num"))%>% as.data.frame()->samples_to_loop 

final_data_s3 <- data.frame()

for(s in samples_to_loop$sample){
  print(s)
  X.tmp.s3[s,] <-1
  
  predicted_s3 <- predict(fit_prop_3, newdata=X.tmp.s3, summary=TRUE) %>% 
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
    mutate(size=rep("1-2mm")) %>%
    mutate(cycle_num = rep(30, nrow(.))) %>%
    bind_rows(predicted_s3,.)%>%
    group_by(cycle_num) %>%
    arrange(desc(cycle_num),desc(n_reads))->sample_temp_sel
  
  
  taxa_list=unique(sample_temp_sel$coord)
  taxa_sel=taxa_list[1:3]
  
  
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


beepr::beep(7)

current_date <- format(Sys.Date(), "%m_%d_%Y")
write.csv(final_data_s3,here(paste0("data/predicted_og/predicted_og_18s_",current_date,"_s3_phy_order.csv")))



#Archive older files
# Directory containing the files
file_dir <- here("data/predicted_og")

# Get a list of files in the directory
files <- list.files(file_dir, pattern = "predicted_og_18s_", full.names = TRUE)

# Archive files with an older date
for (file in files) {
  date_str <- str_extract(basename(file), "(?<=predicted_og_18s_)[0-9_]+(?=_s\\d_phy_order\\.csv)")
  if (!is.na(date_str) && date_str != "") {
    # Convert the date string from format mm_dd_yyyy to Date object
    file_date <- as.Date(date_str, format = "%m_%d_%Y")
    if (!is.na(file_date) && file_date < as.Date(current_date, "%m_%d_%Y")) {
      # Check if the "past" directory exists, create it if it doesn't
      past_dir <- file.path(file_dir, "past")
      if (!dir.exists(past_dir)) {
        dir.create(past_dir)
      }
      # Move the file to the "past" directory
      file.rename(file, file.path(past_dir, basename(file)))
    }
  }
}

