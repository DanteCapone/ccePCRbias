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


###Load in the ECDF-filtered data for the coi primer using long format species and hash name so I ca identify taxa
##First Size 1
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_s1.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
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
gamma=400
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
#Sample select
sample_sel="sample_numC1.T7.H9_S1"



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
write.csv(final_data_s1,here(paste0("data/predicted_og/predicted_og_coi_",current_date,"_s1_phy.csv")))




# 0.5-1 mm ----------------------------------------------------------------

fido_input_filt=read.csv(file.path("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  column_to_rownames("Genus")
#Metadata

meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
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
priors <- pibble(NULL, X, Gamma = gamma*diag(nrow(X)), upsilon = upsilon, Theta = Theta, Xi = Xi, n_samples = 10000)
print(priors)
priors <- to_clr(priors)
summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
##Looks ok, centered at zero
##end of added code

##MPN: Note, you had lower case "gamma" the parameter is upper case "Gamma". Fido was using the default here instead of what you supplied.
fit <- pibble(Y_s2, X, Gamma = gamma*diag(nrow(X)), n_samples = 10000)

#Convert to centered log ratio coordinates
fit_s2 <- to_clr(fit)

#Convert to Proportions
fit_prop_2 <- to_proportions(fit_s2)

############
#Predict at cycle 0
X.tmp.s2 <- matrix(0, nrow(X), 1) #Create fake covariate data to predict the regression line based on 
rownames(X.tmp.s2) <- rownames(X)


#Samples to loop thru
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
  predicted_s2 <- predict(fit_prop_2, newdata=X.tmp.s2, summary=TRUE) %>% 
    mutate(cycle_num = c(0)[sample])%>%
    mutate(size=rep("0.2-0.5mm"))%>%
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
write.csv(final_data_s2,here(paste0("data/predicted_og/predicted_og_coi_",current_date,"_s2_phy.csv")))


# 1-2 mm  -----------------------------------------------------------------

fido_input_filt=read.csv(file.path("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1)%>%
  column_to_rownames("Genus")


#Metadata
meta_coi=read.csv(file.path("data/fido","meta_coi_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##MPN: Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
##MPN: To be clear, this will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))

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
write.csv(final_data_s3,here(paste0("data/predicted_og/predicted_og_coi_",current_date,"_s3_phy.csv")))








## Plot Amplification efficiencies by taxa


counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(-"Sample_ID_short")
  meta_vars=names(sample_dat[,1:ncol(sample_dat)])
  
  X <- t(model.matrix(~PC1, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  
  tax <- tax_table(dat)[,c("Family","Species")]
  hash <- rownames(tax_table(dat))
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,hash,sep="_")
  names_categories(posterior) <- tax
  
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  assign(paste0("posterior_summary_S", counter), posterior_summary)
  
  ##Let's examine this more
  head(posterior_summary)
  ##Mean is the mean of the posterior samples. You can think of it as the estimated beta for the regression model of that specific taxa.
  ##Covariate: we are fitting y = \beta_0 + \beta_1 * potemp2. So we have estimates both for the intercept and slope.
  ## p2.5, p25, etc. these are the 2.5th, 25th, etc. quantiles of the posterior distribution
  ## NOte that p2.5 and p97.5 would give a 95% credible interval.
  ##So for the first taxa, the intercept has a 95% interval of -2.78,12.5
  ##We assess significance by seeing if zero is in this interval. So, for above, zero is in the inteval, this intercept term isn't significant.
  
  ##Now, we are filtering the posterior summary to significant samples only.
  focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),] 
  # filter(covariate != "(Intercept)")
  # focus <- posterior_summary[sign(posterior_summary$p25) == sign(posterior_summary$p75),]
  
  focus ##Note that there are zero rows --> no evidence of an effect of potemp2 on any of the taxa
  
  ##This code will only work if there is a sig. result returned
  focus_coord <- unique(focus$coord)
  # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5)]]
  focus_cov= "PC1"
  focus_cov
  
  # Create the PDF filename string
  save_path="plots/model_multiple_vars/"
  filename <- paste0(save_path,"S_", dat_name, "_multiple_vars_coi.png")
  
  # For demonstration purposes, print the filename
  cat("Saving to:", filename, "\n")
  
  
  png(filename, width = 1000, height = 600)
  
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  pp=plot(posterior, par="Lambda", focus.coord = focus_coord, focus.cov = focus_cov)+
    labs(title=paste(sizes[counter]))
  pp
  Sys.sleep(2)
  # Close the PDF device
  dev.off()
  
}


#####
#Add size column
posterior_summary_S1 =posterior_summary_S1 %>% mutate(size="0.2-0.5")
posterior_summary_S2 =posterior_summary_S2 %>% mutate(size="0.5-1")
posterior_summary_S3 =posterior_summary_S3 %>% mutate(size="1-2")

posterior_summary_all=bind_rows(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3)%>%
  filter(covariate != "(Intercept)")

#Plot
ggplot(posterior_summary_all, aes(x = mean, y = coord, color=size)) + 
  geom_point(size = 3) + 
  geom_vline(aes(xintercept = 0), color = "black")+
  geom_errorbarh(aes(xmin = p2.5, xmax = p97.5), height = 0.2) +
  labs(x = "Mean Log-Ratio PC1",y="", color = "Size Class") +
  theme_minimal()+
  facet_wrap(~size)
