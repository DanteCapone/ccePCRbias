library(fido)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(gridExtra)

set.seed(5903)



# All Pools ---------------------------------------------------------------

## 18S ---------------------------------------------------------------------


###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

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

X[,1:5]
Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_all_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num") %>%
  mutate(pool="All_S1_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1 <- plot(fit_s1, par = "Lambda", focus.cov = "cycle_num") +
  labs(title = "0.2-0.5 mm") + 
  theme(legend.position = "none",  # Remove legend
        axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),   # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))   # Increase y-axis tick label size
p_s1



#Size 2

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/meta_18s_unaveraged_s2.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_all_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="All_S2_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s_s3=read.csv(file.path("data/fido/meta_18s_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s_s3 <- meta_18s_s3[match(colnames(fido_input_filt_s3), meta_18s_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_all_18s=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="All_S3_18S")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size

p_s3
#Plot all 
# Combine legends
amp_effs_18s=grid.arrange(p_s1,p_s2,p_s3)


#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s.png"),
  plot = amp_effs_18s,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s.pdf"),
  plot = amp_effs_18s,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)






# COI ---------------------------------------------------------------------
###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa

#Phyloseq Filtered


fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_genus_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_all_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="All_S1_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.2-0.5 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
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

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_all_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="All_S2_COI")


# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3
#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_taxa_phy.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi_s3=read.csv(file.path("data/fido/meta_coi_unaveraged_s3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi_s3 <- meta_coi_s3[match(colnames(fido_input_filt_s3), meta_coi_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_all_coi=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="All_S3_COI")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_coi=grid.arrange(p_s1,p_s2,p_s3)


#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi.png"),
  plot = amp_effs_coi,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi.pdf"),
  plot = amp_effs_coi,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)




#  Subpools ---------------------------------------------------------------

# COI Subpools ------------------------------------------------------------

# A1-A3 -------------------------------------------------------------------


fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s1_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s1_a1_a3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_a1_a3_coi=summary(fit_s1, pars = "Lambda") %>% 
  as.data.frame() %>% 
  filter(Lambda.covariate=="cycle_num") %>% 
  mutate(pool="A1_A3_S1_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s2_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s2_a1_a3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_a1_a3_coi=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="A1_A3_S2_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s3_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s3_a1_a3.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi_s3 <- meta_coi_s3[match(colnames(fido_input_filt_s3), meta_coi_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_a1_a3_coi=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="A1_A3_S3_COI")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_coi_a1_a3=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_a1_a3.png"),
  plot = amp_effs_coi_a1_a3,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_a1_a3.pdf"),
  plot = amp_effs_coi_a1_a3,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)




# B3-B5 -------------------------------------------------------------------



fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s1_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s1_b3_b5.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_b3_b5_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S1_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
  labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s2_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s2_b3_b5.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_b3_b5_coi=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S2_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size

p_s2

#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s3_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s3_b3_b5.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s3=fido_input_filt%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_b3_b5_coi=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S3_COI")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_coi_b3_b5=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_b3_b5.png"),
  plot = amp_effs_coi_b3_b5,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_b3_b5.pdf"),
  plot = amp_effs_coi_b3_b5,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)






# C5-C7 -------------------------------------------------------------------

fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s1_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s1_c5_c7.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_c5_c7_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="C5_C7_S1_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
  labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s2_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s2_c5_c7.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_c5_c7_coi=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="C5_C7_S2_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_coi_s3_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/sub_pools/meta_coi_unaveraged_s3_c5_c7.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi_s3 <- meta_coi_s3[match(colnames(fido_input_filt_s3), meta_coi_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_c5_c7_coi=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="C5_C7_S3_COI")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_coi_c5_c7=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_c5_c7.png"),
  plot = amp_effs_coi_c5_c7,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_coi_c5_c7.pdf"),
  plot = amp_effs_coi_c5_c7,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





# 18s ---------------------------------------------------------------------

# A1-A3 -------------------------------------------------------------------

fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_a1_a3_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="A1_A3_S1_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
  labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_a1_a3_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="A1_A3_S2_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_a1_a3.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s3=fido_input_filt%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_a1_a3_18s=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="A1_A3_S3_18S")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_18s_a1_a3=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_a1_a3.png"),
  plot = amp_effs_18s_a1_a3,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_a1_a3.pdf"),
  plot = amp_effs_18s_a1_a3,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





# B3-B5 -------------------------------------------------------------------

fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_b3_b5_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S1_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
  labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_b3_b5_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S2_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_b3_b5.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s3=fido_input_filt%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_b3_b5_18s=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="B3_B5_S3_18S")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")

#Plot all 
# Combine legends
amp_effs_18s_b3_b5=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_b3_b5.png"),
  plot = amp_effs_18s_b3_b5,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_b3_b5.pdf"),
  plot = amp_effs_18s_b3_b5,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)




# C5-C7 -------------------------------------------------------------------
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s1_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y_s1[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_c5_c7_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="C5_C7_S1_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
focus.coord=paste0("clr_",rownames(Y_s1))

# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num",focus.coord=focus.coord)+
  labs(title="0.2-0.5 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s2_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_c5_c7_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="C5_C7_S2_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s2

#Size 3

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/sub_pools/fido_18s_s3_ecdf_taxa_phy_c5_c7.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s<- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s3=fido_input_filt%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_c5_c7_18s=summary(fit_s3, pars = "Lambda") %>% 
  as.data.frame() %>% 
  filter(Lambda.covariate=="cycle_num") %>% 
  mutate(pool="C5_C7_S3_18S")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3

#Plot all 
# Combine legends
amp_effs_18s_c5_c7=grid.arrange(p_s1,p_s2,p_s3)

#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_c5_c7.png"),
  plot = amp_effs_18s_c5_c7,
  width = 10,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/fido_amp_effs_18s_c5_c7.pdf"),
  plot = amp_effs_18s_c5_c7,
  width = 10,  # Width in inches
  height = 8  # Height in inches
)





# Add in Predictions made using all and subpools --------------------------


# 18S ---------------------------------------------------------------------

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s1=fido_input_filt%>% as.matrix() 

X[,1:5]
Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_all_and_subpools_18s=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num") %>%
  mutate(pool="AllandSub_S1_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1 <- plot(fit_s1, par = "Lambda", focus.cov = "cycle_num") +
  labs(title = "0.2-0.5 mm") + 
  theme(legend.position = "none",  # Remove legend
        axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),   # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))   # Increase y-axis tick label size
p_s1



#Size 2

#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s <- meta_18s[match(colnames(fido_input_filt), meta_18s$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s))
Y_s2=fido_input_filt%>% as.matrix() 

X[,1:5]

fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_all_and_subpools_18s=summary(fit_s2, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="AllandSub_S2_18S")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3

#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Family")

#Metadata
meta_18s_s3=read.csv(file.path("data/fido/sub_pools/meta_18s_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_18s and fido_input_filt
meta_18s_s3 <- meta_18s_s3[match(colnames(fido_input_filt_s3), meta_18s_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_18s_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_all_and_subpools_18s=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>%
  mutate(pool="AllandSub_S3_18S")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size

p_s3
#Plot all 
#18S
all_fits_18s=rbind(fit_s1_df_all_and_subpools_18s,fit_s2_df_all_and_subpools_18s,fit_s3_df_all_and_subpools_18s) %>% 
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
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(pool_type)), height = 0.2) +
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
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5,color = as.factor(rank_pool)), width=0.2) +
  labs(title = "COI Genus Amplification Efficiencies",
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




# COI ---------------------------------------------------------------------
###Load in the ECDF-filtered data for the 18S primer using long format species and hash name so I ca identify taxa

#Phyloseq Filtered


fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s1=fido_input_filt%>% as.matrix() 


Y[1:5,1:5]

fit_s1 <- pibble(Y_s1, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s1 <- to_clr(fit_s1)
fit_s1_df_all_and_subpools_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="AllandSub_S1_COI")

# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 


# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s1=plot(fit_s1, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.2-0.5 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s1


#Size 2
#Phyloseq Filtered
fido_input_filt=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi=read.csv(file.path("data/fido/meta_coi_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt))
colnames(fido_input_filt) <- gsub("^X", "", colnames(fido_input_filt))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi <- meta_coi[match(colnames(fido_input_filt), meta_coi$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi))
Y_s2=fido_input_filt%>% as.matrix() 



fit_s2 <- pibble(Y_s2, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s2 <- to_clr(fit_s2)
fit_s2_df_all_and_subpools_coi=summary(fit_s1, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="AllandSub_S2_COI")


# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 





# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s2=plot(fit_s2, par="Lambda", focus.cov="cycle_num")+
  labs(title="0.5-1 mm")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size



#Size 3
#Phyloseq Filtered
fido_input_filt_s3=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv"), header=TRUE, check.names = FALSE, row.names = 1) %>%
  column_to_rownames("Genus")

#Metadata
meta_coi_s3=read.csv(file.path("data/fido/meta_coi_unaveraged_all.csv"), header=TRUE) %>%
  select(-c(X)) %>%
  filter(Sample_name %in% colnames(fido_input_filt_s3))
colnames(fido_input_filt_s3) <- gsub("^X", "", colnames(fido_input_filt_s3))

##Next, we need to make sure that the orders are the same between meta_coi and fido_input_filt
meta_coi_s3 <- meta_coi_s3[match(colnames(fido_input_filt_s3), meta_coi_s3$Sample_name),]

#Model matrix
#This will fit a linear model with an intercept for every sample (no global intercept because of the "-1") and a slope for cycle number
X <- t(model.matrix(~ cycle_num+ sample_num  -1, data = meta_coi_s3))
Y_s3=fido_input_filt_s3%>% as.matrix() 



fit_s3 <- pibble(Y_s3, X, Gamma = 20*diag(nrow(X)), n_samples = 10000)

fit_s3 <- to_clr(fit_s3)
fit_s3_df_all_and_subpools_coi=summary(fit_s3, pars = "Lambda") %>% as.data.frame() %>% filter(Lambda.covariate=="cycle_num")%>% 
  mutate(pool="AllandSub_S3_COI")




# pull out indices for random intercepts corresponding to `sample_num`
focus.covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

# Also just so the plot fits nicely in Rmarkdown we are also going to just 
# plot a few of the taxa
# focus.coord <- paste0("clr_", c("S.gallolyticus", "R.intestinalis", "L.ruminis")) 

# Also to make the plot fit nicely, I just flip the orientation of the plot 
plot(fit_s3, par="Lambda", focus.cov=focus.covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate~.)



# Also to make the plot fit nicely, I just flip the orientation of the plot 
p_s3=plot(fit_s3, par="Lambda", focus.cov="cycle_num")+
  labs(title="1-2 mm") + theme(legend.position = "none")+  # Set axis labels
  theme(axis.title.x = element_text(size = 14),  # Increase x-axis label size
        axis.text.x = element_text(size = 14),    # Increase x-axis tick label size
        axis.text.y = element_text(size = 14))    # Increase y-axis tick label size
p_s3+scale_color_brewer("Set2")


### 

# 18S Join all pools: All, A1-A3, B3-B5, C5-C7 ----------------------------------------------------------

#18S
all_fits_18s=rbind(fit_s1_df_all_18s,fit_s2_df_all_18s,fit_s3_df_all_18s,
               fit_s1_df_a1_a3_18s,fit_s2_df_a1_a3_18s,fit_s3_df_a1_a3_18s,
               fit_s1_df_b3_b5_18s,fit_s2_df_b3_b5_18s,fit_s3_df_b3_b5_18s,
               fit_s1_df_c5_c7_18s,fit_s2_df_c5_c7_18s,fit_s3_df_c5_c7_18s,
               fit_s1_df_all_and_subpools_18s,fit_s2_df_all_and_subpools_18s,
               fit_s3_df_all_and_subpools_18s)%>%
  mutate(
    pool_type = str_extract(pool, "^[^_]+"),
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
    pool_type == "AllandSub"~5,
    TRUE ~ NA_integer_  # Handle any other cases
  ))

write.csv(all_fits_18s,here("data/amp_effs/all_amp_effs_18s_all_sub.csv"))


#Load data
all_fits_18s=read.csv(here("data/amp_effs/all_amp_effs_18s_all_sub.csv"))
# Plot
all_fits_18s %>%
ggplot(., aes(x = Lambda.mean, y = Lambda.coord)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(pool_type)), height = 0.2) +
  geom_point(aes(color = as.factor(pool_type)),size = 8, alpha=0.8) + 
  labs(title="18S Family Amplification Efficiencies",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Pool") +
  facet_wrap(~size_fraction, nrow=3)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 

#By taxa
amp_effs_all_by_taxa_18s <- all_fits_18s %>%
  filter(rank_pool != "5") %>% 
  ggplot(aes(x = as.factor(size_fraction), y = Lambda.mean)) + 
  geom_hline(aes(yintercept = 0), color = "black", alpha = 1, size = 2) +
  geom_point(data = . %>% filter(pool_type != "AllandSub"), 
             aes(color = as.factor(rank_pool)), size = 8, alpha = 0.9, stroke = 2) +
  geom_point(data = . %>% filter(pool_type == "AllandSub"),
             aes(), size = 4, shape = 18, fill = "black",stroke=1, alpha=0.9) +
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5, color = as.factor(rank_pool)), width = 0.2) +
  labs(title = "18S Family Amplification Efficiencies",
       x = "Size Fraction",
       y = "Centered Log Ratio",
       color = "Pool") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-0.25, 0.25, by = 0.05)) +
  scale_x_discrete(labels = c("0.2-0.5 mm", "0.5-1 mm", "1-2 mm")) +
  facet_wrap(~Lambda.coord) +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 14)) +
  scale_color_manual(values = c("#FFABAB", "#c996d4", "#93d182", "#a8bbe3", "#000000"), 
                     labels = c("All Pools","Pool 1", "Pool 2", "Pool 3", "All+Subpools"))

amp_effs_all_by_taxa_18s
#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_18s.png"),
  plot = amp_effs_all_by_taxa_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa.pdf"),
  plot = amp_effs_all_by_taxa_18s,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)




# COI Join all pools:  All, A1-A3, B3-B5, C5-C7 ------------------------------------------------------
all_fits_coi=rbind(fit_s1_df_all_coi,fit_s2_df_all_coi,fit_s3_df_all_coi,
               fit_s1_df_a1_a3_coi,fit_s2_df_a1_a3_coi,fit_s3_df_a1_a3_coi,
               fit_s1_df_b3_b5_coi,fit_s2_df_b3_b5_coi,fit_s3_df_b3_b5_coi,
               fit_s1_df_c5_c7_coi,fit_s2_df_c5_c7_coi,fit_s3_df_c5_c7_coi,
               fit_s1_df_all_and_subpools_coi,fit_s2_df_all_and_subpools_coi,fit_s3_df_all_and_subpools_coi)%>%
  mutate(
    pool_type = str_extract(pool, "^[^_]+"),
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
  )) %>% 
  mutate(Lambda.coord = ifelse(Lambda.coord == "clr_Sagittidae", "clr_unidentified Sagittidae", Lambda.coord)) 
  

write.csv(all_fits_coi,here("data/amp_effs/all_amp_effs_coi_all_sub.csv"))



# Plot
all_fits_coi %>%
  ggplot(., aes(x = Lambda.mean, y = Lambda.coord)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_errorbarh(aes(xmin = Lambda.p2.5, xmax = Lambda.p97.5,color = as.factor(rank_pool)), height = 0.2) +
  geom_point(aes(color = as.factor(rank_pool)),size = 8, alpha=0.8) + 
  labs(title="18S Family Amplification Efficiencies",x = "Offfshore \u2190 Centered Log-Ratio(PC1) \u2192 Onshore",
       y="", color = "Pool") +
  facet_wrap(~size_fraction, nrow=3)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 12)) 

#By taxa
amp_effs_all_by_taxa_coi=all_fits_coi %>%
  filter(rank_pool != "5") %>% 
  ggplot(., aes(x = as.factor(size_fraction), y = Lambda.mean)) + 
  geom_hline(aes(yintercept = 0), color = "black", alpha=0.5, size=2)+
  geom_point(aes(color = as.factor(rank_pool),), size = 8, alpha=1.4, stroke=2) + 
  geom_errorbar(aes(ymin = Lambda.p2.5, ymax = Lambda.p97.5,color = as.factor(rank_pool)), width=0.2) +
  labs(title = "COI Genus Amplification Efficiencies",
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
  scale_color_manual(values=c("#FFABAB","#c996d4", "#93d182", "#a8bbe3"),labels=c("All","Pool 1","Pool 2","Pool 3", "All and Subpools"))
# Adjust font size for y-axis tick labels
amp_effs_all_by_taxa_coi
#PNG & PDF Save
ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_coi.png"),
  plot = amp_effs_all_by_taxa_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)

ggsave(
  filename = here("plots/pre_processing/amp_effs_all_by_taxa_coi.pdf"),
  plot = amp_effs_all_by_taxa_coi,
  width = 16,  # Width in inches
  height = 12  # Height in inches
)




