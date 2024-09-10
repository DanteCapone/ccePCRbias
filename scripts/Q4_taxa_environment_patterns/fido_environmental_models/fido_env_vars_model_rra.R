#Script to loop through different environmental variables and compare against taxa from quantitative analysis



# Load libraries ----------------------------------------------------------
library(phyloseq)
library(tidyverse)
library(fido)
library(here)
library(ggplot2)
library(gridExtra)
source(here("scripts/helpful_functions/fido_missing_funs.R"))

# Raw Reads ---------------------------------------------------------------


# 18S ---------------------------------------------------------------------

otu18s1=read.csv(here("data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv")) %>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otu18s2=read.csv(here("data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv"))%>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otu18s3=read.csv(here("data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv"))%>%
  select(starts_with("C"),-X,Family)%>%
  pivot_longer(-Family, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Family) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Family") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)


zhan_taxa=read.csv(here("data/phyloseq_bio_data/18S/fido_18s_family_tax_table.csv"))  %>% 
  #Filter weird taxa
  filter(!(Family == "Corycaeidae" & Order =="Cyclopoida")) %>% 
  select(-X) %>% 
  distinct() %>% 
  mutate(Family2=Family) %>%
  column_to_rownames("Family2")
zhan_taxa=tax_table(as.matrix(zhan_taxa))



meta18s=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(c(-Sizefractionmm,-offshore_onshore,-clust_group,-cycle, -max_size)) %>%
  sample_data(.)

#Taking the mean of replicates
dat_1=phyloseq(otu18s1,zhan_taxa,meta18s) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otu18s2,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otu18s3,zhan_taxa,meta18s)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)






set.seed(899)

#Loop thru fido models for each size
counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-mixedlayerdepths,
             -beam_depth,-chl_max,-intergrated_chl,-PAR_1_depth_adj,-day_night_0_1,
             -PC1)) %>% 
    mutate(distance_from_shore=log(distance_from_shore),
           chl_max_depth=log(chl_max_depth),
           nitracline_depth=log(nitracline_depth),
           hypoxia_depth=log(hypoxia_depth))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(formula_obj, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  #Determine Gamma
  gamma <- c(1,2,3,5,8,10,15,20,50,100,500,700,1000)
  logML <- rep(NA, length(gamma))
  for(i in 1:length(gamma)){
    fit <- pibble(Y, X, Gamma = gamma[i]*diag(nrow(X)), n_samples=5000)
    logML[i] <- fit$logMarginalLikelihood
    print(i)
  }
  plot(gamma, logML, type = "l")
  points(gamma, logML)
  gamma=20
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- gamma*diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  
  tax <- tax_table(dat)[,c("Family")] %>% as.data.frame() %>%
    rownames_to_column("Family2")%>% select(Family)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  
  #Add size column
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  posterior_summary=posterior_summary %>% 
    mutate(size=sizes[counter])
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
  assign(paste0("focus_", counter), focus)
  
  #Plot model wiht fido plot fxn
  fido.p=plot(posterior, par="Lambda", focus.cov = rownames(X)[2:nrow(X)])
  assign(paste0("fido.p_S", counter), fido.p)
  if (nrow(focus) == 0) {
    # Start next loop
    # Your code for the next loop goes here
  } else {
    # Continue with the rest of your code using the 'focus' dataframe
    ##This code will only work if there is a sig. result returned
    focus <- unique(focus$coord)
    # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p25) == sign(posterior_summary$p75)]]
    focus_cov=rownames(X)
    focus_cov[2]
    
    
    
    
  }
}

# grid.arrange(fido.p_S1,fido.p_S2,fido.p_S3, nrow=3)
# focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
# focus <- unique(focus$coord)
# plot(posterior, par="Lambda", focus.cov = rownames(X)[2:8])

posterior_summary_all=rbind(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3) 


# Modify the facet title
facet_titles <- posterior_summary_all %>%
  filter(sign(p25) == sign(p75)) %>%
  select(covariate) %>%
  distinct() %>%
  arrange(covariate) %>%
  mutate(label = case_when(
    covariate == "nitracline_depth" ~ "log(nitracline_depth)",
    covariate == "distance_from_shore" ~ "log(distance_from_shore)",
    covariate == "chl_max_depth" ~ "log(chl_max_depth)",
    covariate == "hypoxia_depth" ~ "log(hypoxia_depth)",
    TRUE ~ covariate
  )) %>%
  deframe()


posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  filter(covariate != "(Intercept)") %>%
  ggplot(aes(x = mean, y = coord, color = size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2) +
  geom_errorbarh(aes(xmin = p25, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape = size), size = 3, alpha = 0.8) + 
  facet_wrap(~covariate, nrow = 4, labeller = as_labeller(facet_titles)) +
  labs(
    title = "Pibble Model for PCR Bias-Corrected 18S Taxa",
    x = "Centered Log-Ratio",
    y = "",
    color = "Size Fraction (mm)",
    shape = "Size Fraction (mm)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) -> all_plot_18s_rra

all_plot_18s_rra


## As a heatmap
facet_titles <- posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  select(covariate) %>%
  distinct() %>%
  arrange(covariate) %>%
  mutate(label = case_when(
    covariate == "nitracline_depth" ~ "log(nitracline_depth)",
    covariate == "distance_from_shore" ~ "log(distance_from_shore)",
    covariate == "chl_max_depth" ~ "log(chl_max_depth)",
    TRUE ~ covariate
  )) %>%
  deframe()

# Filter the posterior summary to find significant cells
significant_cells <- posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  mutate(coord_size = paste(coord, size, sep = "_"))


# Modify the data to combine coord and size
data_modified <- posterior_summary_all %>%
  filter(sign(p25) == sign(p75)) %>%
  filter(covariate != "(Intercept)") %>%
  mutate(coord_size = paste(coord, size, sep = "_")) %>%
  mutate(covariate = fct_relevel(covariate, names(facet_titles)))

# Get unique coord values for horizontal lines
coord_levels <- data_modified %>%
  arrange(coord) %>%
  pull(coord) %>%
  unique()

# Plot the heatmap with black lines below each coord group
ggplot(data_modified, aes(x = covariate, y = coord_size, fill = mean)) + 
  geom_tile(color = NA) +
  geom_hline(yintercept = seq(0.5, length(unique(data_modified$coord_size)) - 0.5, by = 1), color = "black", size = 0.5) +
  geom_vline(xintercept = seq(0.5, length(unique(data_modified$covariate)) - 0.5, by = 1), color = "black", size = 0.5) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "grey") +  geom_hline(yintercept = which(levels(factor(data_modified$coord_size)) %in% paste0(coord_levels, "_", unique(data_modified$size)[1])) - 0.5, color = "black") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "grey") +
  labs(
    title = "Heatmap of CLR Values for PCR Bias-Corrected 18S Taxa",
    x = "Variable",
    y = "Family and Size Class",
    fill = "Centered Log-Ratio"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 20)
  )+
  geom_text(data = significant_cells, aes(x = covariate, y = coord_size, label = "*"), color = "black", size = 10, vjust = 0.5)-> heatmap_plot_18s_rra

heatmap_plot_18s_rra

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa_rra_heatmap_50_sig.pdf"), 
    plot = heatmap_plot_18s_rra,
    width = 16,  # Width in inches
    height = 16  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_18s_taxa_rra_heatmap_50_sig.png"), 
    plot = heatmap_plot_18s_rra,
    width = 16,  # Width in inches
    height = 16  # Height in inches
  )}



# COI ---------------------------------------------------------------------
#Using agglomerated taxa from fido
otucoi1=read.csv(here("data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv")) %>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi2=read.csv(here("data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv"))%>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)

otucoi3=read.csv(here("data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv"))%>%
  select(starts_with("C"),-X,Genus)%>%
  pivot_longer(-Genus, names_to = "sample", values_to = "value") %>%
  mutate(Sample_ID = str_sub(sample, end = -3)) %>%
  group_by(Sample_ID,Genus) %>%
  summarise(value = round(mean(value, na.rm = TRUE), 0)) %>%
  pivot_wider(names_from = Sample_ID, values_from = value) %>%
  column_to_rownames("Genus") %>%
  otu_table(as.matrix(.), taxa_are_rows = TRUE)


#Taxa files

coi_taxa=read.csv(here("data/phyloseq_bio_data/coi/fido_coi_genus_tax_table.csv"))  %>% 
  #Filter weird taxa
  select(-X) %>% 
  distinct() %>% 
  # Remove duplicated Genera
  group_by(Genus) %>%
  filter(n() == 1) %>%
  mutate(Genus2=Genus) %>%
  column_to_rownames("Genus2")
coi_taxa=tax_table(as.matrix(coi_taxa))



metacoi=read.csv(here("data/physical_environmental_data/env_metadata_impute_phyloseq_6.9.2023.csv")) %>%
  dplyr::select(-c("X")) %>%
  column_to_rownames("Sample_ID_dot") %>%
  select(c(-Sizefractionmm,-offshore_onshore,-clust_group,-cycle, -max_size)) %>%
  sample_data(.)

#Taking the mean of replicates
dat_1=phyloseq(otucoi1,coi_taxa,metacoi) %>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_2=phyloseq(otucoi2,coi_taxa,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)
dat_3=phyloseq(otucoi3,coi_taxa,metacoi)%>% merge_samples(.,"Sample_ID_short",fun= mean)

data_list <- list(dat_1=dat_1, dat_2=dat_2, dat_3=dat_3)




set.seed(899)

#Loop thru fido models for each size
counter=0
for (dat_name in names(data_list)) {
  
  #Counter
  counter=counter+1
  print(counter)
  dat <- data_list[[dat_name]]
  
  sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
    select(c(-Sample_ID_short,-oxy_sat,-mixedlayerdepths,
             -beam_depth,-chl_max,-intergrated_chl,-PAR_1_depth_adj,-day_night_0_1,
             -PC1)) %>% 
    mutate(distance_from_shore=log(distance_from_shore),
           chl_max_depth=log(chl_max_depth),
           nitracline_depth=log(nitracline_depth),
           hypoxia_depth=log(hypoxia_depth))
  
  
  formula_string <- paste("~", paste(names(sample_dat), collapse = " + "), sep = "")
  formula_obj <- as.formula(formula_string)
  X <- t(model.matrix(formula_obj, data=sample_dat))
  Y <- otu_table(dat) %>% t(.)
  
  #Determine Gamma
  gamma <- c(1,2,3,5,8,10,15,20,50,100,500,700,1000)
  logML <- rep(NA, length(gamma))
  for(i in 1:length(gamma)){
    fit <- pibble(Y, X, Gamma = gamma[i]*diag(nrow(X)), n_samples=5000)
    logML[i] <- fit$logMarginalLikelihood
    print(i)
  }
  plot(gamma, logML, type = "l")
  points(gamma, logML)
  gamma=20
  
  
  ## This is all prior specification
  upsilon <- ntaxa(dat)+3 
  Omega <- diag(ntaxa(dat))
  G <- cbind(diag(ntaxa(dat)-1), -1)
  Xi <- (upsilon-ntaxa(dat))*G%*%Omega%*%t(G)
  Theta <- matrix(0, ntaxa(dat)-1, nrow(X))
  Gamma <- gamma*diag(nrow(X))
  
  ##This code is used to check priors, not for actual model fitting.
  priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)  
  print(priors)
  
  priors <- to_clr(priors)  
  summary(priors, pars="Lambda", gather_prob=TRUE, as_factor=TRUE, use_names=TRUE)  
  
  names_covariates(priors) <- rownames(X)
  priors$Y <- Y # remember pibblefit objects are just lists
  posterior <- refit(priors, optim_method="lbfgs", jitter = 1e-5)
  
  tax <- tax_table(dat)[,c("Genus")] %>% as.data.frame() %>%
    rownames_to_column("Genus2")%>% select(Genus)
  num <- 1:nrow(tax)
  tax <- unname(apply(tax, 1, paste, collapse="_"))
  tax <- paste(tax,sep="_")
  names_categories(posterior) <- tax
  
  ##This is the "now what?" part. We have our model, what does it tell us?
  posterior_summary <- summary(posterior, pars="Lambda")$Lambda
  
  #Add size column
  # Your plotting code
  sizes=c("0.2-0.5 mm","0.5-1 mm", "1-2 mm")
  posterior_summary=posterior_summary %>% 
    mutate(size=sizes[counter])
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
  assign(paste0("focus_", counter), focus)
  
  #Plot model wiht fido plot fxn
  fido.p=plot(posterior, par="Lambda", focus.cov = rownames(X)[2:nrow(X)])
  assign(paste0("fido.p_S", counter), fido.p)
  if (nrow(focus) == 0) {
    # Start next loop
    # Your code for the next loop goes here
  } else {
    # Continue with the rest of your code using the 'focus' dataframe
    ##This code will only work if there is a sig. result returned
    focus <- unique(focus$coord)
    # focus_cov=rownames(X)[rownames(X)==posterior_summary$covariate[sign(posterior_summary$p25) == sign(posterior_summary$p75)]]
    focus_cov=rownames(X)
    focus_cov[2]
    
    
    
    
  }
}


posterior_summary_all=rbind(posterior_summary_S1,posterior_summary_S2,posterior_summary_S3) 


# Modify the facet title
facet_titles <- posterior_summary_all %>%
  filter(sign(p25) == sign(p75)) %>%
  select(covariate) %>%
  distinct() %>%
  arrange(covariate) %>%
  mutate(label = case_when(
    covariate == "nitracline_depth" ~ "log(nitracline_depth)",
    covariate == "distance_from_shore" ~ "log(distance_from_shore)",
    covariate == "chl_max_depth" ~ "log(chl_max_depth)",
    covariate == "hypoxia_depth" ~ "log(hypoxia_depth)",
    TRUE ~ covariate
  )) %>%
  deframe()


posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  filter(covariate != "(Intercept)") %>%
  ggplot(aes(x = mean, y = coord, color = size)) + 
  geom_vline(aes(xintercept = 0), color = "black", alpha=0.5, size=2) +
  geom_errorbarh(aes(xmin = p25, xmax = p97.5), height = 0.2) +
  geom_point(aes(shape = size), size = 3, alpha = 0.8) + 
  facet_wrap(~covariate, nrow = 4, labeller = as_labeller(facet_titles)) +
  labs(
    title = "Pibble Model for PCR Bias-Corrected coi Taxa",
    x = "Centered Log-Ratio",
    y = "",
    color = "Size Fraction (mm)",
    shape = "Size Fraction (mm)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) -> all_plot_coi_rra

all_plot_coi_rra


## As a heatmap
facet_titles <- posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  select(covariate) %>%
  distinct() %>%
  arrange(covariate) %>%
  mutate(label = case_when(
    covariate == "nitracline_depth" ~ "log(nitracline_depth)",
    covariate == "distance_from_shore" ~ "log(distance_from_shore)",
    covariate == "chl_max_depth" ~ "log(chl_max_depth)",
    TRUE ~ covariate
  )) %>%
  deframe()

# Filter the posterior summary to find significant cells
significant_cells <- posterior_summary_all %>%
  filter(sign(p2.5) == sign(p97.5)) %>%
  mutate(coord_size = paste(coord, size, sep = "_"))


# Modify the data to combine coord and size
data_modified <- posterior_summary_all %>%
  filter(sign(p25) == sign(p75)) %>%
  filter(covariate != "(Intercept)") %>%
  mutate(coord_size = paste(coord, size, sep = "_")) %>%
  mutate(covariate = fct_relevel(covariate, names(facet_titles)))

# Get unique coord values for horizontal lines
coord_levels <- data_modified %>%
  arrange(coord) %>%
  pull(coord) %>%
  unique()

# Plot the heatmap with black lines below each coord group
ggplot(data_modified, aes(x = covariate, y = coord_size, fill = mean)) + 
  geom_tile(color = NA) +
  geom_hline(yintercept = seq(0.5, length(unique(data_modified$coord_size)) - 0.5, by = 1), color = "black", size = 0.5) +
  geom_vline(xintercept = seq(0.5, length(unique(data_modified$covariate)) - 0.5, by = 1), color = "black", size = 0.5) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "grey") +  geom_hline(yintercept = which(levels(factor(data_modified$coord_size)) %in% paste0(coord_levels, "_", unique(data_modified$size)[1])) - 0.5, color = "black") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "grey") +
  labs(
    title = "Heatmap of CLR Values for PCR Bias-Corrected coi Taxa",
    x = "Variable",
    y = "Genus and Size Class",
    fill = "Centered Log-Ratio"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 20)
  )+
  geom_text(data = significant_cells, aes(x = covariate, y = coord_size, label = "*"), color = "black", size = 10, vjust = 0.5)-> heatmap_plot_coi_rra

heatmap_plot_coi_rra

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_rra_heatmap_50_sig.pdf"), 
    plot = heatmap_plot_coi_rra,
    width = 16,  # Width in inches
    height = 16  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_rra_heatmap_50_sig.png"), 
    plot = heatmap_plot_coi_rra,
    width = 16,  # Width in inches
    height = 16  # Height in inches
  )}

# -------------------------------------------------------------------------

# 

saving=1
if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_r_ra_95.pdf"), 
    plot = all_plot_coi_rra,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  )}

if (saving==1) {
  ggsave(
    filename = here("plots/Q4_taxa_vars/fido_env_vars_coi_taxa_r_ra_95.png"), 
    plot = all_plot_coi_rra,
    width = 12,  # Width in inches
    height = 16  # Height in inches
  )}

