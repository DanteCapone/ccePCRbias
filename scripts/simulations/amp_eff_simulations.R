#Fido simulation 

# Load necessary libraries
library(tidyverse)

# Define parameters
set.seed(42)  # For reproducibility

# Number of taxa and PCR cycles
n_taxa <- 3
n_cycles <- 40

# Taxa names
taxa <- c("A", "B", "C")

# Amplification efficiencies for each taxon
amplification_eff <- c(0.65, 0.5, 0.55)  # a values (constant across cycles for each taxon)

initial_concentrations <- runif(n_taxa, min = 1e3, max = 1e4)  # Random initial values

# Create a simulation dataframe
sim_data <- expand.grid(
  Taxon = taxa,
  Cycle = 0:n_cycles
) %>%
  # Map amplification efficiencies and initial concentrations to taxa
  mutate(
    Amplification_Efficiency = case_when(
      Taxon == "A" ~ amplification_eff[1],
      Taxon == "B" ~ amplification_eff[2],
      Taxon == "C" ~ amplification_eff[3]
    ),
    Initial_Concentration = case_when(
      Taxon == "A" ~ initial_concentrations[1],
      Taxon == "B" ~ initial_concentrations[2],
      Taxon == "C" ~ initial_concentrations[3]
    )
  )

# Apply the formula to calculate DNA copies after each cycle
sim_data <- sim_data %>%
  mutate(
    DNA_Copies = Initial_Concentration * (1 + Amplification_Efficiency)^Cycle
  )

# Normalize proportions for plotting proportion plots
prop_data <- sim_data %>%
  group_by(Cycle) %>%
  mutate(Proportion = DNA_Copies / sum(DNA_Copies)) %>%
  ungroup()

# Plot 1: DNA copies across cycles for each taxon
ggplot(sim_data, aes(x = Cycle, y = DNA_Copies, color = Taxon)) +
  geom_line(size = 1.2) +
  scale_y_log10() +
  labs(
    title = "DNA Amplification Simulation",
    x = "PCR Cycle",
    y = "DNA Copies (log scale)",
    color = "Taxon"
  ) +
  theme_minimal()

# Plot 2: Proportions across cycles
ggplot(prop_data, aes(x = Cycle, y = Proportion, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  labs(
    title = "Proportion of DNA Copies Across PCR Cycles",
    x = "PCR Cycle",
    y = "Proportion",
    fill = "Taxon"
  ) +
  theme_minimal()

# Now simulate on my data -------------------------------------------------


# Load necessary libraries
library(tidyverse)
library(readxl)

# Load the data
file_path <- here("PCR_bias_correction/data/amp_effs/amp_effs_18s_processed.csv") # Update with actual path if needed
amp_data <- read.csv(file_path)

# Preview the data
head(amp_data)

# Define parameters
set.seed(42)  # For reproducibility

# Extract unique pools and size fractions
pools <- unique(amp_data$pool_type)
sizes <- unique(amp_data$size_fraction)

# Set arbitrary initial DNA concentrations (random but comparable for taxa)
initial_dna_conc <- runif(n = nrow(amp_data), min = 1e3, max = 1e4)

# Add initial concentrations to the data
amp_data <- amp_data %>%
  mutate(
    Initial_Concentration = initial_dna_conc
  )

# Remove 'clr_' prefix from Taxa names
amp_data <- amp_data %>%
  mutate(Taxa = str_remove(Taxa, "clr_"))

# Assign arbitrary initial concentrations
set.seed(42)
amp_data <- amp_data %>%
  mutate(Initial_Concentration = runif(n(), min = 1e3, max = 1e4))

# Corrected proportion calculation
sim_data <- amp_data %>%
  expand_grid(Cycle = cycles) %>%
  mutate(
    DNA_Copies = Initial_Concentration * (1 + amplification_efficiency)^Cycle
  ) %>%
  group_by(Cycle, pool_type, size_fraction) %>%
  mutate(Proportion = DNA_Copies / sum(DNA_Copies, na.rm = TRUE)) %>%
  ungroup()

# Define PCR cycles
cycles <- seq(0, 30, by = 10)

# Define custom colors for taxa
taxa_colors <- c(
  "other" = "grey",
  "Calanidae" = "#77DD77",       # Pastel green
  "Clausocalanidae" = "#72872d", # Bright pastel green
  "Eucalanidae" = "#89CFF0",     # Baby blue
  "Metridinidae" = "#9370DB",    # Pastel purple
  "Rhincalanidae" = "#4682B4",   # Steel blue
  "Paracalanidae" = "#B0E57C",   # Lime pastel green
  "Oithonidae" = "#f0ca62",      # Light blue
  "Euphausiidae" = "#FF6961",    # Pastel red
  "Salpidae" = "#FFB6C1",        # Pastel pink
  "unidentified Calanoida" = "#2d8087"   # Neutral for unidentified taxa
)

# Apply this palette in ggplot
plot <- sim_data %>%
  filter(pool_type == "AllandSub") %>%
  ggplot(aes(x = Cycle, y = Proportion, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack", width = 4) +
  facet_wrap(~ size_fraction, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = taxa_colors) +
  labs(
    x = "PCR Cycle",
    y = "Proportion",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Larger facet labels
    axis.text = element_text(size = 8),                  # Legible axis text
    axis.title = element_text(size = 10, face = "bold"), # Bold axis titles
    legend.text = element_text(size = 8),                # Smaller legend text
    legend.title = element_text(size = 9, face = "bold"),# Bold legend title
    legend.position = "bottom",                          # Place legend below plot
    legend.key.size = unit(0.5, "cm")                    # Compact legend keys
  )
plot

# Save the plot as a high-resolution PDF
ggsave(
  filename = "Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/amp_effs_simulation.pdf",
  plot = plot,
  device = "pdf",
  width = 112 / 25.4, # Convert 112 mm to inches
  height = 85 / 25.4, # Adjust the height proportionally for the layout
  dpi = 600
)

ggsave(
  filename = "Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/amp_effs_simulation.png",
  plot = plot,
  device = "png",
  width = 112 / 25.4, # Convert 112 mm to inches
  height = 85 / 25.4, # Adjust the height proportionally for the layout
  dpi = 600
)


# Now use real counts -----------------------------------------------------
#Predicted proportions
fido_s1_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s1_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)%>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() 


fido_s2_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s2_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)%>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() 


fido_s3_raw=read.csv(here("PCR_bias_correction/data/fido/phy/fido_18s_s3_ecdf_family_phy_all_subpools.csv")) %>% 
  select(-starts_with("X")) %>% 
  pivot_longer(cols = -Family, names_to = "Sample_ID", values_to = "n_reads") %>%
  mutate(Sample_ID_short= str_extract(Sample_ID, ".*(?=\\.[^.]+$)")) %>%
  group_by(Sample_ID_short, Family) %>%
  summarise(n_reads = sum(n_reads)) %>%
  filter(!grepl("All", Sample_ID_short)) %>% # Filter rows where Sample_ID_short doesn't contain "All"
  pivot_wider(names_from = Sample_ID_short, values_from = n_reads, values_fill = 0)%>%
  #11/2024 sum in Salpidae
  mutate(Family = ifelse(Family == "Salpidae", "other", Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup() 


merge(fido_s1_raw, fido_s2_raw, by = "Family", all = TRUE) %>%
  merge(.,fido_s3_raw, by = "Family", all = TRUE)%>%
  column_to_rownames("Family") %>%
  mutate(across(.cols = everything(), .fns = ~ coalesce(., 0)))-> fido_18s_merged_raw



#Metadata
#Filter to onshore offshore
clusters=read.csv(here("PCR_bias_correction/data/physical_environmental_data/pca_clusters.csv")) %>%
  mutate(Sample_ID=Sample_ID_dot) %>% 
  select(-Sample_ID_dot,-clust_group,-PC1) %>% unique()


env_metadata_phy=env_metadata%>% 
  left_join(.,clusters, by="Sample_ID") %>%
  column_to_rownames("Sample_ID_dot") 

#Make phyloseq objects
OTU = otu_table(as.matrix(fido_18s_merged_raw), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(zhan_taxa))
meta=sample_data(env_metadata_phy)


phy_18s_counts=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Taxa=asv_code) %>%
  select(-asv_code) %>% 
  group_by(Taxa, size_fraction, offshore_onshore) %>% 
  summarize(counts=median(n_reads))


# Normalize counts to proportions within each Cycle and Size Fraction
data_normalized <- phy_18s_counts %>%
  group_by(size_fraction, offshore_onshore) %>%
  mutate(Proportion = counts / sum(counts, na.rm = TRUE)) %>%
  ungroup()%>%
  left_join(amp_data %>% filter(pool_type=="AllandSub") %>% 
              select(-Initial_Concentration,-pool_type), by = c("Taxa", "size_fraction"))%>%
  mutate(counts_at_0 = counts / (1 + amplification_efficiency)^30)

# Define PCR cycles
cycles <- seq(0, 30, by = 10)

# Simulate DNA amplification across cycles
sim_data <- data_normalized %>%
  expand_grid(Cycle = cycles) %>%
  mutate(
    DNA_Copies = counts_at_0 * (1 + amplification_efficiency)^Cycle
  ) %>%
  group_by(offshore_onshore, size_fraction, Cycle) %>%
  mutate(Proportion = DNA_Copies / sum(DNA_Copies, na.rm = TRUE)) %>%
  ungroup()


# Function to create and save plots
create_and_save_plot <- function(data_subset, location) {
  plot <- data_subset %>%
    ggplot(aes(x = Cycle, y = Proportion, fill = Taxa)) +
    geom_bar(stat = "identity", position = "stack", width = 4) +
    facet_wrap(~ size_fraction, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = taxa_colors) +
    labs(
      x = "PCR Cycle",
      y = "Proportion",
      fill = "Taxa",
      title = paste("Stacked Bar Plot of DNA Proportions -", location)
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(10, 10, 20, 10)
    )
  # Display the plot in the plot window
  print(plot)
  
  # Save as PDF
  ggsave(
    filename = paste0("Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/amp_effs_simulation_", location, ".pdf"),
    plot = plot,
    device = "pdf",
    width = 112 / 25.4, # Convert 112 mm to inches
    height = 85 / 25.4, # Adjust the height proportionally for the layout
    dpi = 600
  )
  
  # Save as PNG
  ggsave(
    filename = paste0("Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/amp_effs_simulation_", location, ".png"),
    plot = plot,
    device = "png",
    width = 112 / 25.4, # Convert 112 mm to inches
    height = 85 / 25.4, # Adjust the height proportionally for the layout
    dpi = 600
  )
}

# Create and display plots for offshore and onshore data
create_and_save_plot(sim_data %>% filter(offshore_onshore == "offshore"), "offshore")
create_and_save_plot(sim_data %>% filter(offshore_onshore == "onshore"), "onshore")



# Next by cycle -----------------------------------------------------------

phy_18s_counts=phyloseq(OTU, TAX, meta)%>%
  phyloseq_transform_to_long(.) %>%
  mutate(Taxa=asv_code) %>%
  select(-asv_code) %>% 
  rename(Cycle_name=Cycle) %>% 
  group_by(Taxa, size_fraction, Cycle_name) %>% 
  summarize(counts=median(n_reads))


#Step 2: Normalize counts and join amplification efficiencies
data_normalized <- phy_18s_counts %>%
  group_by(size_fraction, Cycle_name) %>%
  mutate(Proportion = counts / sum(counts, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(
    amp_data %>%
      filter(pool_type == "AllandSub") %>% 
      select(-Initial_Concentration, -pool_type),
    by = c("Taxa", "size_fraction")
  ) %>%
  mutate(counts_at_0 = counts / (1 + amplification_efficiency)^30)

# Step 3: Define PCR cycles
cycles <- seq(0, 30, by = 10)

# Step 4: Simulate DNA amplification across cycles
sim_data <- data_normalized %>%
  expand_grid(Cycle = cycles) %>%
  mutate(
    DNA_Copies = counts_at_0 * (1 + amplification_efficiency)^Cycle
  ) %>%
  group_by(Cycle_name, size_fraction, Cycle) %>%
  mutate(Proportion = DNA_Copies / sum(DNA_Copies, na.rm = TRUE)) %>%
  ungroup()


# Step 6: Function to create, display, and save plots for each Cycle_name
create_and_save_plot <- function(data_subset, cycle_name) {
  plot <- data_subset %>%
    ggplot(aes(x = Cycle, y = Proportion, fill = Taxa)) +
    geom_bar(stat = "identity", position = "stack", width = 4) +
    facet_wrap(~ size_fraction, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = taxa_colors) +
    labs(
      x = "PCR Cycle",
      y = "Proportion",
      fill = "Taxa",
      title = paste("Cycle", cycle_name)
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(10, 10, 20, 10)
    )
  
  # Display the plot in the plot window
  print(plot)
  
  # Save as PDF
  ggsave(
    filename = paste0("Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/cycles/amp_effs_simulation_", cycle_name, ".pdf"),
    plot = plot,
    device = "pdf",
    width = 112 / 25.4, # Convert 112 mm to inches
    height = 85 / 25.4, # Adjust the height proportionally for the layout
    dpi = 600
  )
  
  # Save as PNG
  ggsave(
    filename = paste0("Q:/Dante/ZoopMetaB/New_analysis_plots/amp_effs_simulation/cycles/amp_effs_simulation_", cycle_name, ".png"),
    plot = plot,
    device = "png",
    width = 112 / 25.4, # Convert 112 mm to inches
    height = 85 / 25.4, # Adjust the height proportionally for the layout
    dpi = 600
  )
}

# Step 7: Loop through each unique Cycle_name and create/save plots
unique_cycle_names <- unique(sim_data$Cycle_name)

for (cycle_name in unique_cycle_names) {
  create_and_save_plot(sim_data %>% filter(Cycle_name == cycle_name), cycle_name)
}

