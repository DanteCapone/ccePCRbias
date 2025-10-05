#Script to combine the blast and metazoogene taxa files 


# 18s ---------------------------------------------------------------------

#Taxa Tables 
taxa_18s_meta=read.csv(here("data/past/metazoopruned18s_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_18s_blast=read.csv(here("data/raw_data/BLAST_taxa_class/zhang_taxa.csv")) %>%
  distinct(Hash, .keep_all = TRUE)


# Merge the two taxonomy tables
merged_taxa <- taxa_18s_meta %>%
  full_join(taxa_18s_blast, by = "Hash", suffix = c(".meta", ".blast")) %>% 
  select(-Kingdom)

# Define the taxonomy ranks
taxonomy_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Adjust the function to compare two versions and select the better resolution
select_higher_resolution <- function(meta, blast) {
  non_na_meta <- sum(!is.na(meta), na.rm = TRUE)
  non_na_blast <- sum(!is.na(blast), na.rm = TRUE)
  
  if (non_na_blast > non_na_meta) {
    return(blast)
  } else {
    return(meta)  # Favors meta in case of tie
  }
}

# Apply the function across each taxonomy rank within a rowwise framework
merged_taxa <- merged_taxa %>%
  rowwise() %>%
  mutate(across(all_of(paste0(taxonomy_ranks, ".meta")), 
                .fns = ~ select_higher_resolution(.x, get(gsub("meta", "blast", cur_column()))),
                .names = "{gsub('\\\\.meta$', '', .col)}")) %>%
  ungroup()

# Optional: clean up by removing the original columns with suffixes
final_taxa <- merged_taxa %>%
  select(-matches("\\.meta$|\\.blast$"))

# Define the complete desired order of taxonomy ranks including new ranks
desired_order <- c("Phylum", "Subphylum", "Class", "Subclass", "Order", "Superorder", "Family", "Genus", "Species")

# Create a vector of column names from 'final_taxa' that are not in 'desired_order'
non_taxonomic_cols <- setdiff(names(final_taxa), desired_order)

# Combine the ordered taxonomic ranks with the non-taxonomic columns
# This keeps the non-taxonomic columns at the end, or you can interweave them as needed
new_order <- c(desired_order, non_taxonomic_cols)

# Reorder the columns of 'final_taxa' according to 'new_order'
final_taxa_18s <- final_taxa %>%
  select(all_of(new_order))

# Check the reordered output
print(final_taxa_18s)


write.csv(final_taxa_18s,here("data/taxa_files/blast_metazoo_18s.csv"))



# COI ---------------------------------------------------------------------

#Taxa Tables 
taxa_coi_meta=read.csv(here("data/past/metazooprunedcoi_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_coi_blast=read.csv(here("data/raw_data/BLAST_taxa_class/leray_taxa.csv")) %>%
  distinct(Hash, .keep_all = TRUE)


# Merge the two taxonomy tables
merged_taxa <- taxa_coi_meta %>%
  full_join(taxa_coi_blast, by = "Hash", suffix = c(".meta", ".blast")) %>% 
  select(-Kingdom)

# Define the taxonomy ranks
taxonomy_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Adjust the function to compare two versions and select the better resolution
select_higher_resolution <- function(meta, blast) {
  non_na_meta <- sum(!is.na(meta), na.rm = TRUE)
  non_na_blast <- sum(!is.na(blast), na.rm = TRUE)
  
  if (non_na_blast > non_na_meta) {
    return(blast)
  } else {
    return(meta)  # Favors meta in case of tie
  }
}

# Apply the function across each taxonomy rank within a rowwise framework
merged_taxa <- merged_taxa %>%
  rowwise() %>%
  mutate(across(all_of(paste0(taxonomy_ranks, ".meta")), 
                .fns = ~ select_higher_resolution(.x, get(gsub("meta", "blast", cur_column()))),
                .names = "{gsub('\\\\.meta$', '', .col)}")) %>%
  ungroup()

# Optional: clean up by removing the original columns with suffixes
final_taxa <- merged_taxa %>%
  select(-matches("\\.meta$|\\.blast$"))

# Define the complete desired order of taxonomy ranks including new ranks
desired_order <- c("Phylum", "Subphylum", "Class", "Subclass", "Order", "Superorder", "Family", "Genus", "Species")

# Create a vector of column names from 'final_taxa' that are not in 'desired_order'
non_taxonomic_cols <- setdiff(names(final_taxa), desired_order)

# Combine the ordered taxonomic ranks with the non-taxonomic columns
# This keeps the non-taxonomic columns at the end, or you can interweave them as needed
new_order <- c(desired_order, non_taxonomic_cols)

# Reorder the columns of 'final_taxa' according to 'new_order'
final_taxa_coi <- final_taxa %>%
  select(all_of(new_order))

# Check the reordered output
print(final_taxa_coi)


write.csv(final_taxa_coi,here("data/taxa_files/blast_metazoo_coi.csv"))
