#Script to combine the blast and metazoogene taxa files 
librarian::shelf(tidyverse,here)

# 18s ---------------------------------------------------------------------

#Taxa Tables 
taxa_18s_meta=read.csv(here("data/phyloseq_bio_data/18S/metazoopruned18s_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_18s_blast=read.csv(here("data/taxa_files/zhang_taxa.csv")) %>%
  distinct(Hash, .keep_all = TRUE)


# ---- Keep only BLAST rows whose Hash is NOT in the meta file ----
blast_only <- anti_join(taxa_18s_blast, taxa_18s_meta, by = "Hash")

# ---- Combine: meta first (authoritative), then BLAST-only fills ----
final_taxa <- bind_rows(taxa_18s_meta, blast_only)

# ---- Order columns: desired taxonomy ranks first, then everything else ----
desired_order <- c("Phylum", "Subphylum", "Class", "Subclass", "Order",
                   "Superorder", "Family", "Genus", "Species")

# ensure any missing desired columns exist (to avoid select() errors)
missing_cols <- setdiff(desired_order, names(final_taxa))
if (length(missing_cols) > 0) {
  final_taxa[missing_cols] <- NA_character_
}

non_taxonomic_cols <- setdiff(names(final_taxa), desired_order)
new_order <- c(desired_order, non_taxonomic_cols)

final_taxa_18s <- final_taxa %>%
  select(all_of(new_order))

# ---- Inspect & write ----
print(final_taxa_18s)
write.csv(final_taxa_18s, here("data/taxa_files/blast_metazoo_18s.csv"), row.names = FALSE)

#Make version for realism vetting

# Keep only rows with classification to Family or lower
tax_family_plus <- final_taxa_18s %>%
  filter(!is.na(Family) & Family != "other") %>%
  select(Phylum:Species) %>%
  distinct()

# Write out to CSV
write.csv(
  tax_family_plus,
  here("data/taxa_files/final_taxa_18s_family_prevetting.csv"),
  row.names = FALSE
)

# Preview
head(tax_family_plus, 20)


# Post vetting, rejoin and derank taxa that are not ecologically feasible ------------------------------------------------------------
# ---------- Helpers ----------
is_blank <- function(x) is.na(x) | x == "" | tolower(x) == "other"
nz <- function(x) ifelse(is_blank(x), NA_character_, trimws(x))

# ---------- Load taxonomy (BLAST+Meta merged from the previous step) ----------
tax_path <- here("data/taxa_files/blast_metazoo_18s.csv")
final_taxa_18s <- read.csv(tax_path, check.names = FALSE) %>%
  mutate(
    across(c(Phylum, Subphylum, Class, Subclass, Order, Superorder, Family, Genus, Species), nz)
  )

# ---------- Load vetted file (Cleared column: 1 = OK, 0 = demote) ----------
vetted_path_proj <- here("data/taxa_files/final_taxa_18s_family_vetted.csv")

vetted <- read.csv(vetted_path_proj,
                   check.names = FALSE) %>%
  mutate(across(c(Phylum, Subphylum, Class, Subclass, Order, Superorder, Family, Genus, Species), nz)) %>%
  mutate(Cleared = as.integer(Cleared))

# ---------- Identify taxa to demote by their HIGHEST assigned rank ----------
# Rules: pick highest among Species > Genus > Family (and we allow Family -> Order if Order exists)
v0 <- vetted %>% filter(Cleared == 0)

species_to_demote <- v0 %>%
  filter(!is_blank(Species)) %>%
  pull(Species) %>% unique()

genus_to_demote <- v0 %>%
  filter(is_blank(Species), !is_blank(Genus)) %>%
  pull(Genus) %>% unique()

family_rows <- v0 %>%
  filter(is_blank(Species), is_blank(Genus), !is_blank(Family))

family_to_demote <- family_rows %>% pull(Family) %>% unique()

# (Optional) We’ll keep the target Orders for Family demotions when available in vetted
# If vetted lacks Order for those families, we still demote Family -> (up) Order by simply removing Family.
# No change to Order text itself is needed; we just remove the more specific rank.
# Record families that actually have an Order in the vetted file (for auditing only)
families_with_order <- family_rows %>% filter(!is_blank(Order)) %>% pull(Family) %>% unique()

# ---------- Apply rank-by-rank demotion to ALL matching rows in final_taxa_18s ----------
# Demote Species -> Genus
demoted <- final_taxa_18s %>%
  mutate(
    demote_species_flag = !is_blank(Species) & Species %in% species_to_demote,
    Species = if_else(demote_species_flag, NA_character_, Species)
  ) %>%
  # Demote Genus -> Family (only where current highest rank is Genus)
  mutate(
    is_current_genus = is_blank(Species) & !is_blank(Genus),
    demote_genus_flag = is_current_genus & Genus %in% genus_to_demote,
    Genus = if_else(demote_genus_flag, NA_character_, Genus)
  ) %>%
  # Demote Family -> Order (only where current highest rank is Family)
  mutate(
    is_current_family = is_blank(Species) & is_blank(Genus) & !is_blank(Family),
    demote_family_flag = is_current_family & Family %in% family_to_demote,
    Family = if_else(demote_family_flag, NA_character_, Family)
  )

# ---------- (Optional) Create an audit table of what changed ----------
audit <- final_taxa_18s %>%
  transmute(
    Hash = if ("Hash" %in% names(final_taxa_18s)) Hash else NA_character_,
    old_top_rank = case_when(
      !is_blank(Species) ~ "Species",
      !is_blank(Genus)   ~ "Genus",
      !is_blank(Family)  ~ "Family",
      !is_blank(Order)   ~ "Order",
      TRUE               ~ NA_character_
    ),
    old_Family = Family, old_Genus = Genus, old_Species = Species
  ) %>%
  bind_cols(
    demoted %>%
      transmute(
        new_top_rank = case_when(
          !is_blank(Species) ~ "Species",
          !is_blank(Genus)   ~ "Genus",
          !is_blank(Family)  ~ "Family",
          !is_blank(Order)   ~ "Order",
          TRUE               ~ NA_character_
        ),
        new_Family = Family, new_Genus = Genus, new_Species = Species
      )
  ) %>%
  mutate(changed = old_top_rank != new_top_rank)

# Quick counts:
table(audit$changed, useNA = "ifany")
with(demoted, c(
  species_demoted = sum(demote_species_flag, na.rm = TRUE),
  genus_demoted   = sum(demote_genus_flag,   na.rm = TRUE),
  family_demoted  = sum(demote_family_flag,  na.rm = TRUE)
))

# ---------- Save updated taxonomy ----------
final_taxa_18s_qc = final_taxa_18s %>% 
  select(Phylum:Kingdom)
out_path <- here("data/taxa_files/blast_metazoo_18s_qc.csv")
write.csv(final_taxa_18s_qc, out_path, row.names = FALSE)




# COI ---------------------------------------------------------------------

#Taxa Tables 
taxa_coi_meta=read.csv(here("data/phyloseq_bio_data/COI/metazooprunedcoi_tax.csv"))%>%
  mutate(non_na_count = rowSums(!is.na(select(., -Hash)))) %>%
  group_by(Hash) %>%
  filter(rank(desc(non_na_count)) == 1) %>%
  select(-non_na_count) %>%
  ungroup() 

#BlAST
taxa_coi_blast=read.csv(here("data/taxa_files/leray_taxa.csv")) %>%
  distinct(Hash, .keep_all = TRUE)

# ---- Keep only BLAST rows whose Hash is NOT in the meta file ----
blast_only <- anti_join(taxa_coi_blast, taxa_coi_meta, by = "Hash")

# ---- Combine: meta first (authoritative), then BLAST-only fills ----
final_taxa <- bind_rows(taxa_coi_meta, blast_only)

# ---- Order columns: desired taxonomy ranks first, then everything else ----
desired_order <- c("Phylum", "Subphylum", "Class", "Subclass", "Order",
                   "Superorder", "Family", "Genus", "Species")

# ensure any missing desired columns exist (to avoid select() errors)
missing_cols <- setdiff(desired_order, names(final_taxa))
if (length(missing_cols) > 0) {
  final_taxa[missing_cols] <- NA_character_
}

non_taxonomic_cols <- setdiff(names(final_taxa), desired_order)
new_order <- c(desired_order, non_taxonomic_cols)

final_taxa_coi <- final_taxa %>%
  select(all_of(new_order))

# ---- Inspect & write ----
print(final_taxa_coi)
write.csv(final_taxa_coi, here("data/taxa_files/blast_metazoo_coi.csv"), row.names = FALSE)

#Make version for realism vetting

# Keep only rows with classification to Family or lower
tax_family_plus <- final_taxa_coi %>%
  filter(!is.na(Family) & Family != "other") %>%
  select(Phylum:Species) %>%
  arrange(Species, desc(!is.na(Superorder))) %>%  # rows with Superorder first
  distinct(Species, .keep_all = TRUE)

# Write out to CSV
write.csv(
  tax_family_plus,
  here("data/taxa_files/final_taxa_coi_family_prevetting.csv"),
  row.names = FALSE
)

# Preview
head(tax_family_plus, 20)


# Post vetting, rejoin and derank taxa that are not ecologically feasible ------------------------------------------------------------

# ---------- Load taxonomy (BLAST+Meta merged from the previous step) ----------
tax_path <- here("data/taxa_files/blast_metazoo_coi.csv")
final_taxa_coi <- read.csv(tax_path, check.names = FALSE) %>%
  mutate(
    across(c(Phylum, Subphylum, Class, Subclass, Order, Superorder, Family, Genus, Species), nz)
  )

# ---------- Load vetted file (Cleared column: 1 = OK, 0 = demote) ----------
vetted_path_proj <- here("data/taxa_files/final_taxa_coi_family_vetted.csv")

vetted <- read.csv(vetted_path_proj,
                   check.names = FALSE) %>%
  mutate(across(c(Phylum, Subphylum, Class, Subclass, Order, Superorder, Family, Genus, Species), nz)) %>%
  mutate(Cleared = as.integer(Cleared))

# ---------- Identify taxa to demote by their HIGHEST assigned rank ----------
# Rules: pick highest among Species > Genus > Family (and we allow Family -> Order if Order exists)
v0 <- vetted %>% filter(Cleared == 0)

species_to_demote <- v0 %>%
  filter(!is_blank(Species)) %>%
  pull(Species) %>% unique()

genus_to_demote <- v0 %>%
  filter(is_blank(Species), !is_blank(Genus)) %>%
  pull(Genus) %>% unique()

family_rows <- v0 %>%
  filter(is_blank(Species), is_blank(Genus), !is_blank(Family))

family_to_demote <- family_rows %>% pull(Family) %>% unique()

# (Optional) We’ll keep the target Orders for Family demotions when available in vetted
# If vetted lacks Order for those families, we still demote Family -> (up) Order by simply removing Family.
# No change to Order text itself is needed; we just remove the more specific rank.
# Record families that actually have an Order in the vetted file (for auditing only)
families_with_order <- family_rows %>% filter(!is_blank(Order)) %>% pull(Family) %>% unique()

# ---------- Apply rank-by-rank demotion to ALL matching rows in final_taxa_coi ----------
# Demote Species -> Genus
demoted <- final_taxa_coi %>%
  mutate(
    demote_species_flag = !is_blank(Species) & Species %in% species_to_demote,
    Species = if_else(demote_species_flag, NA_character_, Species)
  ) %>%
  # Demote Genus -> Family (only where current highest rank is Genus)
  mutate(
    is_current_genus = is_blank(Species) & !is_blank(Genus),
    demote_genus_flag = is_current_genus & Genus %in% genus_to_demote,
    Genus = if_else(demote_genus_flag, NA_character_, Genus)
  ) %>%
  # Demote Family -> Order (only where current highest rank is Family)
  mutate(
    is_current_family = is_blank(Species) & is_blank(Genus) & !is_blank(Family),
    demote_family_flag = is_current_family & Family %in% family_to_demote,
    Family = if_else(demote_family_flag, NA_character_, Family)
  )

# ---------- (Optional) Create an audit table of what changed ----------
audit <- final_taxa_coi %>%
  transmute(
    Hash = if ("Hash" %in% names(final_taxa_coi)) Hash else NA_character_,
    old_top_rank = case_when(
      !is_blank(Species) ~ "Species",
      !is_blank(Genus)   ~ "Genus",
      !is_blank(Family)  ~ "Family",
      !is_blank(Order)   ~ "Order",
      TRUE               ~ NA_character_
    ),
    old_Family = Family, old_Genus = Genus, old_Species = Species
  ) %>%
  bind_cols(
    demoted %>%
      transmute(
        new_top_rank = case_when(
          !is_blank(Species) ~ "Species",
          !is_blank(Genus)   ~ "Genus",
          !is_blank(Family)  ~ "Family",
          !is_blank(Order)   ~ "Order",
          TRUE               ~ NA_character_
        ),
        new_Family = Family, new_Genus = Genus, new_Species = Species
      )
  ) %>%
  mutate(changed = old_top_rank != new_top_rank)

# Quick counts:
table(audit$changed, useNA = "ifany")
with(demoted, c(
  species_demoted = sum(demote_species_flag, na.rm = TRUE),
  genus_demoted   = sum(demote_genus_flag,   na.rm = TRUE),
  family_demoted  = sum(demote_family_flag,  na.rm = TRUE)
))

# ---------- Save updated taxonomy ----------
final_taxa_coi_qc = final_taxa_coi %>% 
  select(Phylum:Kingdom)
out_path <- here("data/taxa_files/blast_metazoo_coi_qc.csv")
write.csv(final_taxa_coi_qc, out_path, row.names = FALSE)


