# =============================================================================
# publication_pcr_bias_fido_workflow_cycle10.R
# -----------------------------------------------------------------------------
# PUBLICATION-READY PCR BIAS MITIGATION WORKFLOW (conservative cycle-10 variant)
#
# WHAT THIS SCRIPT DOES
#   Reproduces the Bayesian PCR-bias correction used in this project (the `fido`
#   multinomial logistic-normal "pibble" model) and applies it to the 18S
#   metabarcoding data for three zooplankton size fractions. For every field
#   sample it predicts the corrected relative abundances at a CONFIGURABLE PCR
#   cycle. The published main analysis extrapolated to PCR cycle 0; this script
#   instead extrapolates to PCR cycle 10 as a conservative sensitivity analysis.
#
# WHY CYCLE 10 (the reviewer's concern)
#   The calibration dilution series spans observed PCR cycles ~20-28 and field
#   samples were sequenced at cycle 30. Extrapolating all the way back to cycle 0
#   is the strongest possible extrapolation and lies far outside the observed
#   range, which a reviewer flagged as potentially too extreme. Predicting at
#   cycle 10 still removes most of the later-cycle amplification distortion (the
#   bias accumulates with each cycle) but stays much closer to the observed data,
#   so it is a deliberately CONSERVATIVE estimate. The goal is NOT to claim cycle
#   10 is the "true" starting community; it is a robustness check that reduces
#   amplification bias while avoiding the most aggressive extrapolation.
#
# INPUTS (read-only; already produced by the existing pre-processing pipeline)
#   data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv   (taxa x samples, Family)
#   data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv
#   data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv
#   data/fido/meta_18s_unaveraged_all.csv                   (sample metadata)
#   (optional COI, exploratory) data/fido/phy/fido_coi_s*_ecdf_genus_phy_all_subpools.csv
#                               data/fido/meta_coi_unaveraged_all.csv
#
# OUTPUTS (all written to a NEW directory; existing files are never overwritten)
#   publication_pcr_bias_cycle10_outputs/
#     corrected_relabund_cycle10_18s_long.csv       corrected + raw comp (long)
#     corrected_relabund_cycle10_18s_wide.csv       corrected mean comp (wide)
#     amplification_efficiency_18s.csv              per-taxon cycle_num Lambda
#     filtering_summary_18s.csv                     taxa/sample counts per fraction
#     calibration_pool_coverage_18s.csv             cycles available per pool
#     taxon_calibration_coverage_18s.csv            taxon presence across cycles
#     model_matrix_dimensions_18s.csv               design-matrix sanity checks
#     simplex_check_cycle10_18s.csv                 predictions sum to 1 check
#     value_sanity_check_cycle10_18s.csv            neg/NA/Inf/extreme check
#     fig_raw_vs_cycle10_<fraction>_18s.png         raw vs corrected composition
#     comparison_raw_cycle0_cycle10_18s.csv         (if a cycle-0 output exists)
#     sessionInfo.txt                               package versions for repro
#     README_cycle10_sensitivity.md                 plain-language explanation
#
# HOW TO RUN
#   From the project root (recommended):
#     "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" \
#       scripts/pre-processing/pcr_bias_mitigation_prep_fido/publication_pcr_bias_fido_workflow_cycle10.R
#   or interactively: open in RStudio at the project root and source the file.
#
# This script consults but never modifies the original exploratory scripts:
#   predict_original_proportions_18s_phy_all_and_sub_pools.R (model + cycle-0 logic)
#   predict_original_proportions_18s_phy.R                  (per-fraction loop)
#   fido_amp_effs.R                                         (amplification effs)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fido)
  library(here)
})

# Resolve the helper functions relative to this script's folder.
source(here::here("scripts", "pre-processing", "pcr_bias_mitigation_prep_fido",
                  "publication_pcr_bias_fido_functions.R"))


# -----------------------------------------------------------------------------
# CONFIGURATION  (everything a reviewer might want to change lives here)
# -----------------------------------------------------------------------------

# *** The single key parameter of this sensitivity analysis. ***
# Set to 0 to reproduce the original main analysis; set to 10 for the
# conservative variant requested in revision.
prediction_cycle <- 10

# PCR cycle at which field samples were actually sequenced (used as the observed
# baseline for the raw-vs-corrected comparison).
observed_cycle <- 30

# Number of posterior samples for the production fits (matches the originals).
n_posterior_samples <- 10000

# Whether to re-run the Gamma (prior-scale) log-marginal-likelihood grid search.
# OFF by default because it refits the model ~20x per fraction; the values below
# were selected this way in the original all_and_subpools script.
run_gamma_selection <- FALSE

# Gamma (prior scale on regression coefficients) chosen per size fraction,
# carried over from predict_original_proportions_18s_phy_all_and_sub_pools.R.
gamma_by_fraction <- c(s1 = 10, s2 = 20, s3 = 20)

# Calibration cycles that every taxon should ideally span (18S dilution series).
required_calibration_cycles_18s <- c(20, 24, 28)

# Run the exploratory COI marker too? COI is kept exploratory (Task 6) because
# its calibration series covers fewer PCR cycles (24/28 only, no cycle 20),
# giving a weaker amplification-slope estimate. 18S is the robust main analysis.
run_coi_exploratory <- FALSE

# 18S size fractions: input table, size label, taxon rank column, Gamma.
fractions_18s <- list(
  s1 = list(counts = here("data/fido/phy/fido_18s_s1_family_phy_all_subpools.csv"),
            label = "0.2-0.5mm"),
  s2 = list(counts = here("data/fido/phy/fido_18s_s2_family_phy_all_subpools.csv"),
            label = "0.5-1mm"),
  s3 = list(counts = here("data/fido/phy/fido_18s_s3_family_phy_all_subpools.csv"),
            label = "1-2mm")
)
meta_path_18s <- here("data/fido/meta_18s_unaveraged_all.csv")

# Fold Salpidae into "other" (11/2024 decision in the all_and_subpools script):
# gelatinous salp reads behave erratically across PCR cycles and destabilise the
# logistic-normal fit, so they are merged into the catch-all category.
collapse_18s <- c(Salpidae = "other")

# Create the dedicated output directory (never clobbers an existing one).
out_dir <- pcr_versioned_path(here("publication_pcr_bias_cycle10_outputs"))
pcr_ensure_dir(out_dir)
message("Writing outputs to: ", out_dir)


# -----------------------------------------------------------------------------
# MAIN 18S WORKFLOW  (loop over the three size fractions)
# -----------------------------------------------------------------------------

predictions_18s   <- list()
amp_eff_18s       <- list()
filtering_18s     <- list()
cycle_cov_18s     <- list()
taxon_cov_18s     <- list()
matrix_dims_18s   <- list()
gamma_search_18s  <- list()

for (frac in names(fractions_18s)) {
  cfg <- fractions_18s[[frac]]
  message("\n=== 18S ", frac, " (", cfg$label, ") ===")

  # 1. Load + agglomerate counts, load + align metadata.
  Y_raw  <- pcr_load_counts(cfg$counts, taxon_col = "Family", collapse_rare = collapse_18s)
  aligned <- pcr_load_and_align(Y_raw, meta_path_18s)
  meta <- aligned$meta
  Y    <- aligned$counts
  if (length(aligned$dropped_samples) > 0) {
    message("  Note: ", length(aligned$dropped_samples),
            " count column(s) had no metadata match and were dropped: ",
            paste(aligned$dropped_samples, collapse = ", "))
  }

  # 2. Design matrix (shared cycle slope + per-sample intercepts).
  X <- pcr_build_design(meta)

  # 2b. Diagnostic: model-matrix dimensions and Y/X sample alignment.
  matrix_dims_18s[[frac]] <- pcr_diag_model_matrix(Y, X, cfg$label, "18S")

  # 3. Optional Gamma selection by log marginal likelihood.
  gamma <- unname(gamma_by_fraction[[frac]])
  if (run_gamma_selection) {
    gs <- pcr_select_gamma(Y, X)
    gamma_search_18s[[frac]] <- dplyr::mutate(gs, size = cfg$label, marker = "18S")
    gamma <- gs$gamma[which.max(gs$logML)]
    message("  Gamma selected by logML: ", gamma)
  }

  # 4. Fit the pibble model (CLR + proportion representations).
  fit <- pcr_fit_pibble(Y, X, gamma = gamma, n_samples = n_posterior_samples)

  # 5. Per-taxon amplification efficiency (cycle_num Lambda, CLR coords).
  amp_eff_18s[[frac]] <- pcr_amplification_efficiency(fit$fit_clr, cfg$label, "18S")

  # 6. Predict corrected relative abundances at the configured PCR cycle.
  preds <- pcr_predict_at_cycle(fit$fit_prop, X, Y,
                                prediction_cycle = prediction_cycle,
                                size_label = cfg$label,
                                observed_cycle = observed_cycle)
  predictions_18s[[frac]] <- preds

  # 7. Diagnostics that depend on the data/predictions.
  filtering_18s[[frac]] <- tibble::tibble(
    marker = "18S", size = cfg$label,
    n_taxa = nrow(Y), n_samples = ncol(Y),
    n_calibration_samples = sum(meta$sample_num == "Calibration"),
    n_field_samples = sum(meta$sample_num != "Calibration"),
    gamma_used = gamma
  )
  cycle_cov_18s[[frac]] <- pcr_diag_cycle_coverage(meta, cfg$label, "18S")
  taxon_cov_18s[[frac]] <- pcr_diag_taxon_calibration(
    Y, meta, required_calibration_cycles_18s, cfg$label, "18S")
}

# Bind per-fraction results.
predictions_18s_df <- dplyr::bind_rows(predictions_18s)
amp_eff_18s_df     <- dplyr::bind_rows(amp_eff_18s)
filtering_18s_df   <- dplyr::bind_rows(filtering_18s)
cycle_cov_18s_df   <- dplyr::bind_rows(cycle_cov_18s)
taxon_cov_18s_df   <- dplyr::bind_rows(taxon_cov_18s)
matrix_dims_18s_df <- dplyr::bind_rows(matrix_dims_18s)

# Post-prediction sanity checks (Task 7).
simplex_check_18s <- pcr_diag_simplex(predictions_18s_df)
value_check_18s   <- pcr_diag_values(predictions_18s_df)

# Loud warnings if any prediction fails a sanity check.
if (any(!simplex_check_18s$sums_to_one)) {
  warning("Some predicted compositions do not sum to 1; inspect simplex_check_cycle10_18s.csv")
}
if (value_check_18s$n_negative > 0 || value_check_18s$n_missing > 0 ||
    value_check_18s$n_infinite > 0 || value_check_18s$n_gt_one > 0) {
  warning("Predicted compositions contain negative/missing/infinite/>1 values; inspect value_sanity_check_cycle10_18s.csv")
}


# -----------------------------------------------------------------------------
# COMPARISON: raw reads vs cycle-0 (if available) vs cycle-10 (Task 5)
# -----------------------------------------------------------------------------
# We do NOT refit a cycle-0 model here. Instead we reuse the most recent existing
# cycle-0 output (predicted_og_18s_*_s*_phy_all_and_subpools.csv) if present, so
# the comparison costs no extra model runs. If no cycle-0 output exists the
# comparison gracefully falls back to raw-vs-cycle10 only.

# Mean corrected composition per taxon/fraction at the prediction cycle.
cycle10_mean <- predictions_18s_df %>%
  dplyr::filter(grepl("^predicted ", .data$replicate)) %>%
  dplyr::group_by(.data$size, coord = .data$coord) %>%
  dplyr::summarise(mean_prop = mean(.data$n_reads, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(scenario = paste0("corrected_cycle", prediction_cycle))

# Mean raw (observed cycle-30) composition per taxon/fraction.
raw_mean <- predictions_18s_df %>%
  dplyr::filter(!grepl("^predicted ", .data$replicate)) %>%
  dplyr::group_by(.data$size, coord = .data$coord) %>%
  dplyr::summarise(mean_prop = mean(.data$n_reads, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(scenario = "raw_cycle30")

comparison_18s <- dplyr::bind_rows(raw_mean, cycle10_mean)

# Try to fold in an existing cycle-0 output for a three-way comparison.
find_latest_cycle0 <- function(frac_tag) {
  dir <- here("data/predicted_og")
  if (!dir.exists(dir)) return(NA_character_)
  files <- list.files(
    dir,
    pattern = paste0("^predicted_og_18s_.*_", frac_tag, "_phy_all_and_subpools\\.csv$"),
    full.names = TRUE)
  files <- files[!grepl("clr|nocollodaria", files)]
  if (length(files) == 0) return(NA_character_)
  dates <- as.Date(stringr::str_extract(basename(files), "[0-9]{2}_[0-9]{2}_[0-9]{4}"),
                   format = "%m_%d_%Y")
  files[which.max(dates)]
}

cycle0_files <- c(s1 = find_latest_cycle0("s1"),
                  s2 = find_latest_cycle0("s2"),
                  s3 = find_latest_cycle0("s3"))

if (all(!is.na(cycle0_files))) {
  message("Including existing cycle-0 outputs in the comparison:\n  ",
          paste(basename(cycle0_files), collapse = "\n  "))
  cycle0_mean <- purrr::map_dfr(cycle0_files, function(f) {
    readr::read_csv(f, show_col_types = FALSE) %>%
      dplyr::filter(grepl("^predicted ", .data$replicate)) %>%
      dplyr::group_by(.data$size, coord = .data$coord) %>%
      dplyr::summarise(mean_prop = mean(.data$n_reads, na.rm = TRUE), .groups = "drop")
  }) %>%
    dplyr::mutate(scenario = "corrected_cycle0")
  comparison_18s <- dplyr::bind_rows(comparison_18s, cycle0_mean)
} else {
  message("No complete set of cycle-0 outputs found; comparison uses raw vs cycle-",
          prediction_cycle, " only.")
}


# -----------------------------------------------------------------------------
# FIGURES: raw vs corrected composition per fraction (Task 8, optional figures)
# -----------------------------------------------------------------------------
# Show the top taxa so the figure stays legible; bias is most visible in the
# dominant families. Bars compare observed raw cycle-30 against corrected cycle-10.
make_comp_figure <- function(frac_label, top_n = 12) {
  d <- comparison_18s %>%
    dplyr::filter(.data$size == frac_label,
                  .data$scenario %in% c("raw_cycle30",
                                        paste0("corrected_cycle", prediction_cycle)))
  top_taxa <- d %>% dplyr::group_by(.data$coord) %>%
    dplyr::summarise(m = max(.data$mean_prop), .groups = "drop") %>%
    dplyr::slice_max(.data$m, n = top_n) %>% dplyr::pull(.data$coord)
  ggplot(dplyr::filter(d, .data$coord %in% top_taxa),
         aes(x = reorder(.data$coord, .data$mean_prop),
             y = .data$mean_prop, fill = .data$scenario)) +
    geom_col(position = position_dodge(width = 0.8)) +
    coord_flip() +
    labs(title = paste0("18S ", frac_label,
                        ": raw (cycle 30) vs corrected (cycle ", prediction_cycle, ")"),
         x = "Family", y = "Mean relative abundance", fill = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
}

for (frac in names(fractions_18s)) {
  lab <- fractions_18s[[frac]]$label
  fig <- make_comp_figure(lab)
  ggsave(filename = pcr_versioned_path(file.path(
            out_dir, paste0("fig_raw_vs_cycle", prediction_cycle, "_", frac, "_18s.png"))),
         plot = fig, width = 8, height = 6, dpi = 300)
}


# -----------------------------------------------------------------------------
# OPTIONAL EXPLORATORY COI ANALYSIS (Task 6: clearly marked exploratory)
# -----------------------------------------------------------------------------
predictions_coi_df <- NULL
if (run_coi_exploratory) {
  message("\n### EXPLORATORY COI ANALYSIS (weaker calibration coverage) ###")
  fractions_coi <- list(
    s1 = list(counts = here("data/fido/phy/fido_coi_s1_ecdf_genus_phy_all_subpools.csv"),
              label = "0.2-0.5mm"),
    s2 = list(counts = here("data/fido/phy/fido_coi_s2_ecdf_genus_phy_all_subpools.csv"),
              label = "0.5-1mm"),
    s3 = list(counts = here("data/fido/phy/fido_coi_s3_ecdf_genus_phy_all_subpools.csv"),
              label = "1-2mm")
  )
  meta_path_coi <- here("data/fido/meta_coi_unaveraged_all.csv")
  coi_preds <- list()
  for (frac in names(fractions_coi)) {
    cfg <- fractions_coi[[frac]]
    Yc <- pcr_load_counts(cfg$counts, taxon_col = "Genus")  # COI rank = Genus
    al <- pcr_load_and_align(Yc, meta_path_coi)
    Xc <- pcr_build_design(al$meta)
    fitc <- pcr_fit_pibble(al$counts, Xc, gamma = 20, n_samples = n_posterior_samples)
    coi_preds[[frac]] <- pcr_predict_at_cycle(
      fitc$fit_prop, Xc, al$counts,
      prediction_cycle = prediction_cycle, size_label = cfg$label,
      observed_cycle = observed_cycle) %>%
      dplyr::mutate(marker = "COI_exploratory")
  }
  predictions_coi_df <- dplyr::bind_rows(coi_preds)
}


# -----------------------------------------------------------------------------
# EXPORTS  (never overwrite; pcr_versioned_path() guards every write)
# -----------------------------------------------------------------------------
write_csv_safe <- function(df, name) {
  readr::write_csv(df, pcr_versioned_path(file.path(out_dir, name)))
}

write_csv_safe(predictions_18s_df,
               paste0("corrected_relabund_cycle", prediction_cycle, "_18s_long.csv"))

# Wide table of corrected mean composition (taxa x fraction) for manuscript use.
corrected_wide_18s <- cycle10_mean %>%
  dplyr::select("size", "coord", "mean_prop") %>%
  tidyr::pivot_wider(names_from = "size", values_from = "mean_prop")
write_csv_safe(corrected_wide_18s,
               paste0("corrected_relabund_cycle", prediction_cycle, "_18s_wide.csv"))

write_csv_safe(amp_eff_18s_df,     "amplification_efficiency_18s.csv")
write_csv_safe(filtering_18s_df,   "filtering_summary_18s.csv")
write_csv_safe(cycle_cov_18s_df,   "calibration_pool_coverage_18s.csv")
write_csv_safe(taxon_cov_18s_df,   "taxon_calibration_coverage_18s.csv")
write_csv_safe(matrix_dims_18s_df, "model_matrix_dimensions_18s.csv")
write_csv_safe(simplex_check_18s,
               paste0("simplex_check_cycle", prediction_cycle, "_18s.csv"))
write_csv_safe(value_check_18s,
               paste0("value_sanity_check_cycle", prediction_cycle, "_18s.csv"))
write_csv_safe(comparison_18s, "comparison_raw_cycle0_cycle10_18s.csv")
if (run_gamma_selection && length(gamma_search_18s) > 0) {
  write_csv_safe(dplyr::bind_rows(gamma_search_18s), "gamma_selection_logML_18s.csv")
}
if (!is.null(predictions_coi_df)) {
  write_csv_safe(predictions_coi_df,
                 paste0("EXPLORATORY_corrected_relabund_cycle", prediction_cycle, "_coi_long.csv"))
}

# Record exact package versions for reproducibility.
writeLines(capture.output(sessionInfo()),
           pcr_versioned_path(file.path(out_dir, "sessionInfo.txt")))

# Short README explaining the cycle-10 sensitivity analysis.
readme_lines <- c(
  "# PCR bias mitigation - conservative cycle-10 sensitivity analysis",
  "",
  "## What this is",
  "Corrected 18S relative abundances produced by the `fido` multinomial",
  "logistic-normal (pibble) model, extrapolated to **PCR cycle 10** instead of",
  "PCR cycle 0.",
  "",
  "## Inputs expected",
  "- data/fido/phy/fido_18s_s{1,2,3}_family_phy_all_subpools.csv",
  "- data/fido/meta_18s_unaveraged_all.csv",
  "",
  "## Outputs created (this folder)",
  "- corrected_relabund_cycle10_18s_long.csv / _wide.csv : corrected compositions",
  "- amplification_efficiency_18s.csv : per-taxon cycle_num Lambda (CLR)",
  "- filtering_summary_18s.csv, calibration_pool_coverage_18s.csv,",
  "  taxon_calibration_coverage_18s.csv, model_matrix_dimensions_18s.csv : diagnostics",
  "- simplex_check_cycle10_18s.csv, value_sanity_check_cycle10_18s.csv : QC",
  "- comparison_raw_cycle0_cycle10_18s.csv : raw vs cycle-0 (if available) vs cycle-10",
  "- fig_raw_vs_cycle10_s{1,2,3}_18s.png : composition comparison figures",
  "- sessionInfo.txt : package versions",
  "",
  "## Cycle 10 vs cycle 0",
  "The model fits a per-taxon amplification slope over PCR cycle plus a per-sample",
  "intercept. Predicting at cycle 0 (the original main analysis) is the strongest",
  "possible extrapolation, far below the observed 20-30 cycle range. Predicting at",
  "cycle 10 removes most later-cycle amplification bias while staying closer to the",
  "observed data.",
  "",
  "## Why cycle 10 is conservative",
  "It is NOT a claim that cycle 10 is the true starting community. It is a",
  "robustness check that mitigates amplification bias while avoiding the most",
  "aggressive extrapolation back to cycle 0, addressing the reviewer's concern.",
  "",
  paste0("Generated: ", Sys.time())
)
writeLines(readme_lines,
           pcr_versioned_path(file.path(out_dir, "README_cycle10_sensitivity.md")))


# -----------------------------------------------------------------------------
# CONSOLE SUMMARY
# -----------------------------------------------------------------------------
message("\n================ WORKFLOW COMPLETE ================")
message("Prediction PCR cycle: ", prediction_cycle, " (observed sequencing cycle: ", observed_cycle, ")")
message("18S fractions processed: ", paste(names(fractions_18s), collapse = ", "))
message("Taxa per fraction: ",
        paste(filtering_18s_df$size, filtering_18s_df$n_taxa, sep = "=", collapse = ", "))
message("All predicted compositions sum to 1: ", all(simplex_check_18s$sums_to_one))
message("Value sanity (neg/NA/Inf/>1): ",
        value_check_18s$n_negative, "/", value_check_18s$n_missing, "/",
        value_check_18s$n_infinite, "/", value_check_18s$n_gt_one)
message("Outputs written to: ", out_dir)
message("===================================================")
