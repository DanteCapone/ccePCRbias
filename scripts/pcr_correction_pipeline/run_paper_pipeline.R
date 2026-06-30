# =============================================================================
# run_paper_pipeline.R
# -----------------------------------------------------------------------------
# Reproduces the PCR-bias-corrected data products that feed the paper figures,
# back-extrapolated to `target_cycle` (set in config.R).
#
# Steps:
#   1. predict original proportions, 18S (nocollodaria)  -> data/predicted_og/
#   2. predict original proportions, COI                 -> data/predicted_og/
#
# Amplification efficiencies (fido_amp_effs.R) are the SLOPE of the cycle
# effect and do NOT change with target_cycle, so they are intentionally NOT
# re-run here. Set run_amp_effs <- TRUE below only if you want fresh
# (cycle-invariant) amp-eff outputs; note it regenerates only the non-_nofilt
# files and requires the plots/pre_processing directory.
#
# Each predict step writes DATE-STAMPED outputs, so previous runs are never
# overwritten. The main analysis (PCR_Bias_Correction_Main_Analysis.R) then
# loads the most recent date-stamped file automatically.
#
# Usage (from the PCR_bias_correction/ project root):
#   Rscript scripts/pcr_correction_pipeline/run_paper_pipeline.R
# =============================================================================

# --- Shared configuration (defines target_cycle, here(), helpers) ------------
source(here::here("scripts/pcr_correction_pipeline/config.R"))

predict_dir <- here::here("scripts/pre-processing/pcr_bias_mitigation_prep_fido")

run_step <- function(path, label) {
  message(sprintf("\n=== [%s] sourcing %s ===", label, basename(path)))
  t0 <- Sys.time()
  source(path, local = new.env())
  message(sprintf("=== [%s] done in %.1f min ===",
                  label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# 1. 18S (paper-feeding variant) ----------------------------------------------
run_step(file.path(predict_dir,
                   "predict_original_proportions_18s_phy_all_and_sub_pools_nocollodaria.R"),
         "18S")

# 2. COI (paper-feeding variant) ----------------------------------------------
run_step(file.path(predict_dir, "predict_original_proportions_coi.R"), "COI")

# 3. Amplification efficiencies (slope coefficients; cycle-invariant; optional)
run_amp_effs <- FALSE
if (isTRUE(run_amp_effs)) {
  dir.create(here::here("plots/pre_processing"), showWarnings = FALSE, recursive = TRUE)
  run_step(file.path(predict_dir, "fido_amp_effs.R"), "AMP_EFFS")
}

message("\nAll paper-feeding pipeline steps complete.")
