# =============================================================================
# config.R  --  Shared configuration for the PCR-bias-correction pipeline
# =============================================================================
# Source this file at the top of every pipeline step. It defines the single
# parameter that controls how far back the fido model extrapolates the
# community composition, plus the project-root anchor used by here().
#
# WHY target_cycle MATTERS
# ------------------------
# Each predict_original_proportions_*.R script fits a Bayesian multinomial
# logistic-normal model (fido::pibble) in which the community composition
# changes linearly (in CLR space) with PCR `cycle_num`. To recover the
# pre-amplification composition we evaluate that fitted line at a chosen
# cycle by setting the `cycle_num` row of the prediction design vector
# (X.tmp) to `target_cycle`.
#
#   target_cycle <- 0   -> full extrapolation back to cycle 0 (original paper)
#   target_cycle <- 10  -> more conservative back-extrapolation to cycle 10
#
# NOTE: amplification efficiencies (fido_amp_effs.R) are the SLOPE of that
# line (Lambda for covariate "cycle_num"). A slope does not depend on where
# the line is evaluated, so amp-eff values are unchanged by target_cycle.
# =============================================================================

# --- Back-extrapolation target (PCR cycle) -----------------------------------
target_cycle <- 10

# --- Project-root anchor ------------------------------------------------------
# All pipeline scripts resolve paths relative to the PCR_bias_correction/
# project root via here(). The .here anchor file lives at that root.
if (requireNamespace("here", quietly = TRUE)) {
  library(here)
}

# Each predict_original_proportions_*.R script is self-contained: it defines
# `target_cycle <- 10` via a guard (`if (!exists("target_cycle"))`) so it can
# run standalone in RStudio, and applies it with:
#     X.tmp[["cycle_num", ]] <- target_cycle
# When sourced through run_paper_pipeline.R, this config sets `target_cycle`
# first, so the scripts inherit the value defined here.

message(sprintf("[pcr_correction_pipeline] target_cycle = %d", target_cycle))
