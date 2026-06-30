# =============================================================================
# publication_pcr_bias_fido_functions.R
# -----------------------------------------------------------------------------
# Helper functions for the publication-ready PCR bias mitigation workflow that
# uses the Bayesian multinomial logistic-normal model implemented in `fido`
# (Silverman et al. 2022, the `pibble` model).
#
# WHY THIS FILE EXISTS
#   The original exploratory scripts in this folder
#   (predict_original_proportions_18s_phy*.R, fido_amp_effs.R, etc.) repeat the
#   same data-loading -> model-matrix -> pibble-fit -> predict pipeline three
#   times (once per size fraction) with copy/pasted code. To make the analysis
#   publishable and auditable we factor that shared logic into small, documented
#   functions here, and drive them from a single configurable workflow script
#   (publication_pcr_bias_fido_workflow_cycle10.R).
#
#   The ONLY scientific change relative to the originals is that the prediction
#   target PCR cycle is a parameter. The originals always extrapolated corrected
#   relative abundances back to PCR cycle 0 (the strongest possible
#   extrapolation, well below the observed 20-30 cycle range). Here we expose
#   `prediction_cycle` so we can run a more conservative sensitivity analysis at
#   PCR cycle 10 (see README in the output folder for the rationale).
#
# IMPORTANT
#   These functions do not read or write any of the original input files in
#   place. They only read existing inputs and write to a dedicated, newly
#   created output directory. They never overwrite an existing output file
#   (see pcr_versioned_path()).
#
# CONVENTIONS
#   * All paths are resolved relative to the project root via here::here() so
#     the workflow is reproducible regardless of working directory.
#   * The model design mirrors the originals exactly:
#         X <- t(model.matrix(~ cycle_num + sample_num - 1, data = meta))
#     i.e. a shared linear slope for `cycle_num` (the per-taxon log-ratio
#     amplification efficiency) plus a separate intercept for every biological
#     sample (and one shared intercept for the calibration pools). There is no
#     global intercept because of the `-1`.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr/tidyr/ggplot2/stringr/readr/tibble/purrr
  library(fido)        # Bayesian multinomial logistic-normal (pibble) model
  library(here)        # project-root-relative paths for reproducibility
})


# -----------------------------------------------------------------------------
# 0. Small assertion / safety helpers
# -----------------------------------------------------------------------------

#' Stop with an informative error if a required input file is missing.
#'
#' We fail loudly and early (rather than letting read.csv emit a cryptic error)
#' so that a reviewer re-running the workflow immediately sees which expected
#' input is absent.
pcr_require_file <- function(path, what = "input file") {
  if (!file.exists(path)) {
    stop(sprintf("Required %s not found:\n  %s\nCheck that you are running from the project root and that the intermediate inputs exist.",
                 what, path), call. = FALSE)
  }
  invisible(path)
}

#' Stop if required columns are absent from a data frame.
pcr_require_columns <- function(df, required_cols, df_name = "data frame") {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required column(s): %s",
                 df_name, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(df)
}

#' Return a non-clobbering path. If `path` already exists, append _v2, _v3, ...
#'
#' This enforces the project rule that we never overwrite an existing file.
#' Works for both files (keeps the extension) and directories (no extension).
pcr_versioned_path <- function(path) {
  if (!file.exists(path)) return(path)
  dir_part  <- dirname(path)
  base_part <- basename(path)
  has_ext   <- grepl("\\.[A-Za-z0-9]+$", base_part)
  ext       <- if (has_ext) sub("^.*(\\.[A-Za-z0-9]+)$", "\\1", base_part) else ""
  stem      <- if (has_ext) sub("(\\.[A-Za-z0-9]+)$", "", base_part) else base_part
  i <- 2L
  repeat {
    candidate <- file.path(dir_part, sprintf("%s_v%d%s", stem, i, ext))
    if (!file.exists(candidate)) return(candidate)
    i <- i + 1L
  }
}

#' Create a directory if needed, never overwriting an existing one.
pcr_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}


# -----------------------------------------------------------------------------
# 1. Data loading and alignment
# -----------------------------------------------------------------------------

#' Load a pre-filtered fido count table (taxa x samples).
#'
#' The intermediate inputs produced upstream (phyloseq filtering + ecdf
#' prevalence filtering) store taxa in a column named by their rank (e.g.
#' "Family" for 18S, "Genus" for COI) and samples as the remaining columns.
#'
#' @param path           Path to the count CSV (already filtered upstream).
#' @param taxon_col      Name of the taxonomy-rank column ("Family" or "Genus").
#' @param collapse_rare  Optional named recoding applied BEFORE summing, e.g.
#'                        c(Salpidae = "other"). This reproduces the 11/2024
#'                        decision in the all_and_subpools script to fold the
#'                        Salpidae family into the "other" catch-all, because
#'                        salp reads are gelatinous-tissue artifacts that behave
#'                        erratically across PCR cycles and destabilise the fit.
#' @return A numeric matrix with taxa as rownames and samples as colnames.
pcr_load_counts <- function(path, taxon_col = "Family", collapse_rare = NULL) {
  pcr_require_file(path, "fido count table")
  raw <- read.csv(path, header = TRUE, check.names = FALSE, row.names = 1)
  pcr_require_columns(raw, taxon_col, basename(path))

  df <- raw
  # Optional taxon recoding (e.g. Salpidae -> other), applied before summing so
  # that recoded taxa are merged into their target category.
  if (!is.null(collapse_rare)) {
    df[[taxon_col]] <- ifelse(df[[taxon_col]] %in% names(collapse_rare),
                              collapse_rare[df[[taxon_col]]],
                              df[[taxon_col]])
  }

  df <- df %>%
    dplyr::group_by(.data[[taxon_col]]) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.x, na.rm = TRUE)),
                     .groups = "drop") %>%
    tibble::column_to_rownames(taxon_col)

  # The CSV reader prefixes purely-numeric sample names with "X"; strip it so
  # sample names match the metadata Sample_name field exactly.
  colnames(df) <- gsub("^X", "", colnames(df))

  as.matrix(df)
}

#' Load sample metadata and align it to the columns of a count matrix.
#'
#' Returns metadata rows in EXACTLY the same order as the count-matrix columns,
#' subsetting both objects to their common samples. This is safer than the
#' original `meta[match(colnames(Y), meta$Sample_name), ]` idiom, which silently
#' produced NA rows (and therefore corrupted the model matrix) whenever a sample
#' name failed to match. Any dropped samples are reported so mismatches are
#' visible rather than silent.
#'
#' @return list(meta = aligned metadata, counts = aligned count matrix,
#'              dropped_samples = character vector of unmatched count columns)
pcr_load_and_align <- function(counts, meta_path) {
  pcr_require_file(meta_path, "sample metadata")
  meta <- read.csv(meta_path, header = TRUE)
  if ("X" %in% colnames(meta)) meta <- dplyr::select(meta, -dplyr::any_of("X"))
  pcr_require_columns(meta, c("Sample_name", "cycle_num", "sample_num"),
                      basename(meta_path))

  count_samples <- colnames(counts)
  matched       <- count_samples[count_samples %in% meta$Sample_name]
  dropped       <- setdiff(count_samples, meta$Sample_name)

  if (length(matched) == 0) {
    stop(sprintf("No count-table samples matched Sample_name in %s. Are the right metadata and count files paired?",
                 basename(meta_path)), call. = FALSE)
  }

  counts_aligned <- counts[, matched, drop = FALSE]
  meta_aligned   <- meta[match(matched, meta$Sample_name), , drop = FALSE]
  rownames(meta_aligned) <- NULL

  # Calibration samples (the dilution/cycle series) are mandatory: the cycle_num
  # slope is identified by them. Stop if they are absent.
  if (!any(meta_aligned$sample_num == "Calibration")) {
    stop("No calibration samples (sample_num == 'Calibration') remain after alignment; the amplification-bias slope cannot be estimated.",
         call. = FALSE)
  }

  list(meta = meta_aligned, counts = counts_aligned, dropped_samples = dropped)
}

#' Build the fido design matrix X (covariates x samples).
#'
#' Mirrors the originals exactly: a shared numeric `cycle_num` slope plus a
#' per-sample intercept (and one shared "Calibration" intercept), with no global
#' intercept (`-1`). We transpose because fido expects covariates in rows.
pcr_build_design <- function(meta) {
  X <- t(stats::model.matrix(~ cycle_num + sample_num - 1, data = meta))
  # model.matrix rows follow `meta` row order, which we have already aligned to
  # the count-matrix columns. After transposing, X loses those sample names
  # (columns become 1..n), so we restore them. This keeps X and Y column order
  # explicitly verifiable and does not affect prediction (which uses rownames).
  colnames(X) <- meta$Sample_name
  X
}


# -----------------------------------------------------------------------------
# 2. Model fitting
# -----------------------------------------------------------------------------

#' Optional Gamma (prior scale) selection via the log marginal likelihood.
#'
#' Reproduces the grid search used in the original scripts to justify the chosen
#' Gamma. This is OFF by default in the workflow (it refits the model many
#' times) but is provided so the choice is reproducible and defensible.
#'
#' @return a tibble with columns gamma, logML and the gamma maximising logML.
pcr_select_gamma <- function(Y, X, gamma_grid = c(0.1, 0.25, 0.5, 0.75, 1, 2, 3,
                                                  5, 8, 10, 15, 20, 50, 100, 200,
                                                  300, 400, 500, 700, 1000),
                              n_samples = 5000) {
  logML <- vapply(gamma_grid, function(g) {
    fit <- fido::pibble(Y, X, Gamma = g * diag(nrow(X)), n_samples = n_samples)
    fit$logMarginalLikelihood
  }, numeric(1))
  tibble::tibble(gamma = gamma_grid, logML = logML) %>%
    dplyr::mutate(is_best = .data$logML == max(.data$logML, na.rm = TRUE))
}

#' Fit the pibble (multinomial logistic-normal) model and return CLR + proportion
#' representations.
#'
#' We explicitly construct the default priors (upsilon, Omega, Xi, Theta) so the
#' prior specification is transparent in the manuscript rather than hidden inside
#' fido defaults. `Gamma` controls the prior scale on the regression
#' coefficients and is passed in (see pcr_select_gamma()).
#'
#' @return list(fit_clr = CLR-coordinate fit, fit_prop = proportion-coordinate
#'              fit, gamma = gamma used)
pcr_fit_pibble <- function(Y, X, gamma, n_samples = 10000, seed = 1989) {
  set.seed(seed)  # fido sampling is stochastic; fix the seed for reproducibility

  D <- nrow(Y)  # number of taxa (categories)
  # Default conjugate priors, written out explicitly for transparency:
  upsilon <- D + 3                       # IW degrees of freedom (weakly informative)
  Omega   <- diag(D)
  G       <- cbind(diag(D - 1), -1)      # ALR -> CLR contrast
  Xi      <- (upsilon - D) * G %*% Omega %*% t(G)
  Theta   <- matrix(0, D - 1, nrow(X))   # prior mean of coefficients = 0

  fit <- fido::pibble(Y, X,
                      Gamma     = gamma * diag(nrow(X)),
                      upsilon   = upsilon,
                      Theta     = Theta,
                      Xi        = Xi,
                      n_samples = n_samples)

  fit_clr  <- fido::to_clr(fit)          # centered log-ratio coordinates
  fit_prop <- fido::to_proportions(fit)  # simplex (proportion) coordinates

  list(fit_clr = fit_clr, fit_prop = fit_prop, gamma = gamma)
}


# -----------------------------------------------------------------------------
# 3. Amplification-efficiency summary
# -----------------------------------------------------------------------------

#' Extract the per-taxon amplification-efficiency coefficients.
#'
#' In CLR coordinates the `cycle_num` Lambda coefficient is the per-taxon,
#' per-cycle change in log-ratio abundance, i.e. the estimated PCR amplification
#' bias. Positive values = preferentially amplified taxa; negative = suppressed.
#' This mirrors fido_amp_effs.R but returns a tidy table for export.
pcr_amplification_efficiency <- function(fit_clr, size_label, marker) {
  # fido's summary() returns a one-element list named "Lambda"; as.data.frame()
  # flattens it to columns prefixed "Lambda." (e.g. Lambda.covariate). We strip
  # that prefix for a clean, manuscript-ready table.
  df <- as.data.frame(summary(fit_clr, pars = "Lambda"))
  colnames(df) <- sub("^Lambda\\.", "", colnames(df))
  pcr_require_columns(df, "covariate", "Lambda summary")
  df %>%
    dplyr::filter(.data$covariate == "cycle_num") %>%
    dplyr::mutate(size = size_label, marker = marker) %>%
    tibble::as_tibble()
}


# -----------------------------------------------------------------------------
# 4. Prediction at a configurable PCR cycle
# -----------------------------------------------------------------------------

#' Predict corrected relative abundances for every field sample at a chosen
#' PCR cycle.
#'
#' HOW THE PREDICTION CYCLE ENTERS THE MODEL
#'   The design has one shared `cycle_num` slope and one intercept per sample.
#'   For a given field sample s the linear predictor is
#'       eta(s, cycle) = intercept_s + slope * cycle.
#'   The original scripts predicted at cycle 0 by leaving the `cycle_num` row of
#'   the prediction vector at 0 (intercept only). Here we set that row to
#'   `prediction_cycle`, so corrected abundances are evaluated at, e.g., PCR
#'   cycle 10. Cycle 10 still removes most later-cycle amplification distortion
#'   but stays closer to the observed 20-30 cycle range than cycle 0, making it
#'   a deliberately conservative sensitivity analysis.
#'
#' For each field sample we return the model-predicted composition at
#' `prediction_cycle` together with the observed (raw) cycle-30 composition, so
#' downstream code can plot/compare them directly.
#'
#' @param fit_prop          proportion-coordinate pibble fit.
#' @param X                 design matrix used for the fit (for rownames).
#' @param Y                 count matrix (for the observed cycle-30 comparison).
#' @param prediction_cycle  PCR cycle to extrapolate corrected abundances to.
#' @param size_label        label for the size fraction (e.g. "0.2-0.5mm").
#' @param observed_cycle    PCR cycle at which field samples were sequenced (30).
#' @return long tibble: coord (taxon), replicate, n_reads (proportion),
#'         cycle_num, size, plus posterior quantile columns for predictions.
pcr_predict_at_cycle <- function(fit_prop, X, Y, prediction_cycle,
                                 size_label, observed_cycle = 30) {

  cycle_row <- which(rownames(X) == "cycle_num")
  if (length(cycle_row) != 1) {
    stop("Expected exactly one 'cycle_num' row in the design matrix.", call. = FALSE)
  }

  # Field samples to predict on: every per-sample intercept EXCEPT the shared
  # calibration intercept and the cycle_num slope itself.
  sample_rows <- rownames(X)
  field_samples <- sample_rows[!sample_rows %in% c("sample_numCalibration", "cycle_num")]
  field_samples <- field_samples[grepl("^sample_num", field_samples)]
  if (length(field_samples) == 0) {
    stop("No field samples (sample_num intercepts) found to predict on.", call. = FALSE)
  }

  Y_prop <- as.data.frame(Y) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / sum(.x)))  # per-sample raw proportions

  out <- vector("list", length(field_samples))

  for (i in seq_along(field_samples)) {
    s <- field_samples[i]
    base_name <- stringr::str_replace(s, "^sample_num", "")  # e.g. C1.T7.H9_S1

    # Build the prediction covariate vector: this sample's intercept = 1 and the
    # shared slope evaluated at `prediction_cycle`; everything else 0.
    X.tmp <- matrix(0, nrow(X), 1, dimnames = list(rownames(X), NULL))
    X.tmp[s, 1]         <- 1
    X.tmp[cycle_row, 1] <- prediction_cycle

    predicted <- predict(fit_prop, newdata = X.tmp, summary = TRUE) %>%
      dplyr::mutate(
        coord       = stringr::str_replace(.data$coord, "^prop_", ""),
        cycle_num   = prediction_cycle,
        size        = size_label,
        replicate   = paste("predicted", base_name)
      ) %>%
      dplyr::rename(n_reads = "mean")

    # Observed raw composition at the sequencing cycle (typically 30) for the
    # same sample's replicates, for side-by-side comparison.
    observed <- Y_prop %>%
      tibble::rownames_to_column("coord") %>%
      dplyr::select("coord", dplyr::starts_with(base_name)) %>%
      tidyr::pivot_longer(cols = -"coord", names_to = "replicate", values_to = "n_reads") %>%
      dplyr::mutate(size = size_label, cycle_num = observed_cycle)

    out[[i]] <- dplyr::bind_rows(predicted, observed)
  }

  dplyr::bind_rows(out) %>% tibble::as_tibble()
}


# -----------------------------------------------------------------------------
# 5. Diagnostic checks (Task 7)
# -----------------------------------------------------------------------------

#' Which PCR cycles are available for each pool, by marker.
pcr_diag_cycle_coverage <- function(meta, size_label, marker) {
  meta %>%
    dplyr::group_by(Pool = if ("Pool" %in% names(meta)) .data$Pool else NA,
                    sample_num = .data$sample_num) %>%
    dplyr::summarise(cycles_available = paste(sort(unique(.data$cycle_num)), collapse = ", "),
                     n_samples = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(size = size_label, marker = marker)
}

#' Does each taxon appear (non-zero reads) across all required calibration
#' cycles? Taxa absent from some calibration cycles have a poorly identified
#' amplification slope and are flagged here.
pcr_diag_taxon_calibration <- function(Y, meta, required_cycles, size_label, marker) {
  calib_samples <- meta$Sample_name[meta$sample_num == "Calibration"]
  res <- purrr::map_dfr(required_cycles, function(cyc) {
    cyc_samples <- meta$Sample_name[meta$sample_num == "Calibration" & meta$cycle_num == cyc]
    cyc_samples <- intersect(cyc_samples, colnames(Y))
    if (length(cyc_samples) == 0) {
      return(tibble::tibble(coord = rownames(Y), cycle_num = cyc, present = FALSE))
    }
    present <- rowSums(Y[, cyc_samples, drop = FALSE]) > 0
    tibble::tibble(coord = rownames(Y), cycle_num = cyc, present = present)
  })
  res %>%
    tidyr::pivot_wider(names_from = "cycle_num", values_from = "present",
                       names_prefix = "cycle_") %>%
    dplyr::mutate(
      present_in_all_calibration_cycles = dplyr::if_all(dplyr::starts_with("cycle_"), ~ .x),
      size = size_label, marker = marker
    )
}

#' Confirm the model matrices have the expected, mutually consistent dimensions.
pcr_diag_model_matrix <- function(Y, X, size_label, marker) {
  ok_align <- identical(colnames(Y), colnames(X))
  tibble::tibble(
    size = size_label, marker = marker,
    n_taxa = nrow(Y),
    n_samples_Y = ncol(Y),
    n_samples_X = ncol(X),
    n_covariates = nrow(X),
    samples_aligned = ok_align
  )
}

#' Check that posterior-predicted proportions sum to ~1 within each sample.
pcr_diag_simplex <- function(prediction_df, tol = 1e-6) {
  prediction_df %>%
    dplyr::filter(grepl("^predicted ", .data$replicate)) %>%
    dplyr::group_by(.data$size, .data$replicate, .data$cycle_num) %>%
    dplyr::summarise(prop_sum = sum(.data$n_reads, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(sums_to_one = abs(.data$prop_sum - 1) < tol)
}

#' Flag any unexpected negative, missing, infinite, or extreme corrected values.
pcr_diag_values <- function(prediction_df) {
  pred <- prediction_df %>% dplyr::filter(grepl("^predicted ", .data$replicate))
  tibble::tibble(
    n_predicted_values = nrow(pred),
    n_negative = sum(pred$n_reads < 0, na.rm = TRUE),
    n_missing  = sum(is.na(pred$n_reads)),
    n_infinite = sum(is.infinite(pred$n_reads)),
    n_gt_one   = sum(pred$n_reads > 1 + 1e-6, na.rm = TRUE)
  )
}
