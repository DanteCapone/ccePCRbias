suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

folder <- here("PCR_bias_correction", "data", "dada2_output")

files <- list.files(
  folder,
  pattern = "^dada2_summary_(coi|18s)_run[12]\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(files) == 0) {
  stop("No DADA2 summary CSV files found in ", folder)
}

read_one <- function(fp) {
  nm <- basename(fp) |> tolower()
  m <- stringr::str_match(nm, "^dada2_summary_((coi|18s))_run([12])\\.csv$")
  if (any(is.na(m))) stop("Unexpected filename: ", fp)
  barcode <- m[, 2]
  run <- m[, 4]

  df <- readr::read_csv(fp, show_col_types = FALSE)
  names_l <- tolower(names(df))
  nonchim_idx <- which(stringr::str_detect(names_l, "^nonchim"))
  if (length(nonchim_idx) == 0) stop("No 'nonchim' column in ", fp)
  input_idx <- which(stringr::str_detect(names_l, "^input"))
  if (length(input_idx) == 0) stop("No 'input' column in ", fp)

  tibble(
    barcode = barcode,
    run = run,
    input = as.numeric(df[[input_idx[1]]]),
    nonchim = as.numeric(df[[nonchim_idx[1]]])
  )
}

all_nonchim <- purrr::map_dfr(files, read_one)

stats <- all_nonchim %>%
  group_by(barcode) %>%
  summarise(
    n = sum(!is.na(nonchim)),
    mean_nonchim = mean(nonchim, na.rm = TRUE),
    sd_nonchim = sd(nonchim, na.rm = TRUE),
    .groups = "drop"
  )

cat("DADA2 nonchim summary (combined run1 + run2):\n")
stats %>% arrange(barcode) %>%
  pwalk(function(barcode, n, mean_nonchim, sd_nonchim) {
    cat(sprintf("%s: mean = %.0f, sd = %.0f (n = %d)\n", toupper(barcode), mean_nonchim, sd_nonchim, n))
  })

# Also print the tibble for convenience
print(stats)

# Totals across runs per barcode
totals <- all_nonchim %>%
  group_by(barcode) %>%
  summarise(
    total_input = sum(input, na.rm = TRUE),
    total_nonchim = sum(nonchim, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nTotals (combined run1 + run2):\n")
totals %>% arrange(barcode) %>%
  pwalk(function(barcode, total_input, total_nonchim) {
    cat(sprintf(
      "A total of %s sequences were recovered at the %s marker, %s of which passed quality and chimera filtering and were denoised into ASV.\n",
      format(total_input, big.mark = ",", scientific = FALSE),
      toupper(barcode),
      format(total_nonchim, big.mark = ",", scientific = FALSE)
    ))
  })

