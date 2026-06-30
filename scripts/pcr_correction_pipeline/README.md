# PCR Bias Correction Pipeline

This folder is the clean, reproducible entry point for the fido-based PCR-bias
correction used in the manuscript. It does **not** duplicate the modelling
code; instead it provides shared configuration, documentation, and an
orchestrator that runs the existing pipeline scripts in the correct order.

## The `target_cycle` parameter

The correction fits a Bayesian multinomial logistic-normal model
(`fido::pibble`) in which community composition changes linearly (in CLR
space) with PCR `cycle_num`. The "original" pre-amplification composition is
recovered by evaluating that fitted line at a chosen cycle.

`config.R` sets a single parameter:

```r
target_cycle <- 10   # back-extrapolate to PCR cycle 10 (was 0 in the original submission)
```

Every predict script reads this value (defaulting to `10` if run standalone)
and applies it in two coordinated places:

1. **Prediction** – the `cycle_num` row of the prediction design vector
   `X.tmp` is set to `target_cycle` (`X.tmp["cycle_num", ] <- target_cycle`).
   This is what actually changes the predicted proportions.
2. **Label** – the output column `cycle_num` is tagged with `target_cycle` so
   downstream code can filter the predictions.

The main analysis filters predictions with `filter(cycle_num == target_cycle)`
(currently `== 10`).

### Important: amplification efficiencies are cycle-invariant

`fido_amp_effs.R` reports the **slope** of the cycle effect
(`Lambda.covariate == "cycle_num"`). A slope does not depend on where the
line is evaluated, so amplification-efficiency values are **unchanged** by
`target_cycle`. The script is re-run here only to regenerate a consistent,
date-stamped set of outputs.

## Which scripts feed the paper?

Only two predict variants produce files read by
`PCR_Bias_Correction_Main_Analysis.R`:

| Marker | Script | Output pattern (in `data/predicted_og/`) |
|--------|--------|-------------------------------------------|
| 18S | `predict_original_proportions_18s_phy_all_and_sub_pools_nocollodaria.R` | `predicted_og_18s_<date>_sX_phy_all_and_subpools_nocollodaria.csv` |
| COI | `predict_original_proportions_coi.R` | `predicted_og_coi_<date>_sX.csv` |

All other `predict_original_proportions_*` variants (`_phy`, `_order`,
`_clr`, `_asv_level`, `_nofilt`, COI `_phy*`, and everything in `sub_pools/`)
are also parameterized to `target_cycle` for consistency, but their outputs
are not consumed by the main analysis.

## How to reproduce (paper-feeding products)

From the `PCR_bias_correction/` project root:

```bash
Rscript scripts/pcr_correction_pipeline/run_paper_pipeline.R
```

This runs, in order:

1. 18S prediction (nocollodaria)
2. COI prediction
3. `fido_amp_effs.R`

Outputs are date-stamped, so prior runs are preserved. The main analysis
auto-loads the most recent date-stamped file.

## Notes

- Each predict step is a long Bayesian MCMC fit (`n_samples = 10000`); expect
  several minutes per script.
- The predict scripts call `beepr::beep()` for an audible completion cue. If
  `beepr` is not installed, install it or remove those lines.
- Raw inputs live under `data/fido/`; nothing in `data/raw_data/` is modified.
