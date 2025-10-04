# ccePCRbias
Correcting PCR amplification bias in eDNA metabarcoding for quantitative zooplankton community inference.

## Summary
This repository contains the analysis code and minimal input data for a manuscript focused on quantifying and correcting PCR amplification bias in eDNA metabarcoding of zooplankton. We combine calibration/experimental analyses, simulations of amplification efficiencies, and compositional modeling (fido/pibble) to derive bias-corrected taxonomic proportions. We evaluate improvements against independent benchmarks (e.g., Zooscan biomass) and examine robustness across size fractions and barcodes (18S rRNA, COI).

## What’s included (no FASTA/FASTQ)
We include only the minimal CSV inputs needed to run the main analyses. Raw sequencing files (FASTA/FASTQ) and private data are excluded.

- inst/extdata/physical_environmental_data/
  - env_metadata_impute_phyloseq_6.9.2023.csv (and related env files as used)
- inst/extdata/raw_reads/
  - ASV_table_18s_run1.csv(.gz), ASV_table_18s_run2.csv(.gz)
  - ASV_table_coi_run1.csv(.gz), ASV_table_coi_run2.csv(.gz)
- inst/extdata/taxa_files/
  - blast_metazoo_18s.csv
  - blast_metazoo_coi.csv
- inst/extdata/fido/ and inst/extdata/fido/phy/
  - intermediate tables used by the fido pipelines
- inst/extdata/Zooscan/ (if referenced by methods comparisons)
  - zooscan_by_sample_biomass_esd.csv

Note: Large CSVs may be stored as `.csv.gz`. Analysis code reads gzipped CSVs transparently.

## How to reproduce
1) Clone and open the project root.
2) Start R in the project root (R >= 4.2 recommended).
3) Restore dependencies:
```r
install.packages("renv")
renv::init()          # first time
renv::restore()       # thereafter
install.packages("devtools")
devtools::load_all()  # expose functions in R/
```
4) Run the main analysis:
```r
source("scripts/PCR_Bias_Correction_Main_Analysis.R")
```
Other key scripts:
- `scripts/Experiment_analysis/calibration_experiment_analysis.R`
- `scripts/simulations/amp_eff_simulations.R`
- `scripts/Q5_methods_comparison/*.R`

Outputs are written to:
- `figures/`
- `data/processed/`

## Repository structure
- `R/` … exported helper functions (formerly in `scripts/helpful_functions/`)
- `scripts/`
  - `PCR_Bias_Correction_Main_Analysis.R` (entry point)
- `inst/extdata/` … small, versioned inputs used by scripts
- `data/`
  - `raw/` (gitignored)
  - `processed/` (generated outputs)
- `figures/` (generated outputs)
- `renv/`, `renv.lock` … locked dependency versions
- `DESCRIPTION`, `NAMESPACE`, `.Rbuildignore`, `.Rprofile`, `.gitignore`

## Methods (brief)
- Amplification efficiency estimation and bias-correction workflows.
- Compositional modeling with fido (pibble; hyperparameter scanning; posterior summaries).
- Integration with phyloseq; transformations via microbiome.
- Comparison of bias-corrected outputs with external benchmarks (e.g., Zooscan).

## Data policy
- Raw FASTA/FASTQ, BAM/SAM, and private data are excluded.
- Only small derivative inputs are versioned under `inst/extdata/`.
- Scripts write to `data/processed/` and must never overwrite `data/raw/` (see `.gitignore`).

## Citation
- Paper: “PCR bias correction for zooplankton eDNA metabarcoding” (in prep).
- Please cite this repository and the paper once published. A CITATION.cff will be added upon acceptance.

## License
MIT License (see LICENSE).

## Contact
Questions: open a GitHub issue or contact the corresponding author(s) listed in the manuscript.
