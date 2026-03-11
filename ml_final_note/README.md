# Testing Economic Attitudes with Causal Forests
### Nicholas R. Gonzalez — Northwestern University

Analysis pipeline for the research note examining economic pessimism as a driver of racial resentment, redistribution attitudes, and status threat beliefs among white Americans (2024 ANES). The project compares generalized propensity score OLS (GPS-OLS) against continuous causal forests (GRF) across a 2×2 treatment design crossing economic referent (national vs. personal) with temporal orientation (retrospective vs. prospective).

---

## Requirements

**R packages:** `tidyverse`, `mice`, `psych`, `grf`, `MatchIt`, `cobalt`, `MLbalance`, `rbounds`, `ggplot2`, `patchwork`

**Data:** `anes_timeseries_2024_csv_20250808.csv` — place in working directory before running script 1.

**Execution order:** Scripts must be run sequentially (1 → 10). All scripts after script 1 load from `mpsa_environment.RData`.

---

## Script Overview

### `mpsa_1_data_prep.R`
Loads and cleans the 2024 ANES, restricts to non-Hispanic white respondents (N = 3,946), constructs all treatment variables (binary and continuous 1–5 scale), outcome variables (racial resentment index, welfare opposition, aid-to-poor opposition, white status/jobs threat, govt. favors Blacks), and OLS/CF control sets. Runs MICE multiple imputation (m = 3, maxit = 3) and saves the full environment to `mpsa_environment.RData`.

**Output:** `mpsa_environment.RData`

---

### `mpsa_2_gps_ols.R`
Estimates the benchmark GPS-OLS models following Hirano and Imbens (2004). Fits a generalized propensity score on the continuous 1–5 treatment scale, then estimates a dose-response curve via OLS conditioning on the GPS. Runs all four treatments × five outcomes using the theory-derived OLS control set (8 variables). Pools results across MICE imputations using Rubin's rules.

**Output:** `mpsa_gps_results.csv`, `mpsa_gps_coefplot.png`

---

### `mpsa_3_cf_main.R`
Estimates the primary continuous causal forest models using the GRF package (Wager and Athey 2018; Athey, Tibshirani, and Wager 2019). Targets the local average marginal effect (dY/dW) on the continuous 1–5 treatment scale using sample splitting and recursive partitioning. Uses the expanded CF covariate set (40+ variables). Pools ATEs across MICE imputations using Rubin's rules. This is the primary estimator and main results of the paper.

**Output:** `mpsa_cf_cont_results.csv`, `mpsa_cf_cont_coefplot.png`

---

### `mpsa_4_robustness.R`
Two robustness checks. First, CATE heterogeneity analysis examining treatment effect variation across the distribution for the primary treatment (national retrospective pessimism) across all five outcomes, with partisan media as a potential amplification moderator. Second, covariate balance diagnostics using standardized mean differences (SMD) and love plots via `cobalt`, documenting how much raw imbalance exists before adjustment.

**Output:** `mpsa_cate_[outcome]_nat_retro_[dist/blp/mods].png`, `mpsa_balance_love_[treatment].png`, `mpsa_balance_table.csv`, `mpsa_balance_summary.csv`

---

### `mpsa_5_mlbalance.R`
Implements the ML-based balance test from Fuller and Rametta (working paper) via the `MLbalance` package. Tests whether treated and untreated units remain distinguishable after covariate adjustment — the conditional permutation test (CPT). Runs for all four treatments under both the OLS control set and the expanded CF control set. CPT failure (documented for all treatments) motivates the Rosenbaum sensitivity analysis in script 6.

**Output:** `mpsa_mlbalance_summary.csv`, `mpsa_mlbalance_attenuation.png`

---

### `mpsa_6_rosenbaum.R`
Sensitivity analysis via Rosenbaum bounds (Rosenbaum 2002). For each treatment × outcome combination, estimates the critical Γ — the magnitude of hidden bias (expressed as an odds multiplier on treatment assignment) required to explain away the result entirely as confounding. Uses `rbounds` with nearest-neighbor matching via `MatchIt`. The primary finding is Γ = 3.0 for national retrospective pessimism on racial resentment.

**Output:** `mpsa_sensitivity_critical_gamma.png`, `mpsa_sensitivity_summary.csv`

---

### `mpsa_8_diag.R`
Causal forest validity diagnostics for continuous treatments. Runs `test_calibration()` from the GRF package for each treatment, testing whether the forest's nuisance models are well-specified. Also computes the Ŵ (predicted treatment propensity) distribution for each treatment to assess overlap and identifying variation. National prospective pessimism fails both diagnostics, indicating no credible causal identification for that treatment.

**Output:** `mpsa_cf_cont_calibration.csv`, `mpsa_what_density_[treatment].png`

---

### `mpsa_9_visuals.r`
Consolidated visualization script. Reads all CSV outputs from prior scripts and produces the final publication-ready figures. Requires `mpsa_gps_results.csv`, `mpsa_cf_cont_results.csv`, `mpsa_mlbalance_summary.csv`, `mpsa_sensitivity_summary.csv`, and `mpsa_cf_cont_calibration.csv` to all be present in `outputs_ml/` before running.

**Output:** `fig_01` through `fig_09` PNGs to `final_outputs/`

| Figure | Content |
|--------|---------|
| fig_01 | GPS-OLS vs. CF estimates — primary outcome (racial resentment) |
| fig_02 | GPS-OLS vs. CF estimates — all outcomes |
| fig_03 | DIM → IPW → doubly-robust attenuation by control set |
| fig_04 | Dose-response curve — primary outcome |
| fig_05 | Dose-response curve — all outcomes |
| fig_06 | Rosenbaum sensitivity: critical Γ by treatment and outcome |
| fig_07 | Continuous CF calibration coefficients |
| fig_08 | Ŵ density — propensity support |
| fig_09 | Ŵ support overlap |

---

### `mpsa_10_importance.R`
Variable importance analysis using GRF's `variable_importance()` function, which returns a depth-weighted split frequency for each covariate across all trees. Normalized to [0, 1] within each forest. Runs on the first MICE imputation using the same forest specification as the main results (2,000 trees, min.node.size = 5, honesty = TRUE). Covers all four treatments × five outcomes. The primary finding is that hostile sexism dominates the racial resentment forest for national retrospective pessimism, with psychological distress variables in the second tier.

**Output:** `fig_10_varimp_primary.png`, `fig_11_varimp_outcome_grid.png`, `fig_12_varimp_welfare_flag.png`, `mpsa_varimp_full.csv`

