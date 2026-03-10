# MPSA Pipeline README
## The Racial Economy of the United States: Economic Pessimism, Racecraft, and the Production of Racial Resentment

---

## Overview

This pipeline estimates the causal effect of economic pessimism on racial resentment, redistribution attitudes, and racial threat perceptions among white Americans using ANES 2024 data. The identification strategy escalates from a parametric GPS baseline through nonparametric causal forests, with mediation, CATE heterogeneity, covariate balance, and ML balance diagnostics as stress tests.

**Run scripts in order: 1 → 2 → 3 → 4 → 5.**

Scripts 2 through 5 are fully independent of each other. Once Script 1 has run and saved `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`, any downstream script can be re-run in isolation without re-running MICE.

---

## Dependencies

Install all required packages before running:

```r
install.packages(c(
  "tidyverse", "readr", "mice", "psych",
  "grf", "MatchIt", "cobalt", "mediation",
  "ggplot2", "patchwork"
))

# MLbalance (not on CRAN — install from GitHub)
install.packages("pak")
pak::pak("CetiAlphaFive/MLbalance")
```

---

## Scripts

### `mpsa_1_data_prep.R`
**Run first. All other scripts depend on this.**

Loads the raw ANES 2024 CSV, applies missing value codes, filters to the whites-only analytic sample, and constructs all treatment variables (binary and continuous), outcome variables, and control sets. Runs MICE (m=3, maxit=3) for multiple imputation. Saves the full environment to `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`.

**Inputs:**
- `anes_timeseries_2024_csv_20250808/anes_timeseries_2024_csv_20250808.csv`

**Outputs:**
- `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData` — full environment loaded by all downstream scripts. Contains: `anes_w`, `imp_list`, `df1`, `treats_bin`, `treats_cont`, `outcomes_all`, `controls_ols_raw`, `controls_cf_raw`, `treatment_labels`, `outcome_labels`, and all helper functions.

---

### `mpsa_2_gps_ols.R`
**Appendix baseline estimator. Run after Script 1.**

Generalized propensity score (GPS) → OLS dose-response estimation for all four continuous treatments across all seven outcomes, pooled across MICE imputations. Hirano & Imbens (2004) parametric approach. Serves as the interpretable political science baseline against which the nonparametric CF results are benchmarked. Results go in the appendix, not the paper body.

**Inputs:**
- `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`

**Outputs:**
- `mpsa_gps_results.csv` — full pooled GPS→OLS results table (appendix table)
- `mpsa_gps_coefplot.png` — coefficient plot, all treatments x outcomes (appendix figure)

---

### `mpsa_3_cf_main.R`
**Main results. Run after Script 1.**

Causal forest ATT (binary treatments, expanded CF covariate set) and causal forest continuous ATE (continuous treatments), both pooled across MICE imputations. This is the primary estimator for the paper. The CF ATT results are the main results table and figure in the paper body; the continuous ATE results go in the appendix alongside GPS→OLS for comparison.

**Inputs:**
- `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`

**Outputs:**
- `mpsa_cf_att_results.csv` — pooled CF ATT results, all treatments x outcomes **(main results table)**
- `mpsa_cf_att_coefplot.png` — coefficient plot, binary treatments **(main results figure)**
- `mpsa_cf_cont_results.csv` — pooled CF continuous ATE results (appendix table)
- `mpsa_cf_cont_coefplot.png` — coefficient plot, continuous treatments (appendix figure)

---

### `mpsa_4_robustness.R`
**Stress tests. Run after Script 1.**

Three robustness tests for the main CF ATT results:

**A) Mediation analysis** — tests whether racial resentment partially mediates the economic pessimism → redistribution opposition pathway. Four paths: national retrospective → RR → welfare, national retrospective → RR → aid to poor, personal retrospective → RR → welfare (the sign anomaly check), personal retrospective → RR → aid to poor. Uses `mediation` package (Imai et al.) with probit outcome models for binary outcomes.

**B) CATE heterogeneity plots** — refits causal forests for the two weakest average-effect outcomes (white status threat, immigration threat) and produces: CATE distribution histograms, best linear projection plots (which covariates predict treatment effect heterogeneity), and CATE-by-moderator smooths (ideology, income, education, party registration). Tests whether average nulls mask subgroup effects.

**C) Covariate balance** — propensity score matching (MatchIt, nearest neighbor, caliper 0.2) for all four binary treatments using the OLS control set. Produces love plots (standardized mean differences pre/post matching) and a full SMD table.

**Inputs:**
- `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`

**Outputs:**
- `mpsa_mediation_summary.csv` — ACME, ADE, total effect, % mediated for all four paths
- `mpsa_cate_status_threat_nat_retro_dist.png` — CATE distribution, white status threat
- `mpsa_cate_status_threat_nat_retro_blp.png` — BLP plot, white status threat
- `mpsa_cate_status_threat_nat_retro_mods.png` — moderator smooths, white status threat
- `mpsa_cate_imm_threat_nat_retro_dist.png` — CATE distribution, immigration threat
- `mpsa_cate_imm_threat_nat_retro_blp.png` — BLP plot, immigration threat
- `mpsa_cate_imm_threat_nat_retro_mods.png` — moderator smooths, immigration threat
- `mpsa_balance_love_[treatment].png` — love plot per treatment (4 files)
- `mpsa_balance_table.csv` — full SMD table, all treatments x covariates
- `mpsa_balance_summary.csv` — % covariates achieving |SMD| < 0.1 post-match per treatment

---

### `mpsa_5_mlbalance.R`
**Identification diagnostics. Run after Script 1.**

MLbalance (Rametta & Fuller 2026) classification permutation test and doubly-robust estimator divergence checks for all four binary treatments, run against both the OLS and CF covariate sets. Directly tests whether a random forest can distinguish treated from control on covariates better than chance, and whether the doubly-robust AIPW estimate converges with the CF ATT main results after adjustment. The key diagnostic: a large, significant DiM → DR collapse in the same direction as the CF ATT estimates confirms covariate adjustment is working correctly.

**Inputs:**
- `~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData`

**Outputs:**
- `mpsa_mlbalance_[treatment]_ols.png` — CPT + estimator plot, OLS controls (4 files)
- `mpsa_mlbalance_[treatment]_cf.png` — CPT + estimator plot, CF controls (4 files)
- `mpsa_mlbalance_summary.csv` — CPT statistic, p-value, PS SD ratio, all four estimator estimates and SEs, DiM vs DR divergence test per treatment x control set

---

## Variable Labels

### Treatments
| Script variable | Paper label |
|---|---|
| `T7_nat_econ_last_worse` | National retrospective economic pessimism |
| `T6_nat_econ_next_worse` | National prospective economic pessimism |
| `T8_pers_fin_last_worse` | Personal financial retrospective pessimism |
| `T3_fin_nextyear_worried` | Personal financial prospective worry |
| `C_nat_econ_last` | National retrospective economic pessimism (continuous) |
| `C_nat_econ_next` | National prospective economic pessimism (continuous) |
| `C_pers_fin_last` | Personal financial retrospective pessimism (continuous) |
| `C_fin_nextyear` | Personal financial prospective worry (continuous) |

### Outcomes
| Script variable | Paper label |
|---|---|
| `O1_rr_index` | Racial resentment index |
| `O4_welfare_binary` | Welfare spending opposition |
| `O5_aidpoor_binary` | Aid to poor opposition |
| `O6_status_jobs_threat` | White status / jobs threat |
| `O7_gov_favors_blacks_amt` | Government favors Blacks (amount) |
| `O8_imm_jobs_threat` | Immigration / jobs threat |
| `O9_gov_favors_blacks_dir` | Government favors Blacks (direction) |

---

## Paper → Script Mapping

| Paper section | Script | Key outputs |
|---|---|---|
| Data & Methods | `mpsa_1_data_prep.R` | Environment, sample description |
| Appendix baseline | `mpsa_2_gps_ols.R` | `mpsa_gps_coefplot.png` |
| Main results | `mpsa_3_cf_main.R` | `mpsa_cf_att_coefplot.png`, `mpsa_cf_att_results.csv` |
| Appendix CF continuous | `mpsa_3_cf_main.R` | `mpsa_cf_cont_coefplot.png` |
| Robustness: mediation | `mpsa_4_robustness.R` | `mpsa_mediation_summary.csv` |
| Robustness: heterogeneity | `mpsa_4_robustness.R` | CATE plot files |
| Robustness: balance | `mpsa_4_robustness.R` | Love plots, SMD tables |
| Identification diagnostics | `mpsa_5_mlbalance.R` | `mpsa_mlbalance_summary.csv`, CPT plots |
