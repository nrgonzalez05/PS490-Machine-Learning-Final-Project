# ============================================================
# mpsa_5_mlbalance.R
# Sam Fuller stuff 
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MLbalance)
})

select <- dplyr::select
load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")
cat("Outcomes:", length(outcomes_all), "| Primary: O1_rr_index\n")

# ============================================================
# HELPER: run balance() for one treatment x covariate set
# ============================================================

run_mlbalance <- function(df, treat_var, ctrls,
                           outcome_var = "O1_rr_index",
                           seed = 123, ctrl_label = "") {

  treat_name   <- treatment_labels[treat_var]
  outcome_name <- outcome_labels[outcome_var]

  d <- df %>%
    select(all_of(c(treat_var, outcome_var, ctrls))) %>%
    drop_na()

  cat("\n", strrep("-", 60), "\n")
  cat("Treatment:", treat_name, "\n")
  cat("Outcome:  ", outcome_name, "\n")
  cat("Controls: ", ctrl_label, "| n =", nrow(d), "\n")
  cat(strrep("-", 60), "\n")

  W <- as.numeric(d[[treat_var]])
  Y <- as.numeric(d[[outcome_var]])

  X_raw <- model.matrix(
    as.formula(paste("~", paste(ctrls, collapse = " + "))),
    data = d
  )[, -1, drop = FALSE]
  keep <- apply(X_raw, 2, function(z) sd(z, na.rm = TRUE) > 0)
  X    <- X_raw[, keep, drop = FALSE]

  cat("Covariate matrix: n =", nrow(X), "| p =", ncol(X), "\n\n")

  set.seed(seed)
  tryCatch(
    balance(Y = Y, W = W, X = X),
    error = function(e) {
      cat("ERROR in balance():", conditionMessage(e), "\n")
      NULL
    }
  )
}

# ============================================================
# all 5 analysis outcomes
# ============================================================

mlb_outcomes <- outcomes_all

# ============================================================
# 4 treatments x 5 outcomes x 2 control sets = 40 cells
# ============================================================

treatments <- list(
  list(var = "T7_nat_econ_last_worse",
       label = treatment_labels["T7_nat_econ_last_worse"]),
  list(var = "T8_pers_fin_last_worse",
       label = treatment_labels["T8_pers_fin_last_worse"]),
  list(var = "T3_fin_nextyear_worried",
       label = treatment_labels["T3_fin_nextyear_worried"]),
  list(var = "T6_nat_econ_next_worse",
       label = treatment_labels["T6_nat_econ_next_worse"])
)

all_results <- list()

for (tr in treatments) {
  for (out_var in mlb_outcomes) {

    out_name <- outcome_labels[out_var]

    cat("\n\n", strrep("=", 60), "\n")
    cat(tr$label, "\n")
    cat("Outcome:", out_name, "\n")
    if (out_var == "O1_rr_index") cat("[PRIMARY OUTCOME]\n")
    cat(strrep("=", 60), "\n")

    stem <- paste0(
      gsub("[^a-zA-Z0-9]", "_", gsub("\\s+", "_", trimws(tr$label))),
      "_",
      gsub("[^a-zA-Z0-9]", "_", gsub("\\s+", "_", trimws(out_name)))
    )

    # (a) OLS controls — primary diagnostic
    b_ols <- run_mlbalance(df1, tr$var, controls_ols_raw,
                            outcome_var = out_var,
                            ctrl_label  = "OLS controls (small, interpretable)")
    if (!is.null(b_ols)) {
      cat("\n--- Summary:", tr$label, "| OLS controls | Outcome:", out_name, "---\n")
      print(summary(b_ols))
      png(paste0("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_mlbalance_",
                 stem, "_ols.png"),
          width = 900, height = 700, res = 120)
      plot(b_ols, main = paste0(tr$label, "\n", out_name, " (OLS controls)"))
      dev.off()
      all_results[[paste0(tr$var, "_", out_var, "_ols")]] <- list(
        treatment     = tr$var,
        outcome       = out_var,
        label         = tr$label,
        outcome_label = out_name,
        controls      = "OLS (small)",
        obj           = b_ols
      )
    }

    # (b) CF controls — secondary diagnostic
    b_cf <- run_mlbalance(df1, tr$var, controls_cf_raw,
                           outcome_var = out_var,
                           ctrl_label  = "CF controls (expanded ML set)")
    if (!is.null(b_cf)) {
      cat("\n--- Summary:", tr$label, "| CF controls | Outcome:", out_name, "---\n")
      print(summary(b_cf))
      png(paste0("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_mlbalance_",
                 stem, "_cf.png"),
          width = 900, height = 700, res = 120)
      plot(b_cf, main = paste0(tr$label, "\n", out_name, " (CF controls)"))
      dev.off()
      all_results[[paste0(tr$var, "_", out_var, "_cf")]] <- list(
        treatment     = tr$var,
        outcome       = out_var,
        label         = tr$label,
        outcome_label = out_name,
        controls      = "CF (expanded)",
        obj           = b_cf
      )
    }
  }
}

# ============================================================
# sum
# ============================================================

extract_mlb_safe <- function(res_entry) {
  obj <- res_entry$obj
  if (is.null(obj)) return(NULL)
  tryCatch({
    est_dr  <- obj$aipw$estimate
    se_dr   <- obj$aipw$std.err
    est_dim <- obj$dim$estimate

    tibble(
      treatment           = treatment_labels[res_entry$treatment],
      outcome             = res_entry$outcome_label,
      primary             = res_entry$outcome == "O1_rr_index",
      controls            = res_entry$controls,
      n                   = obj$n,
      n_treated           = obj$n_treated,
      cpt_stat            = obj$balance_test$teststat,
      cpt_p               = obj$balance_test$pval[["ferns"]],
      balance_result      = ifelse(obj$balance_test$pval[["ferns"]] < 0.05,
                                   "FAIL", "PASS"),
      est_dim             = est_dim,
      se_dim              = obj$dim$std.err,
      est_ipw             = obj$ipw$estimate,
      se_ipw              = obj$ipw$std.err,
      est_outcome         = obj$aipw_const$estimate,
      se_outcome          = obj$aipw_const$std.err,
      est_dr              = est_dr,
      se_dr               = se_dr,
      overlap_flag        = obj$overlap_flag,
      n_extreme           = obj$n_extreme,
      dr_z                = est_dr / se_dr,
      dr_p                = 2 * pnorm(-abs(est_dr / se_dr)),
      dr_sig              = (2 * pnorm(-abs(est_dr / se_dr))) < 0.05,
      dim_attenuation_pct = round((est_dim - est_dr) / est_dim * 100, 1)
    )
  }, error = function(e) {
    cat("Note: auto-extraction failed for", res_entry$treatment,
        "/", res_entry$controls, "—", conditionMessage(e), "\n")
    tibble(treatment = treatment_labels[res_entry$treatment],
           outcome   = res_entry$outcome_label,
           controls  = res_entry$controls,
           note      = conditionMessage(e))
  })
}

mlb_summary <- map_dfr(all_results, extract_mlb_safe)

cat("\n\n", strrep("=", 60), "\n")
cat("MLBALANCE SUMMARY TABLE\n")
cat(strrep("=", 60), "\n")
cat("\nPRIMARY OUTCOME (O1):\n")
print(mlb_summary %>% filter(primary == TRUE) %>%
        select(-primary), n = Inf)
cat("\nSECONDARY OUTCOMES:\n")
print(mlb_summary %>% filter(primary == FALSE | is.na(primary)) %>%
        select(-primary), n = Inf)

write_csv(mlb_summary,
          "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_mlbalance_summary.csv")
cat("\nSaved: mpsa_mlbalance_summary.csv\n")

cat("\n", strrep("=", 60), "\n")
cat("INTERPRETATION NOTES\n")
cat(strrep("=", 60), "\n")
cat("
CPT FAIL (p < .05): Expected for observational data.
  Key result is DiM -> DR collapse and whether DR remains
  significant after full adjustment.

OLS vs CF DR divergence: OLS = primary diagnostic.
  CF attenuation expected for national treatments (bad
  controls problem in linear DR). Does not invalidate
  CF ATT forest estimates.

Prop-Adj vs DR divergence: warning sign of model
  disagreement. Flag in methods if significant.

Overlap warnings: OW estimates are conservative bound.
  High SD ratio treatments will show more extreme scores.

Outcome hierarchy: O1 (racial resentment) is primary.
  O4/O5/O6/O7 are secondary per pre-specified plan.
")

cat("Done. All outputs saved.\n")
