
suppressPackageStartupMessages({
  library(tidyverse)
  library(grf)
  library(ggplot2)
  library(patchwork)
  library(broom)   # for tidy() on lm_robust objects from test_calibration()
})

select <- dplyr::select

load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")

# ============================================================
# DIAGNOSTIC FOREST FUNCTION
# Same spec as cf_cont_fn but returns the full forest object
# plus diagnostics rather than just the ATE estimate
# ============================================================

cf_cont_diag_fn <- function(df, y, t, ctrls, seed = 123, min_n = 500) {
  
  d <- df %>% select(all_of(c(y, t, ctrls))) %>% drop_na()
  
  if (nrow(d) < min_n) {
    warning(paste("Skipping", y, "/", t, "— n =", nrow(d)))
    return(NULL)
  }
  
  Y <- d[[y]]
  W <- as.numeric(d[[t]])
  X <- drop_zero_var(make_X(d, ctrls))
  
  set.seed(seed)
  forest <- causal_forest(
    X             = X,
    Y             = Y,
    W             = W,
    num.trees     = 2000,
    min.node.size = 5,
    honesty       = TRUE,
    seed          = seed
  )
  
  # ----------------------------------------------------------
  # test_calibration() fits:
  #   Y ~ mean.forest.pred + differential.forest.pred * (W - W.hat)
  # mean.forest.pred coef ~ 1 means the average CATE is well-calibrated
  # differential.forest.pred coef ~ 1 means heterogeneity is well-captured
  # Both significantly != 0 is good; mean coef far from 1 is a red flag
  # ----------------------------------------------------------
  cal <- tryCatch(
    test_calibration(forest),
    error = function(e) NULL
  )
  
  cal_df <- if (!is.null(cal)) {
    tidy_cal <- tryCatch(
      broom::tidy(cal),
      error = function(e) NULL
    )
    
    if (!is.null(tidy_cal)) {
      mean_row <- tidy_cal[tidy_cal$term == "mean.forest.prediction", ]
      diff_row <- tidy_cal[tidy_cal$term == "differential.forest.prediction", ]
      
      tibble(
        mean_pred_coef         = mean_row$estimate,
        mean_pred_se           = mean_row$std.error,
        mean_pred_p            = mean_row$p.value,
        diff_pred_coef         = diff_row$estimate,
        diff_pred_se           = diff_row$std.error,
        diff_pred_p            = diff_row$p.value,
        mean_calibrated        = abs(mean_row$estimate - 1) < 0.3,
        heterogeneity_detected = diff_row$p.value < 0.05
      )
    } else {
      tibble(
        mean_pred_coef         = cal$coefficients[["mean.forest.prediction"]],
        mean_pred_se           = cal$std.error[["mean.forest.prediction"]],
        mean_pred_p            = cal$p.value[["mean.forest.prediction"]],
        diff_pred_coef         = cal$coefficients[["differential.forest.prediction"]],
        diff_pred_se           = cal$std.error[["differential.forest.prediction"]],
        diff_pred_p            = cal$p.value[["differential.forest.prediction"]],
        mean_calibrated        = abs(cal$coefficients[["mean.forest.prediction"]] - 1) < 0.3,
        heterogeneity_detected = cal$p.value[["differential.forest.prediction"]] < 0.05
      )
    }
  } else {
    tibble(
      mean_pred_coef = NA_real_, mean_pred_se = NA_real_, mean_pred_p = NA_real_,
      diff_pred_coef = NA_real_, diff_pred_se = NA_real_, diff_pred_p = NA_real_,
      mean_calibrated = NA, heterogeneity_detected = NA
    )
  }
  
  # ----------------------------------------------------------
  #  W.HAT dist
  # For a continuous treatment this is the conditional mean of W
  # given covariates. We want it to:
  #   (a) not pile up at the extremes of the [1,5] scale
  #   (b) have spread across covariate subgroups (not all units
  #       predicted to be the same treatment level)
  # ----------------------------------------------------------
  what      <- forest$W.hat
  what_sd   <- sd(what)
  what_min  <- min(what)
  what_max  <- max(what)
  # Piling flag: if >20% of units have W.hat within 0.25 of
  # the scale endpoints (1 or 5), extrapolation risk is high
  piling_low  <- mean(what < (min(W) + 0.25))
  piling_high <- mean(what > (max(W) - 0.25))
  piling_flag <- piling_low > 0.20 | piling_high > 0.20
  
  # ----------------------------------------------------------
  #  support check 
  # Within each quintile of W.hat, does the observed treatment
  # W vary? 
  # ----------------------------------------------------------
  support_df <- tibble(W = W, W_hat = what) %>%
    mutate(what_quintile = ntile(W_hat, 5)) %>%
    group_by(what_quintile) %>%
    summarise(
      W_mean   = mean(W),
      W_sd     = sd(W),
      W_range  = max(W) - min(W),
      n        = n(),
      # Flag if treatment SD within quintile is < 0.5
      # (very little variation to identify effect)
      low_variation = sd(W) < 0.5,
      .groups  = "drop"
    ) %>%
    mutate(outcome = y, treatment = t)
  
  list(
    outcome    = y,
    treatment  = t,
    n          = nrow(d),
    cal        = cal_df,
    what       = tibble(W = W, W_hat = what,
                        outcome = y, treatment = t),
    what_stats = tibble(
      outcome      = y,
      treatment    = t,
      what_sd      = what_sd,
      what_min     = what_min,
      what_max     = what_max,
      piling_low   = piling_low,
      piling_high  = piling_high,
      piling_flag  = piling_flag
    ),
    support    = support_df
  )
}

# ============================================================
# running 
# ============================================================

grid <- expand_grid(
  outcome   = outcomes_all,
  treatment = treats_cont
)

cat("Running continuous CF diagnostics on", nrow(grid),
    "treatment-outcome combinations...\n")
cat("Using first MICE imputation (df1) for diagnostics.\n\n")

diag_results <- pmap(grid, function(outcome, treatment) {
  cat(" ", treatment, "->", outcome, "\n")
  cf_cont_diag_fn(df1, outcome, treatment, controls_cf_raw)
})

# Remove any NULLs from skipped combinations
diag_results <- compact(diag_results)

# ============================================================
# sum table
# ============================================================

cal_summary <- map_dfr(diag_results, function(r) {
  bind_cols(
    tibble(
      outcome           = r$outcome,
      treatment         = r$treatment,
      n                 = r$n,
      outcome_label     = outcome_labels[r$outcome],
      treatment_label   = str_remove(
        treatment_labels[r$treatment], " \\(continuous\\)"
      )
    ),
    r$cal,
    r$what_stats %>% select(-outcome, -treatment)
  )
})

# Flag overall quality
cal_summary <- cal_summary %>%
  mutate(
    calibration_quality = case_when(
      is.na(mean_calibrated)     ~ "Unknown",
      mean_calibrated & !piling_flag ~ "Good",
      mean_calibrated &  piling_flag ~ "Good (overlap caution)",
      !mean_calibrated & !piling_flag ~ "Recalibrate",
      TRUE                            ~ "Recalibrate + overlap caution"
    ),
    primary = outcome == "O1_rr_index"
  )

cat("\n=== CALIBRATION SUMMARY ===\n")
cat("Primary outcome:\n")
cal_summary %>%
  filter(primary) %>%
  select(treatment_label, mean_pred_coef, mean_pred_p,
         diff_pred_coef, diff_pred_p,
         piling_flag, calibration_quality) %>%
  print(n = Inf)

cat("\nAll outcomes:\n")
cal_summary %>%
  select(treatment_label, outcome_label,
         mean_pred_coef, mean_pred_p,
         piling_flag, calibration_quality) %>%
  print(n = Inf)

write_csv(
  cal_summary,
  "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_cf_cont_calibration.csv"
)
cat("\nSaved: mpsa_cf_cont_calibration.csv\n")

# ============================================================
# figure 1
# ============================================================

what_all <- map_dfr(diag_results, ~ .$what) %>%
  mutate(
    treatment_label = str_remove(
      treatment_labels[treatment], " \\(continuous\\)"
    ),
    outcome_label = outcome_labels[outcome]
  )

# One panel per treatment (W.hat distribution is outcome-invariant
# since it only depends on X and W, not Y — so just show primary)
what_primary <- what_all %>%
  filter(outcome == "O1_rr_index")

p_what <- ggplot(what_primary,
                 aes(x = W_hat, fill = treatment_label,
                     color = treatment_label)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_rug(alpha = 0.05, linewidth = 0.3) +
  # Observed treatment distribution for reference
  geom_density(aes(x = W), linetype = "dashed",
               fill = NA, color = "gray40", linewidth = 0.6) +
  scale_x_continuous(breaks = 1:5,
                     labels = c("1\nMuch better", "2", "3\nNeutral",
                                "4", "5\nMuch worse")) +
  facet_wrap(~ treatment_label, ncol = 2,
             labeller = label_wrap_gen(width = 30)) +
  scale_fill_brewer(palette = "Set2")  +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "Continuous Treatment Propensity: W.hat Distribution",
    subtitle = paste(
      "Solid = E[W | X] (propensity model); Dashed = observed W distribution",
      "Primary outcome (racial resentment) shown; W.hat is outcome-invariant",
      sep = "\n"
    ),
    x       = "Treatment value (economic pessimism scale)",
    y       = "Density",
    caption = paste(
      "W.hat = forest's estimate of E[W | X]. Good support: spread across [1,5],",
      "not piling at endpoints. Dashed line = observed treatment distribution.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption     = element_text(size = 8, color = "gray50",
                                    hjust = 0, lineheight = 1.3)
  )

ggsave(
  "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_cf_cont_what_density.png",
  p_what, width = 10, height = 7, dpi = 150
)
cat("Saved: mpsa_cf_cont_what_density.png\n")

# ============================================================
# fig 2
# ============================================================

support_all <- map_dfr(diag_results, ~ .$support) %>%
  mutate(
    treatment_label = str_remove(
      treatment_labels[treatment], " \\(continuous\\)"
    ),
    outcome_label = outcome_labels[outcome],
    quintile_label = paste0("Q", what_quintile)
  ) %>%
  filter(outcome == "O1_rr_index")   # outcome-invariant, show primary

p_support <- ggplot(support_all,
                    aes(x = quintile_label, y = W_sd,
                        fill = low_variation)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  scale_fill_manual(
    values = c("FALSE" = "#2ecc71", "TRUE" = "#e74c3c"),
    labels = c("FALSE" = "Adequate variation",
               "TRUE"  = "Low variation (< 0.5 SD)")
  ) +
  facet_wrap(~ treatment_label, ncol = 2,
             labeller = label_wrap_gen(width = 30)) +
  labs(
    title    = "Local Treatment Support by W.hat Quintile",
    subtitle = paste(
      "SD of observed W within each quintile of E[W | X]",
      "Low SD = limited treatment variation for identification in that region",
      sep = "\n"
    ),
    x       = "W.hat quintile (Q1 = lowest predicted treatment)",
    y       = "SD of observed W within quintile",
    fill    = NULL,
    caption = paste(
      "Dashed line = SD threshold of 0.5.",
      "Red bars indicate regions where the forest has limited variation",
      "in the treatment to identify the causal effect.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption     = element_text(size = 8, color = "gray50",
                                    hjust = 0, lineheight = 1.3)
  )

ggsave(
  "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_cf_cont_support_check.png",
  p_support, width = 10, height = 7, dpi = 150
)
cat("Saved: mpsa_cf_cont_support_check.png\n")

# ============================================================
# key items
# ============================================================

cat("\n=== KEY DIAGNOSTIC TAKEAWAYS ===\n")

cat("\nCalibration (mean coef should be ~ 1.0):\n")
cal_summary %>%
  filter(primary) %>%
  mutate(
    deviation = round(abs(mean_pred_coef - 1), 3),
    flag = ifelse(deviation > 0.3, "FLAG", "ok")
  ) %>%
  select(treatment_label, mean_pred_coef, deviation, flag) %>%
  print()

cat("\nHeterogeneity detected (diff pred coef sig at p < .05):\n")
cal_summary %>%
  filter(primary) %>%
  select(treatment_label, diff_pred_coef, diff_pred_p,
         heterogeneity_detected) %>%
  print()

cat("\nPiling / support issues:\n")
cal_summary %>%
  filter(primary) %>%
  select(treatment_label, what_sd, piling_low, piling_high,
         piling_flag) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  print()

cat("\nOverall quality ratings:\n")
cal_summary %>%
  filter(primary) %>%
  select(treatment_label, calibration_quality) %>%
  print()

cat("\nDone. All diagnostics saved.\n")