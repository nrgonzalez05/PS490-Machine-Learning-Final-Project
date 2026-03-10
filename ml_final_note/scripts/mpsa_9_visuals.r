# ============================================================
# mpsa_visuals.R
# Consolidated visualization script — MPSA 2026
#
# REQUIRES (in outputs_ml/):
#   mpsa_gps_results.csv
#   mpsa_cf_cont_results.csv
#   mpsa_mlbalance_summary.csv
#   mpsa_sensitivity_summary.csv
#   mpsa_dose_response_data.csv
#   mpsa_cf_cont_calibration.csv
#
# OUTPUTS (to final_outputs/):
#   fig_01_estimator_comparison_primary.png
#   fig_02_estimator_comparison_all.png
#   fig_03_balance_attenuation.png
#   fig_04_dose_response_primary.png
#   fig_05_dose_response_full.png
#   fig_06_sensitivity_gamma.png
#   fig_07_calibration_primary.png
#   fig_08_what_density.png   [NOTE: requires forest objects —
#                              copied from outputs_ml if present]
#   fig_09_what_support.png   [same as above]
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

select <- dplyr::select

base_dir  <- "~/Documents/mpsa_dubois/ml_final_note"
in_dir    <- file.path(base_dir, "outputs_ml")
out_dir   <- file.path(base_dir, "final_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", out_dir, "\n\n")

# ============================================================
# lablez
# ============================================================

treatment_labels_clean <- c(
  "C_nat_econ_last"  = "National retrospective\neconomic pessimism",
  "C_nat_econ_next"  = "National prospective\neconomic pessimism",
  "C_pers_fin_last"  = "Personal financial\nretrospective pessimism",
  "C_fin_nextyear"   = "Personal financial\nprospective worry"
)

treatment_order <- c(
  "National retrospective\neconomic pessimism",
  "Personal financial\nretrospective pessimism",
  "Personal financial\nprospective worry",
  "National prospective\neconomic pessimism"
)

outcome_labels_clean <- c(
  "O1_rr_index"              = "Racial resentment index",
  "O4_welfare_binary"        = "Welfare spending opposition",
  "O5_aidpoor_binary"        = "Aid to poor opposition",
  "O6_status_jobs_threat"    = "White status / jobs threat",
  "O7_gov_favors_blacks_amt" = "Govt. favors Blacks (amount)"
)

outcome_order <- c(
  "Racial resentment index",
  "Aid to poor opposition",
  "White status / jobs threat",
  "Govt. favors Blacks (amount)",
  "Welfare spending opposition"
)


theme_mpsa <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title        = element_text(size = base_size + 2,
                                       face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1,
                                       color = "gray40", hjust = 0),
      plot.caption      = element_text(size = base_size - 2,
                                       color = "gray50", hjust = 0,
                                       lineheight = 1.3),
      strip.text        = element_text(size = base_size - 1,
                                       face = "bold"),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = base_size - 1)
    )
}

sig_colors <- c("p < .05" = "#c0392b", "p \u2265 .05" = "#7f8c8d")

clean_treat <- function(x) {
  x %>%
    str_remove(" \\(continuous\\)") %>%
    str_replace("National retrospective economic pessimism",
                "National retrospective\neconomic pessimism") %>%
    str_replace("National prospective economic pessimism",
                "National prospective\neconomic pessimism") %>%
    str_replace("Personal financial retrospective pessimism",
                "Personal financial\nretrospective pessimism") %>%
    str_replace("Personal financial prospective worry",
                "Personal financial\nprospective worry")
}

# ============================================================
# data
# ============================================================

cat("Loading CSVs...\n")

gps     <- read_csv(file.path(in_dir, "mpsa_gps_results.csv"),
                    show_col_types = FALSE) %>%
  select(-any_of(c("p_fmt", "stars", "model"))) %>%
  mutate(
    estimator       = "GPS-OLS",
    treatment_label = clean_treat(treatment_labels[treatment]),
    outcome_label   = outcome_labels_clean[outcome],
    ci_lo           = est - 1.96 * se,
    ci_hi           = est + 1.96 * se,
    sig             = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
    stars           = case_when(
      p < .001 ~ "***", p < .01 ~ "**",
      p < .05  ~ "*",   p < .1  ~ ".",
      TRUE     ~ ""
    ),
    primary         = outcome == "O1_rr_index"
  )

cf_cont <- read_csv(file.path(in_dir, "mpsa_cf_cont_results.csv"),
                    show_col_types = FALSE) %>%
  select(-any_of(c("p_fmt", "stars", "model"))) %>%
  mutate(
    estimator       = "CF Continuous",
    treatment_label = clean_treat(treatment_labels[treatment]),
    outcome_label   = outcome_labels_clean[outcome],
    ci_lo           = est - 1.96 * se,
    ci_hi           = est + 1.96 * se,
    sig             = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
    stars           = case_when(
      p < .001 ~ "***", p < .01 ~ "**",
      p < .05  ~ "*",   p < .1  ~ ".",
      TRUE     ~ ""
    ),
    primary         = outcome == "O1_rr_index"
  )


tryCatch({
  load(file.path(in_dir, "mpsa_environment.RData"))
  cat("Environment loaded for label lookups.\n")
}, error = function(e) {
  cat("Environment not loaded — using inline label map.\n")
  treatment_labels <<- c(
    "T7_nat_econ_last_worse"  = "National retrospective economic pessimism",
    "T6_nat_econ_next_worse"  = "National prospective economic pessimism",
    "T8_pers_fin_last_worse"  = "Personal financial retrospective pessimism",
    "T3_fin_nextyear_worried" = "Personal financial prospective worry",
    "C_nat_econ_last"         = "National retrospective economic pessimism (continuous)",
    "C_nat_econ_next"         = "National prospective economic pessimism (continuous)",
    "C_pers_fin_last"         = "Personal financial retrospective pessimism (continuous)",
    "C_fin_nextyear"          = "Personal financial prospective worry (continuous)"
  )
  outcome_labels <<- c(
    "O1_rr_index"              = "Racial resentment index",
    "O4_welfare_binary"        = "Welfare spending opposition",
    "O5_aidpoor_binary"        = "Aid to poor opposition",
    "O6_status_jobs_threat"    = "White status / jobs threat",
    "O7_gov_favors_blacks_amt" = "Government favors Blacks (amount)"
  )
})


gps <- read_csv(file.path(in_dir, "mpsa_gps_results.csv"),
                show_col_types = FALSE) %>%
  select(-any_of(c("p_fmt", "stars", "model"))) %>%
  mutate(
    estimator       = "GPS-OLS",
    treatment_label = clean_treat(treatment_labels[treatment]),
    outcome_label   = outcome_labels_clean[outcome],
    ci_lo           = est - 1.96 * se,
    ci_hi           = est + 1.96 * se,
    sig             = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
    stars           = case_when(
      p < .001 ~ "***", p < .01 ~ "**",
      p < .05  ~ "*",   p < .1  ~ ".",
      TRUE     ~ ""
    ),
    primary         = outcome == "O1_rr_index"
  )

cf_cont <- read_csv(file.path(in_dir, "mpsa_cf_cont_results.csv"),
                    show_col_types = FALSE) %>%
  select(-any_of(c("p_fmt", "stars", "model"))) %>%
  mutate(
    estimator       = "CF Continuous",
    treatment_label = clean_treat(treatment_labels[treatment]),
    outcome_label   = outcome_labels_clean[outcome],
    ci_lo           = est - 1.96 * se,
    ci_hi           = est + 1.96 * se,
    sig             = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
    stars           = case_when(
      p < .001 ~ "***", p < .01 ~ "**",
      p < .05  ~ "*",   p < .1  ~ ".",
      TRUE     ~ ""
    ),
    primary         = outcome == "O1_rr_index"
  )

balance  <- read_csv(file.path(in_dir, "mpsa_mlbalance_summary.csv"),
                     show_col_types = FALSE)
sens     <- read_csv(file.path(in_dir, "mpsa_sensitivity_summary.csv"),
                     show_col_types = FALSE)
dr_data  <- read_csv(file.path(in_dir, "mpsa_dose_response_data.csv"),
                     show_col_types = FALSE)
calib    <- read_csv(file.path(in_dir, "mpsa_cf_cont_calibration.csv"),
                     show_col_types = FALSE) %>%
  mutate(treatment_label = str_remove(treatment_label, " \\(continuous\\)"))

cat("All CSVs loaded.\n\n")

# ============================================================
# FIG 01: ESTIMATOR COMPARISON — PRIMARY OUTCOME
# Side-by-side GPS-OLS vs CF Continuous for racial resentment
# ============================================================

cat("Building Fig 01: Estimator comparison (primary)...\n")

comp_primary <- bind_rows(
  gps     %>% filter(primary),
  cf_cont %>% filter(primary)
) %>%
  mutate(
    estimator       = factor(estimator,
                             levels = c("GPS-OLS", "CF Continuous")),
    treatment_label = factor(treatment_label, levels = treatment_order),
    outcome_label   = factor(outcome_label,   levels = outcome_order)
  )

fig01 <- ggplot(comp_primary,
                aes(x = est, xmin = ci_lo, xmax = ci_hi,
                    y = estimator, color = sig, shape = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.3, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = sig_colors, name = NULL) +
  scale_shape_manual(values = c("GPS-OLS" = 16, "CF Continuous" = 17),
                     name = NULL) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  facet_wrap(~ treatment_label, ncol = 2,
             labeller = label_wrap_gen(width = 35)) +
  labs(
    title    = "GPS-OLS vs. Continuous Causal Forest: Primary Outcome",
    subtitle = "Outcome: Racial resentment index | ATE on continuous 1\u20135 treatment scale | 95% CI",
    x        = "Estimated ATE (slope on 1\u20135 scale)",
    y        = NULL,
    caption  = paste(
      "GPS-OLS uses OLS control set (8 variables); CF uses expanded set (40+ variables).",
      "Both estimators target the continuous ATE. Red = p < .05.",
      sep = "\n"
    )
  ) +
  theme_mpsa() +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "fig_01_estimator_comparison_primary.png"),
       fig01, width = 11, height = 6, dpi = 200)
cat("  Saved: fig_01_estimator_comparison_primary.png\n")

# ============================================================
# FIG 02: ESTIMATOR COMPARISON — ALL OUTCOMES
# GPS-OLS vs CF Continuous across all 5 outcomes
# Outcomes as rows, treatments as facets
# ============================================================

cat("Building Fig 02: Estimator comparison (all outcomes)...\n")

comp_all <- bind_rows(
  gps,
  cf_cont %>% filter(outcome != "O4_welfare_binary")  # welfare appendix
) %>%
  mutate(
    estimator       = factor(estimator,
                             levels = c("GPS-OLS", "CF Continuous")),
    treatment_label = factor(treatment_label, levels = treatment_order),
    outcome_label   = factor(outcome_label,   levels = outcome_order),
    # Flag primary outcome for visual emphasis
    primary_label   = ifelse(primary,
                             paste0(outcome_label, " \u2605"),
                             as.character(outcome_label))
  )

fig02 <- ggplot(comp_all,
                aes(x = est, xmin = ci_lo, xmax = ci_hi,
                    y = estimator, color = sig, shape = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_errorbarh(height = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = sig_colors, name = NULL) +
  scale_shape_manual(values = c("GPS-OLS" = 16, "CF Continuous" = 17),
                     name = NULL) +
  facet_grid(outcome_label ~ treatment_label,
             labeller = labeller(
               treatment_label = label_wrap_gen(20),
               outcome_label   = label_wrap_gen(22)
             )) +
  labs(
    title    = "GPS-OLS vs. Continuous CF: All Outcomes",
    subtitle = "ATE on continuous 1\u20135 treatment scale | 95% CI | \u2605 = primary outcome",
    x        = "Estimated ATE",
    y        = NULL,
    caption  = paste(
      "Welfare spending opposition excluded from CF column (calibration failure — see appendix).",
      "Red = p < .05. GPS-OLS uses OLS control set; CF uses expanded 40+ variable set.",
      sep = "\n"
    )
  ) +
  theme_mpsa(base_size = 9) +
  theme(
    strip.text.y    = element_text(angle = 0, size = 8),
    axis.text.y     = element_text(size = 8),
    panel.spacing   = unit(0.5, "lines")
  )

ggsave(file.path(out_dir, "fig_02_estimator_comparison_all.png"),
       fig02, width = 13, height = 10, dpi = 200)
cat("  Saved: fig_02_estimator_comparison_all.png\n")

# ============================================================
# FIG 03: BALANCE ATTENUATION
# DIM → IPW → DR under OLS vs CF controls
# Panel per treatment, primary outcome only
# ============================================================

cat("Building Fig 03: Balance attenuation...\n")

bal_primary <- balance %>%
  filter(outcome == "Racial resentment index") %>%
  select(treatment, controls, est_dim, est_ipw, est_dr) %>%
  pivot_longer(cols = c(est_dim, est_ipw, est_dr),
               names_to = "estimator_type",
               values_to = "estimate") %>%
  mutate(
    estimator_type = factor(estimator_type,
                            levels = c("est_dim", "est_ipw", "est_dr"),
                            labels = c("DIM\n(no adjustment)",
                                       "IPW\n(propensity weighted)",
                                       "DR\n(doubly robust)")),
    controls = factor(controls,
                      levels = c("OLS (small)", "CF (expanded)"),
                      labels = c("OLS controls\n(8 variables)",
                                 "CF controls\n(40+ variables)")),
    treatment_label = clean_treat(treatment)
  )

fig03 <- ggplot(bal_primary,
                aes(x = estimator_type, y = estimate,
                    color = controls, group = controls)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray60", linewidth = 0.5) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3.5) +
  scale_color_manual(
    values = c("OLS controls\n(8 variables)"  = "#e67e22",
               "CF controls\n(40+ variables)" = "#2980b9"),
    name = "Control set"
  ) +
  facet_wrap(~ treatment_label, ncol = 2,
             labeller = label_wrap_gen(width = 32)) +
  labs(
    title    = "Confounding Attenuation: DIM \u2192 IPW \u2192 Doubly-Robust",
    subtitle = "Primary outcome (racial resentment) | OLS vs. CF expanded control sets",
    x        = "Adjustment stage",
    y        = "Estimated effect",
    caption  = paste(
      "DIM = raw difference in means; IPW = inverse probability weighted;",
      "DR = doubly robust. Steeper decline = more confounding absorbed.",
      "OLS controls absorb ~67% of raw confounding; CF controls absorb ~93%.",
      sep = "\n"
    )
  ) +
  theme_mpsa()

ggsave(file.path(out_dir, "fig_03_balance_attenuation.png"),
       fig03, width = 11, height = 7, dpi = 200)
cat("  Saved: fig_03_balance_attenuation.png\n")


# ============================================================
# FIG 06: ROSENBAUM SENSITIVITY — CRITICAL GAMMA SUMMARY
# ============================================================

cat("Building Fig 06: Sensitivity (critical gamma)...\n")

sens_plot <- sens %>%
  filter(outcome %in% c("Racial resentment index",
                        "Aid to poor opposition",
                        "White status / jobs threat",
                        "Government favors Blacks (amount)")) %>%
  mutate(
    treatment_label = clean_treat(treatment),
    treatment_label = factor(treatment_label, levels = treatment_order),
    outcome_label   = factor(outcome, levels = c(
      "Racial resentment index",
      "Aid to poor opposition",
      "White status / jobs threat",
      "Government favors Blacks (amount)"
    )),
    robust     = as.logical(robust),
    bar_color  = case_when(
      critical_gamma >= 2.5 ~ "#27ae60",
      critical_gamma >= 1.5 ~ "#f39c12",
      TRUE                  ~ "#e74c3c"
    ),
    bar_label  = paste0("\u0393 = ", critical_gamma)
  )

fig06 <- ggplot(sens_plot,
                aes(x = treatment_label, y = critical_gamma,
                    fill = bar_color, label = bar_label)) +
  geom_col(width = 0.65) +
  geom_text(vjust = -0.4, size = 3.2, fontface = "bold") +
  geom_hline(yintercept = c(1.5, 2.0, 3.0),
             linetype = "dashed", color = "gray60",
             linewidth = 0.5) +
  annotate("text", x = 0.4, y = 1.55,
           label = "\u0393 = 1.5\n(moderate)", size = 2.8,
           color = "gray40", hjust = 0) +
  annotate("text", x = 0.4, y = 2.55,
           label = "\u0393 = 2.5\n(strong)", size = 2.8,
           color = "gray40", hjust = 0) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 3.4),
                     breaks = seq(0, 3, by = 0.5)) +
  facet_wrap(~ outcome_label, ncol = 2,
             labeller = label_wrap_gen(width = 28)) +
  labs(
    title    = "Rosenbaum Sensitivity Analysis: Critical \u0393 by Treatment and Outcome",
    subtitle = "Critical \u0393 = unmeasured confounding required to explain away the result",
    x        = NULL,
    y        = "Critical \u0393",
    caption  = paste(
      "Green = \u0393 \u2265 2.5 (strong robustness); Orange = 1.5\u20132.5 (moderate); Red = < 1.5 (weak).",
      "\u0393 = 3.0 indicates the result survives even when an unmeasured confounder triples treatment odds.",
      sep = "\n"
    )
  ) +
  theme_mpsa() +
  theme(
    axis.text.x = element_text(size = 8, angle = 15, hjust = 1),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(file.path(out_dir, "fig_06_sensitivity_gamma.png"),
       fig06, width = 11, height = 8, dpi = 200)
cat("  Saved: fig_06_sensitivity_gamma.png\n")

# ============================================================
# FIG 07: CALIBRATION SUMMARY — PRIMARY OUTCOME
# Mean calibration coefficient ± SE, colored by quality
# ============================================================

cat("Building Fig 07: Calibration summary...\n")

calib_primary <- calib %>%
  filter(primary == TRUE) %>%
  mutate(
    treatment_label = factor(
      treatment_label,
      levels = c(
        "National retrospective economic pessimism",
        "Personal financial retrospective pessimism",
        "Personal financial prospective worry",
        "National prospective economic pessimism"
      )
    ),
    ci_lo = mean_pred_coef - 1.96 * mean_pred_se,
    ci_hi = mean_pred_coef + 1.96 * mean_pred_se,
    quality_color = case_when(
      calibration_quality == "Good"        ~ "#27ae60",
      calibration_quality == "Recalibrate" ~ "#e74c3c",
      TRUE                                 ~ "#f39c12"
    )
  )

fig07 <- ggplot(calib_primary,
                aes(x = mean_pred_coef, xmin = ci_lo, xmax = ci_hi,
                    y = fct_rev(treatment_label),
                    color = quality_color)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  geom_vline(xintercept = c(0.7, 1.3), linetype = "dotted",
             color = "gray60", linewidth = 0.5) +
  geom_errorbarh(height = 0.25, linewidth = 1.0) +
  geom_point(size = 4) +
  geom_text(aes(label = round(mean_pred_coef, 3)),
            hjust = -0.35, size = 3.5, fontface = "bold") +
  scale_color_identity() +
  scale_x_continuous(limits = c(-5, 6),
                     breaks = c(0, 0.7, 1, 1.3, 2, 3)) +
  annotate("rect", xmin = 0.7, xmax = 1.3,
           ymin = 0.4, ymax = 4.6,
           alpha = 0.07, fill = "#27ae60") +
  labs(
    title    = "Continuous CF Calibration: Mean Forest Prediction Coefficient",
    subtitle = "Primary outcome (racial resentment) | Coefficient should be \u2248 1.0 with tight SE",
    x        = "Mean forest prediction coefficient (95% CI)",
    y        = NULL,
    caption  = paste(
      "Green shaded band = acceptable calibration zone [0.7, 1.3].",
      "Dashed line = ideal value of 1.0.",
      "Green points = Good calibration; Red = Recalibrate (insufficient identifying power).",
      "Note: CI width reflects identifying power — wide CI = low W.hat variation.",
      sep = "\n"
    )
  ) +
  theme_mpsa() +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(out_dir, "fig_07_calibration_primary.png"),
       fig07, width = 10, height = 5, dpi = 200)
cat("  Saved: fig_07_calibration_primary.png\n")

# ============================================================
# FIG 08 / 09: W.HAT DENSITY + SUPPORT CHECK
# These figures require forest objects to re-generate.
# If they already exist in outputs_ml/, copy them to final_outputs/.
# ============================================================

cat("Checking for W.hat figures...\n")

what_src     <- file.path(in_dir, "mpsa_cf_cont_what_density.png")
support_src  <- file.path(in_dir, "mpsa_cf_cont_support_check.png")

if (file.exists(what_src)) {
  file.copy(what_src,
            file.path(out_dir, "fig_08_what_density.png"),
            overwrite = TRUE)
  cat("  Copied: fig_08_what_density.png\n")
} else {
  cat("  NOTE: fig_08_what_density.png not found in outputs_ml/.",
      "Run mpsa_8_diag.R first.\n")
}

if (file.exists(support_src)) {
  file.copy(support_src,
            file.path(out_dir, "fig_09_what_support.png"),
            overwrite = TRUE)
  cat("  Copied: fig_09_what_support.png\n")
} else {
  cat("  NOTE: fig_09_what_support.png not found in outputs_ml/.",
      "Run mpsa_8_diag.R first.\n")
}

