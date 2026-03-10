# ============================================================
# mpsa_6_rosenbaum.R
# Rosenbaum bounds sensitivity analysis
#
# Tests how strong unmeasured confounding would need to be to
# explain away main treatment effect estimates across all four
# binary treatments and all 5 analysis outcomes.
#
# OUTPUTS:
#   mpsa_sensitivity_bounds_[treatment]_[outcome].png
#   mpsa_sensitivity_critical_gamma.png
#   mpsa_sensitivity_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rbounds)
  library(MatchIt)
})

select <- dplyr::select
load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")
cat("Outcomes:", length(outcomes_all), "| Primary: O1_rr_index\n")

# ============================================================
# config
# ============================================================

# All 5  outcomes
sens_outcomes <- outcomes_all

sens_treatments <- c(
  "T7_nat_econ_last_worse",
  "T6_nat_econ_next_worse",
  "T8_pers_fin_last_worse",
  "T3_fin_nextyear_worried"
)

gamma_seq <- seq(1.0, 3.0, by = 0.25)

# ============================================================
# HELPER: Run matching + Rosenbaum bounds for one cell
# ============================================================

run_rosenbaum <- function(df, treat_var, outcome_var,
                           ctrls, gamma_seq,
                           treat_label, outcome_label,
                           is_primary = FALSE) {

  cat("\n", strrep("-", 60), "\n")
  cat("Treatment:", treat_label, "\n")
  cat("Outcome:  ", outcome_label,
      if (is_primary) " [PRIMARY]" else "", "\n")
  cat(strrep("-", 60), "\n")

  d <- df %>%
    select(all_of(c(treat_var, outcome_var, ctrls))) %>%
    drop_na()

  cat("n =", nrow(d), "| treated =", sum(d[[treat_var]]), "\n")

  match_formula <- as.formula(
    paste(treat_var, "~", paste(ctrls, collapse = " + "))
  )

  m_out <- tryCatch(
    matchit(match_formula, data = d, method = "nearest",
            distance = "glm", ratio = 1, replace = FALSE),
    error = function(e) {
      cat("Matching failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(m_out)) return(NULL)

  matched_data <- match.data(m_out)
  cat("Matched n =", nrow(matched_data), "\n")

  treated_ids <- matched_data %>% filter(!!sym(treat_var) == 1)
  control_ids <- matched_data %>% filter(!!sym(treat_var) == 0)

  pairs <- inner_join(
    treated_ids %>% select(subclass, y_t = !!sym(outcome_var)),
    control_ids %>% select(subclass, y_c = !!sym(outcome_var)),
    by = "subclass"
  )

  cat("Matched pairs:", nrow(pairs), "\n")

  if (nrow(pairs) < 10) {
    cat("Too few pairs, skipping.\n")
    return(NULL)
  }

  hl_est <- wilcox.test(pairs$y_t, pairs$y_c,
                         paired = TRUE,
                         conf.int = TRUE,
                         exact = FALSE)$estimate

  cat("Hodges-Lehmann estimate:", round(hl_est, 4), "\n")

  # Use psens for continuous, binarysens for binary
  is_binary <- all(d[[outcome_var]] %in% c(0, 1), na.rm = TRUE)

  bounds_result <- tryCatch({
    if (is_binary) {
      binarysens(pairs$y_t, pairs$y_c,
                 Gamma = max(gamma_seq), GammaInc = 0.25)
    } else {
      psens(pairs$y_t, pairs$y_c,
            Gamma = max(gamma_seq), GammaInc = 0.25)
    }
  }, error = function(e) {
    cat("rbounds failed:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(bounds_result)) return(NULL)

  bounds_df <- bounds_result$bounds %>%
    as_tibble() %>%
    mutate(
      treatment   = treat_label,
      outcome     = outcome_label,
      hl_estimate = hl_est,
      is_binary   = is_binary
    )

  critical_gamma <- bounds_df %>%
    filter(`Upper bound` > 0.05) %>%
    slice(1) %>%
    pull(Gamma)

  if (length(critical_gamma) == 0) critical_gamma <- max(gamma_seq)

  cat("Critical Gamma (upper p > 0.05):", critical_gamma, "\n")

  # Individual bounds plot
  plot_df <- bounds_df %>%
    rename(gamma   = Gamma,
           lower_p = `Lower bound`,
           upper_p = `Upper bound`)

  p <- ggplot(plot_df, aes(x = gamma)) +
    geom_ribbon(aes(ymin = lower_p, ymax = upper_p),
                fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = lower_p), color = "steelblue", linewidth = 1) +
    geom_line(aes(y = upper_p), color = "steelblue4", linewidth = 1,
              linetype = "dashed") +
    geom_hline(yintercept = 0.05, color = "red", linetype = "dashed",
               linewidth = 0.8) +
    geom_vline(xintercept = critical_gamma, color = "darkred",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = critical_gamma + 0.05, y = 0.15,
             label = paste0("Critical \u0393 = ", critical_gamma),
             hjust = 0, size = 3.5, color = "darkred") +
    scale_x_continuous(breaks = gamma_seq) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
    labs(
      title    = paste0("Rosenbaum Sensitivity Bounds",
                        if (is_primary) " [PRIMARY OUTCOME]" else ""),
      subtitle = paste0(treat_label, "\nOutcome: ", outcome_label),
      x        = expression(Gamma ~ "(unmeasured confounding odds ratio)"),
      y        = "p-value bounds",
      caption  = paste0(
        "Solid = lower bound; Dashed = upper bound (adverse confounding).\n",
        "Red dashed = \u03b1 = 0.05. Dotted vertical = critical \u0393.\n",
        "Matched pairs: ", nrow(pairs),
        " | HL estimate: ", round(hl_est, 3)
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption  = element_text(size = 8, color = "gray40")
    )

  stem <- paste0(
    gsub("[^a-zA-Z0-9]", "_", gsub("\\s+", "_", trimws(treat_label))),
    "_",
    gsub("[^a-zA-Z0-9]", "_", gsub("\\s+", "_", trimws(outcome_label)))
  )

  fname <- paste0(
    "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/",
    "mpsa_sensitivity_bounds_", stem, ".png"
  )
  ggsave(fname, p, width = 8, height = 5, dpi = 150)
  cat("Saved:", basename(fname), "\n")

  tibble(
    treatment       = treat_label,
    outcome         = outcome_label,
    primary         = is_primary,
    n_matched_pairs = nrow(pairs),
    hl_estimate     = round(hl_est, 4),
    is_binary       = is_binary,
    critical_gamma  = critical_gamma,
    robust          = critical_gamma >= 1.5
  )
}

# ============================================================
# all
# ============================================================

sensitivity_results <- list()

for (tr in sens_treatments) {
  for (out in sens_outcomes) {

    key <- paste0(tr, "_", out)

    res <- run_rosenbaum(
      df            = df1,
      treat_var     = tr,
      outcome_var   = out,
      ctrls         = controls_ols_raw,
      gamma_seq     = gamma_seq,
      treat_label   = treatment_labels[tr],
      outcome_label = outcome_labels[out],
      is_primary    = out == "O1_rr_index"
    )

    if (!is.null(res)) {
      sensitivity_results[[key]] <- res
    }
  }
}

# ============================================================
# sum
# ============================================================

sens_summary <- bind_rows(sensitivity_results) %>%
  arrange(primary %>% desc(), outcome, treatment) %>%
  mutate(
    robustness_rating = case_when(
      critical_gamma >= 2.0  ~ "Strong (Gamma >= 2.0)",
      critical_gamma >= 1.5  ~ "Moderate (Gamma 1.5-2.0)",
      critical_gamma >= 1.25 ~ "Weak (Gamma 1.25-1.5)",
      TRUE                   ~ "Fragile (Gamma < 1.25)"
    )
  )

cat("\n\n", strrep("=", 60), "\n")
cat("ROSENBAUM BOUNDS SENSITIVITY SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("\nPRIMARY OUTCOME:\n")
print(sens_summary %>% filter(primary), n = Inf)
cat("\nSECONDARY OUTCOMES:\n")
print(sens_summary %>% filter(!primary), n = Inf)

write_csv(sens_summary,
          "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_sensitivity_summary.csv")
cat("\nSaved: mpsa_sensitivity_summary.csv\n")

# ============================================================
# plot
# ============================================================

outcome_levels <- unname(outcome_labels)

p_summary <- sens_summary %>%
  mutate(
    treatment = str_wrap(treatment, width = 30),
    outcome   = factor(outcome, levels = outcome_levels),
    primary_label = ifelse(primary, paste0(outcome, "\n[PRIMARY]"), outcome)
  ) %>%
  ggplot(aes(x = treatment, y = critical_gamma,
             fill = robustness_rating)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1.5, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 2.0, linetype = "dotted",
             color = "darkred", linewidth = 0.8) +
  geom_text(aes(label = round(critical_gamma, 2)),
            vjust = -0.4, size = 3) +
  facet_wrap(~ outcome, ncol = 3) +
  scale_fill_manual(values = c(
    "Strong (Gamma >= 2.0)"    = "#2166ac",
    "Moderate (Gamma 1.5-2.0)" = "#74add1",
    "Weak (Gamma 1.25-1.5)"    = "#fdae61",
    "Fragile (Gamma < 1.25)"   = "#d73027"
  )) +
  scale_y_continuous(
    limits = c(0, max(sens_summary$critical_gamma, na.rm = TRUE) * 1.15),
    breaks = c(1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
  ) +
  labs(
    title    = "Rosenbaum Sensitivity: Critical Gamma by Treatment and Outcome",
    subtitle = paste0(
      "Critical Gamma = confounding strength needed to make upper p-bound > 0.05\n",
      "Red dashed = 1.5 (moderate); Red dotted = 2.0 (strong). ",
      "O1 (racial resentment) = primary outcome."
    ),
    x       = NULL,
    y       = expression("Critical " * Gamma),
    fill    = "Robustness",
    caption = "1:1 nearest-neighbor PSM on OLS control set. Wilcoxon signed-rank test."
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 8),
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(
  "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_sensitivity_critical_gamma.png",
  p_summary, width = 14, height = 8, dpi = 150
)
cat("Saved: mpsa_sensitivity_critical_gamma.png\n")

cat("\n", strrep("=", 60), "\n")
cat("INTERPRETATION NOTES\n")
cat(strrep("=", 60), "\n")
cat("
Critical Gamma interpretation:
  Gamma = 1.0  No unmeasured confounding (naive)
  Gamma = 1.25 Weak — small confounder breaks result
  Gamma = 1.5  Moderate — conventional minimum bar
  Gamma = 2.0  Strong — confounder must double treatment odds
  Gamma > 2.0  Very strong — rarely achieved in survey data

Upper bound = conservative case (adverse confounding).
Critical Gamma based on upper bound.

Results with critical Gamma >= 1.5 reported with moderate
confidence. Below 1.5 flagged as fragile in paper.

Outcome hierarchy: O1 (racial resentment) primary.
O4/O5/O6/O7 secondary per pre-specified plan.
")

cat("Done. All outputs saved.\n")
