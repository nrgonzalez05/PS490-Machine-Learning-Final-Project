# ============================================================
# mpsa_3_cf_main.R
# Causal Forest ATT — main results (cont treatments)
#
# REQUIRES: mpsa_environment.RData (from mpsa_1_data_prep.R)
#
# OUTPUTS:
#   mpsa_cf_cont_results.csv   — continuous ATE results
#   mpsa_cf_cont_coefplot.png  — continuous ATE figure 
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(grf)
  library(ggplot2)
  library(patchwork)
})

select <- dplyr::select
load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")
cat("Outcomes:", length(outcomes_all), "| Primary: O1_rr_index\n")

# ============================================================
# CF functions 
# ============================================================

cf_att_fn <- function(df, y, t, ctrls, seed = 123, min_n = 400) {
  d <- df %>% select(all_of(c(y, t, ctrls))) %>% drop_na()
  if (nrow(d) < min_n || length(unique(d[[t]])) < 2)
    return(tibble(outcome = y, treatment = t, model = "CF(ATT)",
                  est = NA_real_, se = NA_real_, p = NA_real_,
                  n_used = nrow(d), ok = FALSE))
  tryCatch({
    Y <- d[[y]]; W <- d[[t]]
    X <- drop_zero_var(make_X(d, ctrls))
    set.seed(seed)
    f <- causal_forest(X = X, Y = Y, W = W,
                       num.trees = 2000, min.node.size = 5,
                       honesty = TRUE, seed = seed)
    a <- average_treatment_effect(f, target.sample = "treated")
    e <- as.numeric(a["estimate"]); s <- as.numeric(a["std.err"])
    tibble(outcome = y, treatment = t, model = "CF(ATT)",
           est = e, se = s, p = 2 * (1 - pnorm(abs(e / s))),
           n_used = nrow(d), ok = TRUE)
  }, error = function(e)
    tibble(outcome = y, treatment = t, model = "CF(ATT)",
           est = NA_real_, se = NA_real_, p = NA_real_,
           n_used = nrow(d), ok = FALSE))
}

cf_cont_fn <- function(df, y, t, ctrls, seed = 123, min_n = 500) {
  d <- df %>% select(all_of(c(y, t, ctrls))) %>% drop_na()
  if (nrow(d) < min_n)
    return(tibble(outcome = y, treatment = t, model = "CF(cont ATE)",
                  est = NA_real_, se = NA_real_, p = NA_real_,
                  n_used = nrow(d), ok = FALSE))
  tryCatch({
    Y <- d[[y]]; W <- as.numeric(d[[t]])
    X <- drop_zero_var(make_X(d, ctrls))
    set.seed(seed)
    f <- causal_forest(X = X, Y = Y, W = W,
                       num.trees = 2000, min.node.size = 5,
                       honesty = TRUE, seed = seed)
    a <- average_treatment_effect(f)
    e <- as.numeric(a["estimate"]); s <- as.numeric(a["std.err"])
    tibble(outcome = y, treatment = t, model = "CF(cont ATE)",
           est = e, se = s, p = 2 * (1 - pnorm(abs(e / s))),
           n_used = nrow(d), ok = TRUE)
  }, error = function(e)
    tibble(outcome = y, treatment = t, model = "CF(cont ATE)",
           est = NA_real_, se = NA_real_, p = NA_real_,
           n_used = nrow(d), ok = FALSE))
}

pool_cf <- function(imp_list, grid, fn, ctrls) {
  map_dfr(imp_list, function(df) {
    pmap_dfr(grid, \(outcome, treatment) fn(df, outcome, treatment, ctrls))
  }, .id = "imp") %>%
    group_by(model, outcome, treatment) %>%
    summarise(
      est    = mean(est, na.rm = TRUE),
      se     = sqrt(mean(se^2, na.rm = TRUE)),
      p      = 2 * (1 - pnorm(abs(est / se))),
      n_used = NA_integer_,
      ok     = TRUE,
      .groups = "drop"
    ) %>%
    mutate(
      stars           = stars_p(p),
      p_fmt           = fmt_p(p),
      outcome_label   = outcome_labels[outcome],
      treatment_label = treatment_labels[treatment],
      primary         = outcome == "O1_rr_index"
    )
}

# ============================================================
# run CF continous ATE 
# ============================================================

cat("\nRunning CF(cont ATE) — appendix...\n")
cf_cont_pooled <- pool_cf(imp_list, grid_cont, cf_cont_fn, controls_cf_raw)

cat("\n=== CF CONTINUOUS ATE POOLED RESULTS ===\n")
print(cf_cont_pooled %>%
        select(treatment_label, outcome_label, est, se, p_fmt, stars),
      n = Inf)

write_csv(cf_cont_pooled,
          "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_cf_cont_results.csv")
cat("Saved: mpsa_cf_cont_results.csv\n")

# ============================================================
#  plot 
# ============================================================

make_coefplot <- function(results_df, title_str, subtitle_str,
                           caption_str = "") {
  plot_df <- results_df %>%
    filter(ok) %>%
    mutate(
      ci_lo           = est - 1.96 * se,
      ci_hi           = est + 1.96 * se,
      sig             = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
      outcome_label   = factor(outcome_label,
                                levels = rev(unname(outcome_labels))),
      treatment_label = factor(treatment_label,
                                levels = unname(treatment_labels[
                                  treatment_labels %in% unique(treatment_label)
                                ]))
    )

  ggplot(plot_df,
         aes(x = est, y = outcome_label, color = sig,
             xmin = ci_lo, xmax = ci_hi)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_errorbarh(height = 0.25, linewidth = 0.7) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c("p < .05" = "#c0392b", "p \u2265 .05" = "#7f8c8d")
    ) +
    facet_wrap(~ treatment_label, ncol = 2) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x        = "Estimated treatment effect (ATT)",
      y        = NULL,
      color    = NULL,
      caption  = caption_str
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position  = "bottom",
          strip.text       = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank())
}

# ============================================================
# SAVE 
# ============================================================


p_cont <- make_coefplot(
  cf_cont_pooled,
  title_str    = "Causal Forest Continuous ATE: Dose-Response of Economic Pessimism",
  subtitle_str = "Continuous treatments (full 1-5 scale) | 5 outcomes | Appendix",
  caption_str  = "Points = ATE slope estimate; error bars = 95% CI. Appendix figure."
)

ggsave("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_cf_cont_coefplot.png",
       p_cont, width = 12, height = 8, dpi = 150)
cat("Saved: mpsa_cf_cont_coefplot.png\n")

