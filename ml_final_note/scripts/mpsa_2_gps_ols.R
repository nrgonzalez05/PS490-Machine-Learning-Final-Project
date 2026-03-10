# ============================================================
# mpsa_2_gps_ols.R

# REQUIRES: mpsa_environment.RData (from mpsa_1_data_prep.R)

# OUTPUTS:
#   mpsa_gps_results.csv       — full results table
#   mpsa_gps_coefplot.png      — coefficient plot (appendix fig)

# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

select <- dplyr::select
load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")
cat("Outcomes:", length(outcomes_all), "| Treatments:", length(treats_cont), "\n")

# ============================================================
# GPS -> OLS 
# ============================================================

gps_ols <- function(df, y, t, ctrls, min_n = 400) {

  d <- df %>% select(all_of(c(y, t, ctrls))) %>% drop_na()

  if (nrow(d) < min_n)
    return(tibble(outcome = y, treatment = t,
                  est = NA_real_, se = NA_real_, p = NA_real_,
                  n_used = nrow(d), ok = FALSE, reason = "too_few_rows"))

  tryCatch({
    w_form <- as.formula(paste(t, "~", paste(ctrls, collapse = " + ")))
    w_fit  <- lm(w_form, data = d)
    mu     <- predict(w_fit)
    sigma  <- sd(resid(w_fit), na.rm = TRUE)

    if (is.na(sigma) || sigma <= 0) stop("sigma_nonpositive")

    W   <- d[[t]]
    gps <- dnorm(W, mean = mu, sd = sigma)
    gps <- pmin(pmax(gps,
                     quantile(gps, 0.01, na.rm = TRUE)),
                     quantile(gps, 0.99, na.rm = TRUE))

    d2     <- d %>% mutate(.W = W, .gps = gps)
    y_form <- as.formula(paste(
      y, "~ .W + .gps + I(.W * .gps) +", paste(ctrls, collapse = " + ")
    ))
    fit <- lm(y_form, data = d2)
    sm  <- summary(fit)$coefficients

    tibble(outcome = y, treatment = t,
           est    = sm[".W", "Estimate"],
           se     = sm[".W", "Std. Error"],
           p      = sm[".W", "Pr(>|t|)"],
           n_used = nrow(d), ok = TRUE, reason = "")

  }, error = function(e)
    tibble(outcome = y, treatment = t,
           est = NA_real_, se = NA_real_, p = NA_real_,
           n_used = nrow(d), ok = FALSE,
           reason = paste0("gps_error: ", conditionMessage(e)))
  )
}

# ============================================================
# run across the mice
# ============================================================

grid_cont <- expand_grid(outcome = outcomes_all, treatment = treats_cont)

cat("Running GPS->OLS across", length(imp_list), "MICE imputations...\n")

gps_imp <- map_dfr(imp_list, function(df) {
  pmap_dfr(grid_cont, \(outcome, treatment)
    gps_ols(df, outcome, treatment, controls_ols_raw))
}, .id = "imp")

gps_pooled <- gps_imp %>%
  group_by(outcome, treatment) %>%
  summarise(
    est    = mean(est, na.rm = TRUE),
    se     = sqrt(mean(se^2, na.rm = TRUE)),
    p      = 2 * (1 - pnorm(abs(est / se))),
    n_used = NA_integer_,
    ok     = TRUE,
    .groups = "drop"
  ) %>%
  mutate(
    model           = "GPS->OLS (ATE, continuous)",
    stars           = stars_p(p),
    p_fmt           = fmt_p(p),
    outcome_label   = outcome_labels[outcome],
    treatment_label = treatment_labels[treatment],
    primary         = outcome == "O1_rr_index"
  )

cat("\n=== GPS->OLS POOLED RESULTS ===\n")
cat("PRIMARY OUTCOME (O1 — racial resentment):\n")
print(gps_pooled %>%
        filter(primary) %>%
        select(treatment_label, outcome_label, est, se, p_fmt, stars),
      n = Inf)
cat("\nSECONDARY OUTCOMES:\n")
print(gps_pooled %>%
        filter(!primary) %>%
        select(treatment_label, outcome_label, est, se, p_fmt, stars),
      n = Inf)

write_csv(gps_pooled,
          "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_gps_results.csv")
cat("\nSaved: mpsa_gps_results.csv\n")

# ============================================================
# plot
# ============================================================

plot_df <- gps_pooled %>%
  filter(ok) %>%
  mutate(
    ci_lo = est - 1.96 * se,
    ci_hi = est + 1.96 * se,
    sig   = ifelse(p < 0.05, "p < .05", "p \u2265 .05"),
    outcome_label = factor(outcome_label,
                           levels = rev(unname(outcome_labels))),
    # Bold label for primary outcome
    outcome_label_display = ifelse(
      outcome == "O1_rr_index",
      paste0(outcome_label, " \u2605"),  # star marker for primary
      as.character(outcome_label)
    )
  )

p <- ggplot(plot_df,
            aes(x = est, y = outcome_label, color = sig,
                xmin = ci_lo, xmax = ci_hi)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(height = 0.25, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("p < .05" = "#c0392b", "p \u2265 .05" = "#7f8c8d")) +
  facet_wrap(~ treatment_label, ncol = 2, scales = "free_x") +
  labs(
    title    = "GPS\u2192OLS: Continuous ATE of Economic Pessimism",
    subtitle = "Pooled across 3 MICE imputations | Controls: OLS set | 5 outcomes (\u2605 = primary)",
    x        = "Estimated ATE (dose-response slope)",
    y        = NULL,
    color    = NULL,
    caption  = "Error bars = 95% CI. Appendix results."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9))

ggsave("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_gps_coefplot.png",
       p, width = 11, height = 7, dpi = 150)
cat("Saved: mpsa_gps_coefplot.png\n")
