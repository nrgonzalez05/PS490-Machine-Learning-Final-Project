# ============================================================
# mpsa_4_robustness.R
# Robustness checks:
#   A) CATE heterogeneity — all 5 outcomes x T7 (primary
#      treatment), with partisan media as amplification moderator
#   B) Covariate balance table (cobalt SMD, love plots)
#
# REQUIRES: mpsa_environment.RData (from mpsa_1_data_prep.R)
#
# OUTPUTS:
#   mpsa_cate_[outcome]_nat_retro_[dist/blp/mods].png
#   mpsa_balance_love_[treatment].png
#   mpsa_balance_table.csv
#   mpsa_balance_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(grf)
  library(MatchIt)
  library(cobalt)
  library(ggplot2)
  library(patchwork)
})

select <- dplyr::select
load("~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData")
cat("Environment loaded.\n")

# ============================================================
# CATE HETEROGENEITY
# ============================================================

cat("\n", strrep("=", 60), "\n")
cat("A) CATE HETEROGENEITY (T7 x all 5 outcomes)\n")
cat(strrep("=", 60), "\n\n")

fit_cf_cate <- function(df, y, t, ctrls, seed = 123) {
  extra_cols <- intersect("media_high", names(df))
  d <- df %>% select(all_of(c(y, t, ctrls, extra_cols))) %>% drop_na()
  Y <- d[[y]]; W <- as.numeric(d[[t]])
  X <- drop_zero_var(make_X(d, ctrls))
  set.seed(seed)
  f <- causal_forest(X = X, Y = Y, W = W,
                     num.trees = 2000, min.node.size = 5,
                     honesty = TRUE, seed = seed)
  list(forest = f, data = d, X = X)
}

plot_cate <- function(cf_obj, outcome_label, treatment_label,
                      moderators = c("V241177", "V241567x",
                                     "V241465x", "V241025",
                                     "media_high"),
                      mod_labels = c("Ideology\n(1=lib, 7=con)",
                                     "Household income",
                                     "Education",
                                     "Party registration",
                                     "Partisan media\n(0=low 0-2, 1=high 3+)")) {
  forest  <- cf_obj$forest
  d       <- cf_obj$data
  tau_hat <- predict(forest)$predictions

  # A: CATE distribution
  p_dist <- ggplot(data.frame(tau = tau_hat), aes(x = tau)) +
    geom_histogram(bins = 40, fill = "#2c6fad", color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_vline(xintercept = mean(tau_hat), linetype = "solid",
               color = "darkblue", linewidth = 0.8) +
    labs(title    = paste("CATE distribution:", outcome_label),
         subtitle = paste("Treatment:", treatment_label),
         x = "Individual treatment effect estimate", y = "Count",
         caption = "Red dashed = 0; Blue solid = mean CATE") +
    theme_minimal(base_size = 12)

  # B: Best linear projection
  blp     <- best_linear_projection(forest, cf_obj$X)
  blp_mat <- as.data.frame(unclass(blp))
  blp_mat$covariate <- rownames(blp_mat)
  rownames(blp_mat) <- NULL

  blp_df <- blp_mat %>%
    filter(covariate != "(Intercept)") %>%
    mutate(
      sig       = ifelse(`Pr(>|t|)` < 0.05, "p < .05", "p \u2265 .05"),
      covariate = str_replace_all(covariate, "TRUE|1$", "")
    ) %>%
    arrange(desc(abs(Estimate))) %>%
    slice_head(n = 20)

  p_blp <- ggplot(blp_df,
                  aes(x = reorder(covariate, abs(Estimate)),
                      y = Estimate, fill = sig)) +
    geom_col(alpha = 0.85) +
    geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`,
                      ymax = Estimate + 1.96 * `Std. Error`),
                  width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    scale_fill_manual(
      values = c("p < .05" = "#c0392b", "p \u2265 .05" = "#7f8c8d")
    ) +
    labs(title    = "Best linear projection of CATEs (top 20 covariates)",
         subtitle = paste(outcome_label, "|", treatment_label),
         x = NULL, y = "Coefficient on CATE", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  # C: CATE by moderator smooths
  mod_plots <- map2(moderators, mod_labels, function(mod, lab) {

    if (mod == "media_high") {
      if (!"media_high" %in% names(d)) {
        cat("  NOTE: media_high not available for moderator plot\n")
        return(NULL)
      }
      df_mod <- data.frame(moderator = factor(d$media_high),
                           tau = tau_hat) %>% drop_na()
      return(
        ggplot(df_mod, aes(x = moderator, y = tau, fill = moderator)) +
          geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
          scale_x_discrete(labels = c("0" = "Low\n(0-2 outlets)",
                                      "1" = "High\n(3+ outlets)")) +
          scale_fill_manual(values = c("0" = "#7f8c8d", "1" = "#c0392b"),
                            guide = "none") +
          labs(x = lab, y = "CATE",
               title = "Partisan media amplification moderator",
               subtitle = "High = 3+ right-wing outlets") +
          theme_minimal(base_size = 10)
      )
    }

    if (!mod %in% names(d)) return(NULL)
    df_mod <- data.frame(moderator = d[[mod]], tau = tau_hat) %>% drop_na()
    ggplot(df_mod, aes(x = moderator, y = tau)) +
      geom_point(alpha = 0.15, size = 0.6, color = "#2c6fad") +
      geom_smooth(method = "loess", se = TRUE,
                  color = "#c0392b", linewidth = 1.1,
                  fill = "#f1948a", alpha = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(x = lab, y = "CATE", title = mod) +
      theme_minimal(base_size = 10)
  }) %>% compact()

  list(p_dist = p_dist, p_blp = p_blp, mod_plots = mod_plots)
}

save_cate <- function(plots, stem) {
  ggsave(paste0(stem, "_dist.png"), plots$p_dist, width = 7,  height = 4.5, dpi = 150)
  ggsave(paste0(stem, "_blp.png"),  plots$p_blp,  width = 9,  height = 6,   dpi = 150)
  if (length(plots$mod_plots) > 0)
    ggsave(paste0(stem, "_mods.png"),
           wrap_plots(plots$mod_plots, ncol = 2),
           width = 12, height = 10, dpi = 150)
  cat("Saved CATE plots:", basename(stem), "\n")
}

base_path <- "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/"

# Run CATE for T7 (primary treatment) x all 5 outcomes
cate_targets <- list(
  list(outcome = "O1_rr_index",
       label   = "Racial resentment index [PRIMARY]",
       stem    = "mpsa_cate_rr_nat_retro"),
  list(outcome = "O4_welfare_binary",
       label   = "Welfare spending opposition",
       stem    = "mpsa_cate_welfare_nat_retro"),
  list(outcome = "O5_aidpoor_binary",
       label   = "Aid to poor opposition",
       stem    = "mpsa_cate_aidpoor_nat_retro"),
  list(outcome = "O6_status_jobs_threat",
       label   = "White status / jobs threat",
       stem    = "mpsa_cate_status_threat_nat_retro"),
  list(outcome = "O7_gov_favors_blacks_amt",
       label   = "Government favors Blacks (amount)",
       stem    = "mpsa_cate_gov_favors_blacks_nat_retro")
)

blp_results <- list()

for (ct in cate_targets) {
  cat("\nFitting CF CATE:", ct$label, "x National retrospective (T7)...\n")
  cf_obj <- fit_cf_cate(df1, ct$outcome, "T7_nat_econ_last_worse",
                         controls_cf_raw)
  save_cate(
    plot_cate(cf_obj, ct$label,
              "National retrospective economic pessimism"),
    paste0(base_path, ct$stem)
  )
  blp_results[[ct$outcome]] <- best_linear_projection(cf_obj$forest,
                                                        cf_obj$X)
}

cat("\n", strrep("=", 60), "\n")
cat("BLP SUMMARIES (T7 x all outcomes)\n")
cat(strrep("=", 60), "\n")
for (nm in names(blp_results)) {
  cat("\n--- BLP:", outcome_labels[nm], "x National retrospective ---\n")
  print(blp_results[[nm]])
}

# ============================================================
# balance 
# ============================================================

cat("\n", strrep("=", 60), "\n")
cat("B) COVARIATE BALANCE\n")
cat(strrep("=", 60), "\n\n")

get_balance <- function(df, treat_var, ctrls, seed = 123) {
  d <- df %>% select(all_of(c(treat_var, ctrls))) %>% drop_na()
  cat("Balance:", treatment_labels[treat_var], "| n =", nrow(d), "\n")
  ps_form <- as.formula(paste(treat_var, "~", paste(ctrls, collapse = " + ")))
  tryCatch({
    set.seed(seed)
    m <- matchit(ps_form, data = d, method = "nearest",
                 caliper = 0.2, estimand = "ATT")
    list(matchit = m, treat = treat_var, ok = TRUE)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    list(matchit = NULL, treat = treat_var, ok = FALSE)
  })
}

bal_objects <- setNames(
  map(treats_bin, ~ get_balance(df1, .x, controls_ols_raw)),
  treats_bin
)

iwalk(bal_objects, function(obj, nm) {
  if (!obj$ok) return(invisible(NULL))
  p <- love.plot(
    obj$matchit,
    stats          = "mean.diffs",
    thresholds     = c(m = 0.1),
    abs            = TRUE,
    var.order      = "unadjusted",
    drop.distance  = TRUE,
    colors         = c("Pre-match" = "#e74c3c", "Matched" = "#2980b9"),
    shapes         = c("Pre-match" = 19, "Matched" = 17),
    title          = paste("Covariate Balance:\n", treatment_labels[nm]),
    sample.names   = c("Pre-match", "Matched")
  ) + theme_minimal(base_size = 11)

  fname <- paste0(
    base_path, "mpsa_balance_love_",
    gsub("[^a-zA-Z0-9]", "_", treatment_labels[nm]), ".png"
  )
  ggsave(fname, p, width = 8, height = 6, dpi = 150)
  cat("Saved love plot:", treatment_labels[nm], "\n")
})

extract_smd <- function(bt) {
  df  <- as.data.frame(bt$Balance)
  nms <- names(df)
  cat("  bal.tab columns:", paste(nms, collapse = ", "), "\n")

  un_col  <- intersect(c("Diff.Un",  "M.Threshold.Un",  "SMD.Un"),  nms)[1]
  adj_col <- intersect(c("Diff.Adj", "M.Threshold.Adj", "SMD.Adj"), nms)[1]

  if (is.na(un_col))
    un_col  <- nms[grepl("Un$|Unadj", nms, ignore.case = TRUE)][1]
  if (is.na(adj_col))
    adj_col <- nms[grepl("Adj$", nms, ignore.case = TRUE)][1]

  if (is.na(un_col) || is.na(adj_col)) {
    cat("  WARNING: could not identify SMD columns.\n")
    print(head(df))
    return(NULL)
  }

  cat("  Using:", un_col, "(unadjusted) |", adj_col, "(adjusted)\n")
  df %>%
    rownames_to_column("covariate") %>%
    dplyr::select(covariate,
                  SMD_unadjusted = all_of(un_col),
                  SMD_matched    = all_of(adj_col))
}

balance_table <- map_dfr(bal_objects, function(obj) {
  if (!obj$ok) return(NULL)
  bt  <- bal.tab(obj$matchit, stats = "mean.diffs", un = TRUE)
  smd <- extract_smd(bt)
  if (is.null(smd)) return(NULL)
  smd %>% mutate(treatment        = treatment_labels[obj$treat],
                 balance_achieved = abs(SMD_matched) < 0.1)
}) %>%
  dplyr::select(treatment, covariate,
                SMD_unadjusted, SMD_matched, balance_achieved) %>%
  mutate(across(c(SMD_unadjusted, SMD_matched), ~ round(., 4))) %>%
  arrange(treatment, desc(abs(SMD_unadjusted)))

balance_summary <- balance_table %>%
  group_by(treatment) %>%
  summarise(
    n_covariates      = n(),
    n_balanced        = sum(balance_achieved, na.rm = TRUE),
    pct_balanced      = round(mean(balance_achieved, na.rm = TRUE) * 100, 1),
    max_SMD_prematch  = round(max(abs(SMD_unadjusted), na.rm = TRUE), 4),
    max_SMD_postmatch = round(max(abs(SMD_matched),    na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\n--- Balance Table ---\n");   print(as_tibble(balance_table),   n = Inf)
cat("\n--- Balance Summary ---\n"); print(as_tibble(balance_summary), n = Inf)

write_csv(balance_table,
          paste0(base_path, "mpsa_balance_table.csv"))
write_csv(balance_summary,
          paste0(base_path, "mpsa_balance_summary.csv"))
cat("Saved: mpsa_balance_table.csv, mpsa_balance_summary.csv\n")

cat("\n", strrep("=", 60), "\n")
cat("DONE. All robustness outputs saved.\n")
cat(strrep("=", 60), "\n")
