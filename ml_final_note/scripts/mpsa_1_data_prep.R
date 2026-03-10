# ============================================================
# mpsa_1_data_prep.R
# Load, clean, construct all treatments/outcomes/controls,
# run MICE, and save environment for downstream scripts.
#
# OUTPUT: mpsa_environment.RData
# ============================================================



suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(mice)
  library(psych)
})

select <- dplyr::select

# ============================================================
# load & clean
# ============================================================

anes_raw <- read_csv(
  "anes_timeseries_2024_csv_20250808/anes_timeseries_2024_csv_20250808.csv",
  show_col_types = FALSE
)

miss_codes <- c(-9, -8, -7, -6, -5, -2, -1, 98, 99, 998, 999)

anes <- anes_raw %>%
  mutate(across(where(is.numeric), ~ ifelse(. %in% miss_codes, NA, .)))

anes_w <- anes %>% filter(V241501x == 1)
cat("N full sample:", nrow(anes), "\n")
cat("N whites-only:", nrow(anes_w), "\n")

# ============================================================
# helpers
# ============================================================

bin_worse_drop3 <- function(x) {
  case_when(
    x %in% c(4, 5) ~ 1L,
    x %in% c(1, 2) ~ 0L,
    TRUE           ~ NA_integer_
  )
}

stars_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < .001 ~ "***", p < .01 ~ "**",
    p < .05  ~ "*",   p < .1  ~ ".", TRUE ~ ""
  )
}

fmt_p <- function(p) format.pval(p, digits = 3, eps = 1e-16)

make_X <- function(df, vars) {
  mm <- model.matrix(
    as.formula(paste0("~ ", paste(vars, collapse = " + "))),
    data = df
  )
  mm[, -1, drop = FALSE]
}

drop_zero_var <- function(X) {
  if (ncol(X) == 0) return(X)
  keep <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 0)
  X[, keep, drop = FALSE]
}

# ============================================================
# treatment set up 
# For binary, neutral (3) droppe -- re-coded where 1=worse, 0=better
# For ontinuous, full 1-5 scale kept
# ============================================================

anes_w <- anes_w %>%
  mutate(
    T7_nat_econ_last_worse  = bin_worse_drop3(V241294x),
    T6_nat_econ_next_worse  = bin_worse_drop3(V241297x),
    T8_pers_fin_last_worse  = bin_worse_drop3(V241451),
    T3_fin_nextyear_worried = bin_worse_drop3(V241452),
    C_nat_econ_last = as.numeric(V241294x),
    C_nat_econ_next = as.numeric(V241297x),
    C_pers_fin_last = as.numeric(V241451),
    C_fin_nextyear  = as.numeric(V241452)
  )

treats_bin <- c(
  "T7_nat_econ_last_worse",
  "T6_nat_econ_next_worse",
  "T8_pers_fin_last_worse",
  "T3_fin_nextyear_worried"
)

treats_cont <- c(
  "C_nat_econ_last",
  "C_nat_econ_next",
  "C_pers_fin_last",
  "C_fin_nextyear"
)

# ============================================================
# outcomes
# ============================================================

anes_w <- anes_w %>%
  mutate(

    # resentment index (higher = more RR)
    RR1 = 6 - V242300,
    RR2 = V242301,
    RR3 = V242302,
    RR4 = 6 - V242303,
    O1_rr_index = rowMeans(across(c(RR1, RR2, RR3, RR4)), na.rm = TRUE),

    # Redistribution
    O4_welfare_binary = case_when(
      V241275x %in% c(1, 2) ~ 0L,
      V241275x %in% c(4, 5) ~ 1L,
      TRUE ~ NA_integer_
    ),
    O5_aidpoor_binary = case_when(
      V241281x %in% c(1, 2) ~ 0L,
      V241281x %in% c(4, 5) ~ 1L,
      TRUE ~ NA_integer_
    ),

    # Racial threat
    O6_status_jobs_threat = case_when(
      is.na(V242519) ~ NA_real_,
      TRUE ~ 6 - as.numeric(V242519)
    ),
    O7_gov_favors_blacks_amt = case_when(
      is.na(V242522x) ~ NA_real_,
      TRUE ~ as.numeric(V242522x) - 4
    ),
  )

outcomes_all <- c(
  "O1_rr_index",            
  "O4_welfare_binary",
  "O5_aidpoor_binary",
  "O6_status_jobs_threat",
  "O7_gov_favors_blacks_amt"
)

# Cronbach's on RR 
rr_cc <- anes_w %>% select(RR1, RR2, RR3, RR4) %>% drop_na()
if (nrow(rr_cc) > 50) {
  cat("\nCronbach's alpha (RR index):\n")
  print(psych::alpha(rr_cc)$total[c("raw_alpha", "std.alpha")])
}

# ============================================================
# labeling for readability 
# ============================================================

treatment_labels <- c(
  "T7_nat_econ_last_worse"  = "National retrospective economic pessimism",
  "T6_nat_econ_next_worse"  = "National prospective economic pessimism",
  "T8_pers_fin_last_worse"  = "Personal financial retrospective pessimism",
  "T3_fin_nextyear_worried" = "Personal financial prospective worry",
  "C_nat_econ_last"         = "National retrospective economic pessimism (continuous)",
  "C_nat_econ_next"         = "National prospective economic pessimism (continuous)",
  "C_pers_fin_last"         = "Personal financial retrospective pessimism (continuous)",
  "C_fin_nextyear"          = "Personal financial prospective worry (continuous)"
)

outcome_labels <- c(
  "O1_rr_index"              = "Racial resentment index",
  "O4_welfare_binary"        = "Welfare spending opposition",
  "O5_aidpoor_binary"        = "Aid to poor opposition",
  "O6_status_jobs_threat"    = "White status / jobs threat",
  "O7_gov_favors_blacks_amt" = "Government favors Blacks (amount)"
)

# ============================================================
# cOnTrOLs
# ============================================================

controls_ols_raw <- c(
  "V241551",   # Gender
  "V241025",   # Party registration
  "V241177",   # Ideology
  "V241420",   # Religion importance
  "V241465x",  # Education
  "V241483",   # Employment status
  "V241497",   # Union membership
  "V241567x"   # Household income
)

controls_cf_extra_raw <- c(
  "V241004", "V241039", "V241104", "V241117",
  paste0("V241", sprintf("%03d", 118:126)),
  "V241232", "V241236", "V241237",
  "V241239", "V241252", "V241255", "V241266x", "V241269x", "V241272x",
  "V241275x", "V241281x", "V241290x", "V241302",
  "V241319x", "V241320", "V241335", "V241372x", "V241375x", "V241378x",
  "V241381x", "V241385x", "V241387", "V241397", "V241422",
  "V241570", "V241571", "V241576", "V241577", "V241578", "V241579",
  "V241607", "V241608", "V241621",
  "V241539", "V241540", "V241554"
)

controls_ols_raw <- intersect(controls_ols_raw, names(anes_w))
controls_cf_raw  <- intersect(
  unique(c(controls_ols_raw, controls_cf_extra_raw)),
  names(anes_w)
)

cat("\nOLS controls:", length(controls_ols_raw), "\n")
cat("CF controls:", length(controls_cf_raw), "\n")

# ============================================================
# ratatouille (m=3, maxit=3)
# ============================================================

mice_vars <- unique(c(
  outcomes_all, treats_bin, treats_cont,
  controls_ols_raw, controls_cf_raw
))
mice_vars <- intersect(mice_vars, names(anes_w))

cat("\nRunning MICE (m=3, maxit=3) on", length(mice_vars), "variables...\n")
imp      <- mice(anes_w[, mice_vars], m = 3, maxit = 3, printFlag = FALSE)
imp_list <- complete(imp, "all")

stopifnot(all(treats_bin  %in% names(imp_list[[1]])))
stopifnot(all(treats_cont %in% names(imp_list[[1]])))
cat("MICE complete. All treatment columns verified.\n")

imp_list <- map(imp_list, function(df) {
  df$media_high <- anes_w$media_high[match(rownames(df), rownames(anes_w))]
  df
})

df1 <- imp_list[[1]]

# ============================================================
# export/save
# ============================================================

save(
  anes_w, imp_list, df1,
  treats_bin, treats_cont, outcomes_all,
  controls_ols_raw, controls_cf_raw,
  media_vars_rw,
  treatment_labels, outcome_labels,
  bin_worse_drop3, stars_p, fmt_p, make_X, drop_zero_var,
  file = "~/Documents/mpsa_dubois/ml_final_note/outputs_ml/mpsa_environment.RData"
)

cat("\nEnvironment saved.\n")
cat("Outcomes (5): O1 [PRIMARY], O4, O5, O6, O7\n")
cat("Media moderator: high = 3+ outlets\n")
