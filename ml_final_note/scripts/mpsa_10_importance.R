# ============================================================
# mpsa_10_importance.R
# Variable importance from continuous causal forests
#
# Uses GRF's built-in variable_importance() — a weighted
# count of how often each variable is used as a splitting
# variable, adjusted for depth. This is the forest's internal
# measure of which covariates matter most for identifying
# treatment effect heterogeneity and adjusting for confounding.
#
# REQUIRES: mpsa_environment.RData (from mpsa_1_data_prep.R)
#
# OUTPUTS (to final_outputs/):
#   fig_10_varimp_primary.png
#   fig_11_varimp_outcome_grid.png
#   fig_12_varimp_welfare_flag.png
#   mpsa_varimp_full.csv  — full importance scores, all combos
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(grf)
  library(ggplot2)
  library(patchwork)
})

select <- dplyr::select

# ============================================================
# all roads lead to documents 
# ============================================================

in_dir  <- "~/Documents/mpsa_dubois/ml_final_note/outputs_ml"
out_dir <- "~/Documents/mpsa_dubois/ml_final_note/final_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(in_dir, "mpsa_environment.RData"))
cat("Environment loaded.\n")
cat("Using df1 (first MICE imputation) for all importance diagnostics.\n\n")

# ============================================================
# config
# ============================================================

N_TOP   <- 15   
SEED    <- 123
N_TREES <- 2000

# welfre in appendix 
outcomes_main <- c(
  "O1_rr_index",
  "O5_aidpoor_binary",
  "O6_status_jobs_threat",
  "O7_gov_favors_blacks_amt"
)

# ============================================================
# readable
# ============================================================

covariate_labels <- c(
  # V241550: PRE: What is R's sex? [revised]
  "V241550"    = "Sex (female)",
  # V241227x: PRE: SUMMARY: Party ID
  "V241227x"   = "Party ID (summary)",
  # V241177: PRE: 7pt scale liberal-conservative self-placement
  "V241177"    = "Ideology (liberal-conservative)",
  # V241420: PRE: Present religion of R
  "V241420"    = "Religion (denomination)",
  # V241465x: PRE: SUMMARY: Respondent 5 Category level of education
  "V241465x"   = "Education",
  # V241483: PRE: Describe R's employment
  "V241483"    = "Employment status",
  # V241497: PRE: Anyone in HH belong to labor union
  "V241497"    = "Union membership (HH)",
  # V241566x: PRE: SUMMARY: Total (household) income
  "V241566x"   = "Household income",
  
  # V241004: PRE: How often does R pay attention to politics and elections
  "V241004"    = "Attention to politics",
  # V241005: PRE: How interested in following campaigns
  "V241005"    = "Interest in campaigns",
  
  # V241025: PRE: Party of registration
  "V241025"    = "Party of registration",
  # V241221: PRE: Party ID: Does R think of self as Dem, Rep, or Ind
  "V241221"    = "Party ID (3-category)",
  # V241222: PRE: Party Identification strong - Democrat Republican
  "V241222"    = "Party ID strength",
  # V241228: PRE: Party identity importance
  "V241228"    = "Party identity importance",
  
  # V241039: PRE: For whom did R vote for President (pre-election)
  "V241039"    = "Pre-election presidential vote",
  # V241104: PRE: Recall of last (2020) Presidential vote choice
  "V241104"    = "2020 presidential vote recall",
  
  # V241137x: PRE: SUMMARY: Approve or disapprove President handling job
  "V241137x"   = "Presidential approval",
  # V241129x: PRE: SUMMARY: Approval of Congress handling its job
  "V241129x"   = "Congress approval",
  # V241133x: PRE: SUMMARY: Approval of Supreme Court handling its job
  "V241133x"   = "Supreme Court approval",
  
  # V241117: PRE: Are things in the country on right track
  "V241117"    = "Right track / wrong track",
  # V241229: PRE: How often trust government in Washington to do what is right
  "V241229"    = "Trust in government",
  # V241231: PRE: Government run by a few big interests or for benefit of all
  "V241231"    = "Govt run by few big interests",
  # V241232: PRE: Does government waste much tax money
  "V241232"    = "Government waste",
  # V241233: PRE: How many in government are corrupt
  "V241233"    = "Government corruption",
  # V241235: PRE: Elections make government pay attention
  "V241235"    = "Elections make govt pay attention",

  # V241118: PRE: How hopeful R feels about how things are going in the country
  "V241118"    = "Emotion: hopeful",
  # V241119: PRE: How afraid R feels about how things are going in the country
  "V241119"    = "Emotion: afraid",
  # V241120: PRE: How outraged R feels about how things are going in the country
  "V241120"    = "Emotion: outraged",
  # V241121: PRE: How angry R feels about how things are going in the country
  "V241121"    = "Emotion: angry",
  # V241122: PRE: How happy R feels about how things are going in the country
  "V241122"    = "Emotion: happy",
  # V241123: PRE: How worried R feels about how things are going in the country
  "V241123"    = "Emotion: worried",
  # V241124: PRE: How proud R feels about how things are going in the country
  "V241124"    = "Emotion: proud",
  # V241125: PRE: How irritated R feels about how things are going in the country
  "V241125"    = "Emotion: irritated",
  # V241126: PRE: How nervous R feels about how things are going in the country
  "V241126"    = "Emotion: nervous",
  
  # V241607: PRE: Agree/disagree - women interpret innocent remarks as sexist
  "V241607"    = "Sexism: women over-interpret remarks",
  # V241608: PRE: Women seek to gain power by getting control over men
  "V241608"    = "Sexism: women seek control over men",
  
  # V241621: PRE: How satisfied is R with life
  "V241621"    = "Life satisfaction",
  # V241570: PRE: Trouble concentrating in past week
  "V241570"    = "Trouble concentrating (past week)",
  # V241571: PRE: Does R have health insurance
  "V241571"    = "Has health insurance",
  # V241576: PRE: Felt little interest or pleasure in doing things (last 2 weeks)
  "V241576"    = "Little interest/pleasure (PHQ-2)",
  # V241577: PRE: Felt down, depressed, or hopeless (last 2 weeks)
  "V241577"    = "Felt depressed/hopeless (PHQ-2)",
  # V241540: PRE: How likely R able to make all housing payments in next 12 months
  "V241540"    = "Able to make housing payments",
  # V241551: PRE: What is R's gender?
  "V241551"    = "Gender identity",
  
  # V241578: PRE: Need to be more sensitive / people too easily offended
  "V241578"    = "Too sensitive / easily offended",
  # V241579: PRE: How often does R self censor
  "V241579"    = "Self-censorship frequency",
  
  # V241248: PRE: 7pt scale abortion: self-placement
  "V241248"    = "Abortion (7pt scale)",
  # V241302: PRE: STD Abortion: self-placement
  "V241302"    = "Abortion stance (STD)",
  # V241252: PRE: 7pt scale guaranteed job-income scale: self-placement
  "V241252"    = "Guaranteed job/income scale",
  # V241255: PRE: 7pt scale gov assistance to blacks scale: self-placement
  "V241255"    = "Govt assistance to Blacks scale",
  # V241239: PRE: 7pt scale spending & services: self-placement
  "V241239"    = "Spending & services scale",
  # V241242: PRE: 7pt scale defense spending: self-placement
  "V241242"    = "Defense spending scale",
  # V241245: PRE: 7pt scale gov-private medical insurance scale: self-placement
  "V241245"    = "Govt-private health insurance scale",
  # V241258: PRE: 7pt scale environment-business tradeoff: self-placement
  "V241258"    = "Environment vs. business tradeoff",
  # V241236: PRE: Which party better: handling economy
  "V241236"    = "Party better: economy",
  # V241237: PRE: Which party better: handling immigration
  "V241237"    = "Party better: immigration",
  # V241387: PRE: Which party better: handling environment
  "V241387"    = "Party better: environment",
  
  # V241263x: PRE: SUMMARY: Federal Budget Spending: Social Security
  "V241263x"   = "Spending: Social Security",
  # V241266x: PRE: SUMMARY: Federal Budget Spending: public schools
  "V241266x"   = "Spending: public schools",
  # V241269x: PRE: SUMMARY: Federal Budget Spending: border security
  "V241269x"   = "Spending: border security",
  # V241272x: PRE: SUMMARY: Federal Budget Spending: dealing with crime
  "V241272x"   = "Spending: crime",
  # V241275x: PRE: SUMMARY: Federal Budget Spending: welfare programs
  "V241275x"   = "Spending: welfare programs",
  # V241278x: PRE: SUMMARY: Federal Budget Spending: highways
  "V241278x"   = "Spending: highways",
  # V241281x: PRE: SUMMARY: Federal Budget Spending: aid to the poor
  "V241281x"   = "Spending: aid to the poor",
  # V241284x: PRE: SUMMARY: Federal Budget Spending: protecting the environment
  "V241284x"   = "Spending: environment",
  
  # V241288: PRE: Approve or disapprove DEI
  "V241288"    = "DEI approval",
  # V241290x: PRE: SUMMARY Approve/disapprove DEI
  "V241290x"   = "DEI approval (summary)",
  # V241317: PRE: Favor or oppose requiring ID when voting
  "V241317"    = "Favor voter ID requirement",
  # V241319x: PRE: SUMMARY: Favor/oppose requiring ID when voting
  "V241319x"   = "Favor voter ID (summary)",
  # V241320: PRE: Favor or oppose allowing felons to vote
  "V241320"    = "Favor felon voting rights",
  # V241376: PRE: Favor/oppose laws protect gays/lesbians against job discrimination
  "V241376"    = "Gay/lesbian job protection",
  # V241378x: PRE: SUMMARY: Favor/oppose laws protect gays/lesbians job discrim
  "V241378x"   = "Gay/lesbian job protection (summary)",
  # V241383: PRE: Favor or oppose right of gay/lesbian couples to legally marry
  "V241383"    = "Gay marriage",
  # V241385x: PRE: SUMMARY: Favor/oppose right of gay/lesbian couples to legally marry
  "V241385x"   = "Gay marriage (summary)",
  # V241379: PRE: Should gay/lesbian couples be allowed to adopt
  "V241381x"   = "Gay/lesbian couples adopt (summary)",
  # V241373: PRE: Favor or oppose banning transgender girls from K-12 girls sports
  "V241373"    = "Ban transgender girls from sports",
  # V241375x: PRE: SUMMARY: Favor/oppose banning transgender girls from sports
  "V241375x"   = "Ban transgender girls sports (summary)",
  # V241370: PRE: Approve/disapprove transgender bathroom policy
  "V241370"    = "Transgender bathroom policy",
  # V241372x: PRE: SUMMARY: Transgender bathroom policy
  "V241372x"   = "Transgender bathroom policy (summary)",
  
  # V241386: PRE: US government policy toward unauthorized immigrants
  "V241386"    = "Policy toward unauthorized immigrants",
  # V241393: PRE: Favor or oppose building a wall on border with Mexico
  "V241393"    = "Favor border wall",
  # V241395x: PRE: SUMMARY: Favor or oppose building a wall
  "V241395x"   = "Favor border wall (summary)",
  # V241410: PRE: Should children brought illegally be sent back or allowed to stay
  "V241410"    = "Children brought illegally: stay/leave",
  # V241407: PRE: Favor or oppose ending birthright citizenship
  "V241407"    = "End birthright citizenship",
  
  # V241335: PRE: How much trust in news media
  "V241335"    = "Trust in news media",
  # V241397: PRE: Best way to deal with urban unrest
  "V241397"    = "Best way to deal with urban unrest",
  
  # V241294x: PRE: SUMMARY: National economy better or worse in last year
  "V241294x"   = "National economy: last year",
  # V241297x: PRE: SUMMARY: Economy better or worse in next 12 months
  "V241297x"   = "National economy: next year",
  # V241451: PRE: R how much better or worse off financially than 1 year ago
  "V241451"    = "Personal finances: last year",
  # V241452: PRE: R how much better or worse off financially next year
  "V241452"    = "Personal finances: next year",
  # V241539: PRE: How worried is R about current financial situation
  "V241539"    = "Financial worry (current)",
  
  # V241600a: PRE: Media used to follow campaign: tv programs
  "V241600a"   = "Media: TV news",
  # V241600b: PRE: Media used to follow campaign: newspapers
  "V241600b"   = "Media: newspapers",
  # V241600c: PRE: Media used to follow campaign: internet sites
  "V241600c"   = "Media: internet news",
  # V241600d: PRE: Media used to follow campaign: radio news
  "V241600d"   = "Media: radio news",
  
  # V241458x: PRE: SUMMARY: Respondent age on Election Day
  "V241458x"   = "Age (Election Day)",
  # V241501x: PRE: SUMMARY: R self-identified race/ethnicity
  "V241501x"   = "Race/ethnicity",
  # V241461x: PRE: SUMMARY: Marital status (5 category)
  "V241461x"   = "Marital status",
  # V241499: PRE: Are you Spanish, Hispanic, or Latino
  "V241499"    = "Hispanic/Latino identity",
  # V241554: PRE: Money invested in Stock Market
  "V241554"    = "Stock market investment",
  # V241567x: PRE: SUMMARY: Total (household) income [six category]
  "V241567x"   = "Household income (6-cat)"
)

label_covariate <- function(varnames) {
  base <- str_remove(varnames, "(?<=[A-Za-z0-9])\\d+$") %>%
    str_remove("TRUE$")
  
  out <- covariate_labels[varnames]
  fallback <- covariate_labels[base]
  ifelse(is.na(out), ifelse(is.na(fallback), varnames, fallback), out)
}

# ============================================================
#  fits
# ============================================================

get_varimp <- function(df, y, t, ctrls,
                       seed = SEED, n_trees = N_TREES) {
  
  treat_lab   <- str_remove(treatment_labels[t], " \\(continuous\\)")
  outcome_lab <- outcome_labels[y]
  
  cat("  Fitting:", treat_lab, "->", outcome_lab, "\n")
  
  d <- df %>% select(all_of(c(y, t, ctrls))) %>% drop_na()
  
  if (nrow(d) < 500) {
    warning(paste("Skipping — n =", nrow(d)))
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
    num.trees     = n_trees,
    min.node.size = 5,
    honesty       = TRUE,
    seed          = seed
  )
  
  # variable_importance() returns a named numeric vector
  # with one score per column of X
  imp_raw <- variable_importance(forest)
  
  tibble(
    outcome         = y,
    treatment       = t,
    outcome_label   = outcome_lab,
    treatment_label = treat_lab,
    variable        = colnames(X),
    importance      = as.numeric(imp_raw)
  ) %>%
    arrange(desc(importance)) %>%
    mutate(
      rank          = row_number(),
      variable_label = label_covariate(variable),
      # Normalize to [0, 1] within each forest for comparability
      importance_norm = importance / max(importance)
    )
}

# ============================================================
# run
# ============================================================

grid <- expand_grid(
  outcome   = outcomes_all,
  treatment = treats_cont
)

cat("Running variable importance on", nrow(grid),
    "treatment-outcome combinations...\n\n")

varimp_list <- pmap(grid, function(outcome, treatment) {
  get_varimp(df1, outcome, treatment, controls_cf_raw)
})

varimp_all <- compact(varimp_list) %>% bind_rows()

# Save full results
write_csv(
  varimp_all,
  file.path(out_dir, "mpsa_varimp_full.csv")
)
cat("\nSaved: mpsa_varimp_full.csv\n\n")

# ============================================================
# themes
# ============================================================

theme_varimp <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, color = "gray40",
                                      hjust = 0),
      plot.caption     = element_text(size = base_size - 2, color = "gray50",
                                      hjust = 0, lineheight = 1.3),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y      = element_text(size = base_size - 1),
      legend.position  = "bottom"
    )
}

treatment_order <- c(
  "National retrospective economic pessimism",
  "Personal financial retrospective pessimism",
  "Personal financial prospective worry",
  "National prospective economic pessimism"
)

# ============================================================
# FIG 10: PRIMARY OUTCOME — ALL TREATMENTS
# Top N_TOP variables for racial resentment x each treatment
# ============================================================

cat("Building Fig 10: Variable importance — primary outcome...\n")

varimp_primary <- varimp_all %>%
  filter(outcome == "O1_rr_index") %>%
  group_by(treatment) %>%
  slice_max(importance, n = N_TOP) %>%
  ungroup() %>%
  mutate(
    treatment_label = factor(treatment_label, levels = treatment_order),
    # Reorder within each facet by importance
    variable_label  = reorder_within(variable_label, importance, treatment_label)
  )

fig10 <- ggplot(varimp_primary,
                aes(x = importance_norm, y = variable_label,
                    fill = treatment_label)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_text(aes(label = round(importance_norm, 2)),
            hjust = -0.1, size = 2.8) +
  scale_y_reordered() +
  scale_x_continuous(limits = c(0, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", ".25", ".50", ".75", "1.0")) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ treatment_label, ncol = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 32)) +
  labs(
    title    = "Variable Importance: Racial Resentment (Primary Outcome)",
    subtitle = paste0("Top ", N_TOP,
                      " covariates by GRF variable importance | Normalized to [0,1] within each forest"),
    x        = "Normalized importance (1.0 = most important covariate in forest)",
    y        = NULL,
    caption  = paste(
      "GRF variable importance = depth-weighted split frequency.",
      "Higher = covariate used more often for splitting in nuisance and CATE models.",
      "Fit on first MICE imputation (df1). Same forest spec as main results.",
      sep = "\n"
    )
  ) +
  theme_varimp()

ggsave(file.path(out_dir, "fig_10_varimp_primary.png"),
       fig10, width = 12, height = 9, dpi = 200)
cat("  Saved: fig_10_varimp_primary.png\n")

# ============================================================
# FIG 11: OUTCOME GRID — MAIN OUTCOMES ONLY
# Top N_TOP for each outcome x treatment (welfare excluded)
# One panel per outcome, facet columns = treatments
# ============================================================

cat("Building Fig 11: Variable importance — outcome grid...\n")

varimp_main <- varimp_all %>%
  filter(outcome %in% outcomes_main) %>%
  group_by(outcome, treatment) %>%
  slice_max(importance, n = N_TOP) %>%
  ungroup() %>%
  mutate(
    treatment_label = factor(treatment_label, levels = treatment_order),
    outcome_label   = factor(outcome_label,
                             levels = outcome_labels[outcomes_main]),
    variable_label  = reorder_within(variable_label, importance,
                                     list(treatment_label, outcome_label))
  )

fig11 <- ggplot(varimp_main,
                aes(x = importance_norm, y = variable_label,
                    fill = outcome_label)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  scale_y_reordered() +
  scale_x_continuous(limits  = c(0, 1.2),
                     breaks  = c(0, 0.5, 1.0),
                     labels  = c("0", ".5", "1")) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(outcome_label ~ treatment_label,
             scales = "free_y",
             labeller = labeller(
               treatment_label = label_wrap_gen(18),
               outcome_label   = label_wrap_gen(20)
             )) +
  labs(
    title    = "Variable Importance: All Main Analysis Outcomes",
    subtitle = paste0("Top ", N_TOP,
                      " covariates | Normalized importance | Welfare spending excluded (see Fig 12)"),
    x        = "Normalized importance",
    y        = NULL,
    caption  = paste(
      "GRF variable importance = depth-weighted split frequency, normalized to [0,1] within each forest.",
      "Welfare spending opposition excluded — see Fig 12 for comparison with racial resentment.",
      sep = "\n"
    )
  ) +
  theme_varimp(base_size = 8) +
  theme(
    strip.text.y  = element_text(angle = 0, size = 7),
    axis.text.y   = element_text(size = 7),
    panel.spacing = unit(0.4, "lines")
  )

ggsave(file.path(out_dir, "fig_11_varimp_outcome_grid.png"),
       fig11, width = 16, height = 14, dpi = 200)
cat("  Saved: fig_11_varimp_outcome_grid.png\n")

# ============================================================
# FIG 12: WELFARE FLAG — RACIAL RESENTMENT VS WELFARE
# ============================================================

cat("Building Fig 12: Variable importance — welfare flag...\n")

nat_retro_label <- "National retrospective economic pessimism"

varimp_flag <- varimp_all %>%
  filter(
    outcome   %in% c("O1_rr_index", "O4_welfare_binary"),
    treatment == "C_nat_econ_last"
  ) %>%
  # Get top N_TOP from union of both outcomes so same vars appear
  group_by(outcome) %>%
  slice_max(importance, n = N_TOP) %>%
  ungroup() %>%
  # Flag theoretically important "problematic" covariates
  mutate(
    flagged = variable_label %in% c(
      "Ideology (liberal-conservative)",
      "Party ID (summary)",
      "Party ID (3-category)",
      "Party ID strength",
      "Party identity importance",
      "Party of registration",
      "Pre-election presidential vote",
      "2020 presidential vote recall",
      "DEI approval",
      "DEI approval (summary)",
      "Spending: welfare programs",
      "Spending: aid to the poor",
      "Govt assistance to Blacks scale"
    ),
    flag_color = case_when(
      flagged & outcome == "O4_welfare_binary" ~ "Party/ideology\n(welfare only)",
      flagged & outcome == "O1_rr_index"       ~ "Party/ideology\n(racial resentment)",
      TRUE                                      ~ "Other covariate"
    ),
    outcome_label = factor(
      outcome_label,
      levels = c("Racial resentment index", "Welfare spending opposition")
    ),
    variable_label = reorder_within(variable_label, importance_norm, outcome_label)
  )

flag_colors <- c(
  "Party/ideology\n(welfare only)"        = "#e74c3c",
  "Party/ideology\n(racial resentment)"   = "#f39c12",
  "Other covariate"                        = "#95a5a6"
)

fig12 <- ggplot(varimp_flag,
                aes(x = importance_norm, y = variable_label,
                    fill = flag_color)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = round(importance_norm, 2)),
            hjust = -0.1, size = 2.8) +
  scale_y_reordered() +
  scale_x_continuous(limits = c(0, 1.2),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c("0", ".25", ".50", ".75", "1.0")) +
  scale_fill_manual(values = flag_colors, name = NULL) +
  facet_wrap(~ outcome_label, ncol = 2, scales = "free_y") +
  labs(
    title    = "Why Welfare Fails: Party ID and Ideology Dominate the Forest",
    subtitle = paste(
      "Treatment: National retrospective economic pessimism",
      paste0("Top ", N_TOP, " covariates | Normalized importance"),
      sep = " | "
    ),
    x        = "Normalized importance",
    y        = NULL,
    caption  = paste(
      "Red = party/ideology-adjacent covariates in the welfare forest.",
      "Orange = same covariates in the racial resentment forest (lower relative importance).",
      "When party ID and ideology dominate the nuisance model for welfare spending,",
      "the forest absorbs treatment variation into those covariates rather than",
      "attributing it to economic pessimism — producing calibration failure and sign reversals.",
      sep = "\n"
    )
  ) +
  theme_varimp() +
  theme(legend.text = element_text(size = 8))

ggsave(file.path(out_dir, "fig_12_varimp_welfare_flag.png"),
       fig12, width = 12, height = 7, dpi = 200)
cat("  Saved: fig_12_varimp_welfare_flag.png\n")

