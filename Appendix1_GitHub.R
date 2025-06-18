# -----------------------------------------------------------------------------
# Description: 
#   1. Reads an ordinal outcome dataset.
#   2. Cleans and prepares the data.
#   3. Computes win probability via:
#      a) Mann–Whitney U (rank‐based) method with tie correction.
#      b) Rank‐fraction regression (HC2‐robust).
#   4. Performs:
#      - Exact Wilcoxon test.
#      - Cumulative logistic regression (PO and partial PO) via VGAM.
#      - Proportional‐odds assumption tests (Brant and LRT).
#   5. Generates diagnostic plots.
#
# License: MIT
#
# Prerequisites:
#   install.packages(c(
#     "readr", "dplyr", "tidyr", "purrr", "ggplot2", "MASS", "stats",
#     "sandwich", "lmtest", "coin", "brant", "VGAM"
#   ))
#
# Usage:
#   1. Set `file_path` to your CSV file location.
#   2. Run in R/RStudio: source("Appendix1.R")
# -----------------------------------------------------------------------------

# 0. Clear environment and console
rm(list = ls())
cat("\014")

# 1. Load required packages ---------------------------------------------------
library(readr)    # read_csv()
library(dplyr)    # data manipulation
library(tidyr)    # complete(), pivot_longer()
library(purrr)    # list manipulation (if needed)
library(ggplot2)  # plotting
library(MASS)     # polr()
library(stats)    # pnorm(), plogis()
library(sandwich) # vcovHC()
library(lmtest)   # coeftest()
library(coin)     # wilcox_test()
library(brant)    # brant()
library(VGAM)     # vglm()

# 2. User inputs ---------------------------------------------------------------
file_path <- "C:/FILE_PATH/FILE_NAME.csv"  # <-- Update this path



# 3. Read and clean data ------------------------------------------------------
epds <- read_csv(file_path, show_col_types = FALSE) %>%
  rename(
    post  = OUTCOME_VARIABLE_NAME,  # ordinal outcome (e.g., 0–10)
    group = TREATMENT_VARIABLE_NAME # treatment indicator (0 = control, 1 = intervention)
  ) %>%
  filter(!is.na(post), !is.na(group))

# 4. Compute Win Probability via Mann–Whitney U --------------------------------

# 4.1 Assign average ranks to 'post'
epds <- epds %>%
  mutate(rpost = rank(post, ties.method = "average"))

# 4.2 Summarise rank sums and sample sizes by group
rank_sum <- epds %>%
  group_by(group) %>%
  summarise(
    rank_sum = sum(rpost),
    n        = n(),
    .groups  = "drop"
  )

n1 <- rank_sum$n[rank_sum$group == 1]
n2 <- rank_sum$n[rank_sum$group == 0]
N  <- n1 + n2

# 4.3 Compute U statistic for intervention arm
U1 <- rank_sum$rank_sum[rank_sum$group == 1] - n1 * (n1 + 1) / 2

# 4.4 Tie correction terms
#   - total tied pairs across groups
tied_pairs_df <- epds %>%
  count(post, group) %>%
  complete(post, group = c(0, 1), fill = list(n = 0)) %>%
  group_by(post) %>%
  summarise(
    tied_pairs = n[group == 1] * n[group == 0],
    .groups     = "drop"
  )
total_tied_pairs <- sum(tied_pairs_df$tied_pairs)

#   - sum of (t_i^3 - t_i) over all posts
tie_corr <- epds %>%
  count(post) %>%
  mutate(t_i3_minus_t_i = if_else(n > 1, n^3 - n, 0))
sum_ti3_minus_t_i <- sum(tie_corr$t_i3_minus_t_i)

# 4.5 Variance and SE of U with tie correction
Var_U <- (n1 * n2 / 12) * ((N + 1) - sum_ti3_minus_t_i / (N * (N - 1)))
SE_U  <- sqrt(Var_U)

# 4.6 Win probability estimate and 95% CI
WinP    <- U1 / (n1 * n2)
SE_p    <- SE_U / (n1 * n2)
CI_low  <- max(0, WinP - 1.96 * SE_p)
CI_high <- min(1, WinP + 1.96 * SE_p)

# 4.7 z‐statistic and two‐tailed p‐value for H0: WinP = 0.5
z       <- (WinP - 0.5) / SE_p
p_value <- 2 * (1 - pnorm(abs(z)))

# 4.8 Assemble and print results
WinP_result <- tibble::tibble(
  WinP     = WinP,
  SE_p     = SE_p,
  CI_lower = CI_low,
  CI_upper = CI_high,
  p_value  = p_value
)
print(WinP_result)

# 5. Win Probability via Rank‐Fraction Regression ------------------------------

# 5.1 Compute descending ranks overall and within group
epds <- epds %>%
  mutate(
    rpost_desc   = rank(-post, ties.method = "average"),
    grpost_desc  = ave(-post, group, FUN = function(x) rank(x, ties.method = "average"))
  )

# 5.2 Get group sample sizes (swapped to match 'wins for intervention in our 
#     dataset')
freqcnt <- epds %>%
  count(group) %>%
  mutate(group = 1 - group)

# 5.3 Merge and compute win fractions
WinF <- epds %>%
  dplyr::select(ID, group, post, rpost_desc, grpost_desc) %>%  # force dplyr::select()
  left_join(freqcnt, by = "group") %>%
  mutate(winf = (rpost_desc - grpost_desc) / n) %>%
  dplyr::select(ID, group, post, winf)

# 5.4 Linear regression of win fractions with HC2‐robust SEs
fit_lm <- lm(winf ~ group, data = WinF)
robust <- coeftest(fit_lm, vcov. = vcovHC(fit_lm, type = "HC2"))
est  <- robust["group", "Estimate"]
se   <- robust["group", "Std. Error"]

# 5.5 Derive WinP, log‐odds, CIs, Somers’ D, etc.
WinP_tbl <- tibble::tibble(
  Estimate = est,
  HCStdErr = se,
  WinP     = est / 2 + 0.5
) %>%
  mutate(
    lgt        = log(WinP / (1 - WinP)),
    lgtSe      = HCStdErr / (WinP * (1 - WinP)),
    ln_Odds_l  = lgt - 1.96 * lgtSe,
    ln_Odds_u  = lgt + 1.96 * lgtSe,
    Odds       = exp(lgt),
    Odds_l     = exp(ln_Odds_l),
    Odds_u     = exp(ln_Odds_u),
    WinP_l     = plogis(ln_Odds_l),
    WinP_u     = plogis(ln_Odds_u),
    test       = lgt / lgtSe,
    p_value    = 2 * (1 - pnorm(abs(test))),
    SomersD    = 2 * WinP - 1,
    D_L        = 2 * WinP_l - 1,
    D_U        = 2 * WinP_u - 1,
    win_diff   = 2 * WinP   - 0.5,
    win_diff_l = 2 * WinP_l - 0.5,
    win_diff_u = 2 * WinP_u - 0.5,
    win_ratio   = WinP     / (1 - WinP),
    win_ratio_l = WinP_l   / (1 - WinP_l),
    win_ratio_u = WinP_u   / (1 - WinP_u)
  ) %>%
  dplyr::select(
    WinP, WinP_l, WinP_u,
    win_diff, win_diff_l, win_diff_u,
    win_ratio, win_ratio_l, win_ratio_u,
    test, p_value,
    Odds, Odds_l, Odds_u,
    SomersD, D_L, D_U
  )

print(WinP_tbl)

# 6. Exact Wilcoxon Test -------------------------------------------------------
epds <- epds %>%
  mutate(
    group_fct = factor(group),
    post_fct  = factor(post, ordered = TRUE)
  )
wilcox_exact <- wilcox_test(post ~ group_fct, data = epds, distribution = "exact")
print(wilcox_exact)

# 7. Cumulative Logistic Regression -------------------------------------------

# 7.1 Proportional‐odds model via polr()
fit_polr <- polr(post_fct ~ group_fct, data = epds, method = "logistic", Hess = TRUE)
summary(fit_polr)

# 7.2 PO and partial‐PO via VGAM::vglm()
fit_vglm  <- vglm(post_fct ~ group_fct,
                  family = cumulative(link = "logitlink", parallel = TRUE),
                  data   = epds)
fit_ppom  <- vglm(post_fct ~ group_fct,
                  family = cumulative(link = "logitlink", parallel = FALSE),
                  data   = epds)
summary(fit_vglm)
summary(fit_ppom)

# 8. Proportional‐Odds Assumption Tests ---------------------------------------

# 8.1 Brant test on polr()
brant_res <- brant(fit_polr)
print(brant_res)

# 8.2 Likelihood‐ratio test between PO and partial‐PO
LL1 <- logLik(fit_vglm)
LL2 <- logLik(fit_ppom)
LR  <- as.numeric(2 * (LL2 - LL1))
p   <- 1 - pchisq(LR, df = length(coef(fit_ppom)) - length(coef(fit_vglm)))
cat("LRT p‐value:", p, "\n") # p < 0.05 indicates PO violation

# 9. Diagnostic Plots ---------------------------------------------------------

## Global goodness‐of‐fit (Pearson)
R    <- residuals(fit_vglm, type = "pearson")
X2   <- sum(R^2)
df   <- nrow(R) * ncol(R) - length(coef(fit_vglm))
p_glob <- 1 - pchisq(X2, df)
cat("Pearson X^2 =", X2, " on df =", df, " → p =", p_glob, "\n")

## 9.3 Half‐normal Q–Q plot of |Pearson residuals|
res_long <- as.data.frame(R) %>%
  mutate(ID = seq_len(n())) %>%
  pivot_longer(-ID, names_to = "threshold", values_to = "pearson")

ggplot(res_long, aes(sample = abs(pearson))) +
  stat_qq(distribution = qnorm) +
  stat_qq_line(distribution = qnorm) +
  labs(
    title = "Half‐Normal Q–Q Plot of |Pearson Residuals|",
    x     = "Theoretical Quantiles",
    y     = "Sample Quantiles"
  ) +
  theme_minimal()

