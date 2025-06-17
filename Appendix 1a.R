# 0. Load required packages
library(readr)      # read_csv
library(dplyr)      # data manipulation
library(tidyr)      # complete, pivot_wider, pivot_longer
library(purrr)      # for list manipulation
library(ggplot2)    # plotting
library(MASS)       # polr
library(stats)      # rank, pnorm

# 1. Import the Data
epds <- read_csv("C://FILE_PATH//FILE_NAME.csv")

# 2. Clean and Prepare the Data
epds <- epds %>%
  rename(post = OUTCOME_VARIABLE_NAME) %>%  # rename variable
  filter(!is.na(post))                      # drop missing

epds <- epds %>%
  rename(group = TREATMENT_VARIABLE_NAME) %>%  # rename variable
  filter(!is.na(post))                      # drop missing
 
# 3. Assign Ranks to the ordinal 'post' column
epds <- epds %>%
  mutate(rpost = rank(post, ties.method = "average"))

# 4. Calculate rank sums and sample sizes by group
rank_sum <- epds %>%
  group_by(group) %>%
  summarise(
    rank_sum = sum(rpost),
    n        = n(),
    .groups  = "drop"
  )

n1 <- rank_sum$n[ rank_sum$group == 1 ]
n2 <- rank_sum$n[ rank_sum$group == 0 ]
N  <- n1 + n2

# 5. Calculate Mann–Whitney U statistic
U1 <- rank_sum$rank_sum[ rank_sum$group == 1 ] - n1*(n1+1)/2
U2 <- rank_sum$rank_sum[ rank_sum$group == 0 ] - n2*(n2+1)/2
U  <- U1 #U1 is the U statistic for our "intervention" arm.

# 6. Tie correction
## 6.1 Frequency of each post by group, then tied pairs per post
freq_pg <- epds %>%
  count(post, group) %>%
  complete(post, group = c(0,1), fill = list(n = 0))

tied_pairs_df <- freq_pg %>%
  group_by(post) %>%
  summarise(
    n1_post    = n[group == 1],
    n0_post    = n[group == 0],
    tied_pairs = n1_post * n0_post,
    .groups    = "drop"
  )

total_tied_pairs <- sum(tied_pairs_df$tied_pairs)

## 6.2 Overall frequency of post, and t_i^3 - t_i
freq_p <- epds %>% count(post)

tie_corr <- freq_p %>%
  mutate(
    t_i3_minus_t_i = if_else(n > 1, n^3 - n, 0)
  )

sum_ti3_minus_t_i <- sum(tie_corr$t_i3_minus_t_i)

# 7. Variance of U with tie correction
Var_U <- (n1 * n2 / 12) * ( (N + 1) - sum_ti3_minus_t_i / (N * (N - 1)) )
SE_U  <- sqrt(Var_U)

# 8. Compute WinP and its 95% CI
WinP    <- U / (n1 * n2)
SE_p    <- SE_U / (n1 * n2)
CI_low  <- WinP - 1.96 * SE_p
CI_high <- WinP + 1.96 * SE_p
CI_low  <- max(0, CI_low)
CI_high <- min(1, CI_high)

# 9. z-statistic and two‐tailed p‐value for H0: WinP = 0.5
z       <- (WinP - 0.5) / SE_p
p_value <- 2 * (1 - pnorm(abs(z)))

# 10. Assemble results
WinP_result <- tibble(
  WinP     = WinP,
  SE_p     = SE_p,
  CI_lower = CI_low,
  CI_upper = CI_high,
  p_value  = p_value
)

print(WinP_result)