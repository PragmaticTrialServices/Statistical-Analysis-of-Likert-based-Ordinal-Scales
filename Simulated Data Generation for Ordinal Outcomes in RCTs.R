##############################################################################
# Simulated Data Generation for Ordinal Outcomes in RCTs under the PO model  #
# and its violation, following the ADEMP framework as described in the paper. #
##############################################################################


# Instructions to make changes to code:
# 
# 1. **Change the Ordinal Outcome Distribution:**
#    - Replace the current probability vector in `pylist` with the new one:  
#      `(0.20, 0.10, 0.05, 0.05, 0.05, 0.10, 0.05, 0.05, 0.05, 0.10, 0.20)`. Make sure that the new vector sums to 1 to maintain a valid probability distribution.
# 
# 2. **Change the Mode to Run the PO vs. Non-PO Model:**
#    - Find the flag variable that determines the model type (e.g., ).
#    - To run under the Proportional Odds (PO) assumption, set the flag for `nonPOlist` to `0`. For a Non-Proportional Odds (nonPO) model, set `nonPOlist` to `1`.
#    - If running a nonPO model, check if the attenuation or strengthening of the treatment effect is desired by adjusting the `attenuateeffectlist` flag (typically `0` for strengthening or `1` for attenuating).
#    - Optionally, verify that any model-specific options (such as the `runvglm` flag) are set appropriately based on the sample size and chosen model.
# 
# 3. **Change the Total Sample Size:**
#    - Locate the variable (e.g., `nlist`) that specifies the total sample size for the simulated datasets.
#    - Update the value of `nlist` to one of the desired sizes (300, 1000, 1500, or 2000).
#    - If the sample size is adjusted to a smaller number (like 300), ensure that any conditions or flags (e.g., `runvglm`) that depend on the sample size are updated accordingly to prevent model fitting issues.

# ------------------------------------------------------------------------------
# Function: load_libraries
#
# Description:
#   Loads all required R packages for the simulation study.
#   If a package is not installed, it will install the package with its dependencies.
#   For `winprob` package, please install it using devtools: 
#     install.packages("devtools") 
#     devtools::install_github("proshano/winprob")
#
# Details:
#   The function uses 'requireNamespace' to check package availability and
#   'suppressPackageStartupMessages' to load packages quietly.
#
# Returns:
#   A message confirming that the libraries have been successfully loaded.
# ------------------------------------------------------------------------------
load_libraries <- function() {
  packages <- c("rje", "rms", "MASS", "winprob", "dplyr", "lmtest", "sandwich", "VGAM", "beepr")
  
  # Loop through the list of required packages
  for (pkg in packages) {
    # Install the package if it is not available
    if (!requireNamespace(pkg, quietly = TRUE)) {
      suppressMessages(install.packages(pkg, dependencies = TRUE))
    }
    # Load the package without showing startup messages
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  message("Libraries successfully loaded")
}

# Run the library-loading function
load_libraries()


# ------------------------------------------------------------------------------
# Function: runsimfn
#
# Description:
#   Conducts simulation replications to generate ordinal outcome data based on 
#   a proportional odds (PO) model and its violation.
#
# Parameters:
#   beta_1         - The log odds ratio (treatment effect) for the PO model.
#   py             - A vector of probabilities for each ordinal outcome category 
#                    in the control group (e.g., for the skewed distribution).
#   n              - The total sample size for each simulated dataset.
#   n_sim          - The number of simulation replications.
#   ybinhigh       - The cut-point for dichotomizing the ordinal outcome for binary regression.
#   nonPO          - Flag indicating whether to generate data under non-proportional odds 
#                    (0 = PO holds, 1 = nonPO setting).
#   attenuateeffect- Flag to decide whether to attenuate or strengthen the treatment effect 
#                    in the nonPO setting.
#
# Returns:
#   A numeric vector containing summary performance metrics across the simulations.
#
# Simulation Process:
#   1. Solve for the intercepts (α_k) such that the control group probabilities match the specified values.
#   2. Adjust treatment effects for categories when the PO assumption is violated.
#   3. Compute the “approximate win probability” based on the true OR using Harrell's method.
#   4. Calculate the difference in means for linear regression.
#   5. For each simulation:
#      - Generate treatment assignment (random 0/1 with probability 0.5).
#      - Simulate ordinal outcomes (0 to 10) using cumulative probabilities derived 
#        from the PO model with or without treatment effect modifications.
#      - Fit several models:
#          a) Proportional Odds (PO) ordinal regression using orm()
#          b) Non-parametric win probability (winP) method
#          c) Binary logistic regression on dichotomized outcome
#          d) Linear regression on the ordinal outcome
#          e) Optionally, a nonPO cumulative logit model using vglm() when applicable.
#      - Extract parameter estimates, variance estimates, coverage indicators, and p-value indicators.
#   6. Return summary metrics including bias, coverage probabilities, p-values, and bias test statistics.
#
# ------------------------------------------------------------------------------
runsimfn = function(beta_1, py, n, n_sim, ybinhigh, nonPO, attenuateeffect){
  
  # ------------------------------
  # Solve for intercepts (α_k)
  # ------------------------------
  # yval: ordinal category labels (0, 1, …, 10)
  yval = seq(0, length(py)-1, 1)
  params <- data.frame(yval, py)
  params$rownum = 1:nrow(params)
  params$alpha = 999  # placeholder for intercepts
  
  # Iteratively solve for α_k using the inverse logit (logitlink/expit) functions.
  for(j in 1:nrow(params)){
    if (j == 1) {    # For the first category: P(Y=0)
      params$alpha[j] = logitlink(params$py[j])
    }
    if (j > 1 & j <= 10) {  # For cumulative probabilities P(Y<=k) for k = 1 to 9
      params$alpha[j] = logitlink(params$py[j] + expit(params$alpha[j-1]))
    }
  }
  
  # ------------------------------
  # Set differential treatment effects for nonPO scenarios
  # ------------------------------
  if (nonPO == 0) {
    # Under the PO assumption, additional treatment effect parameters are zero.
    betaa6 = 0
    betaa10 = 0
  }  
  if (nonPO == 1 & attenuateeffect == 1) {  # Attenuated treatment effects for higher outcome categories
    if (beta_1 < 0) {
      betaa6 = 0.05
      betaa10 = 0.10
    }
    if (beta_1 == 0) {
      betaa6 = 0
      betaa10 = 0
    }
    if (beta_1 > 0) {
      betaa6 = -0.05
      betaa10 = -0.10
    }
  }
  if (nonPO == 1 & attenuateeffect == 0) {  # Strengthened treatment effects for higher categories
    if (beta_1 < 0) {
      betaa6 = -0.05
      betaa10 = -0.10
    }
    if (beta_1 == 0) {
      betaa6 = 0
      betaa10 = 0
    }
    if (beta_1 > 0) {
      betaa6 = 0.05
      betaa10 = 0.10
    }
  }
  
  # ------------------------------
  # Calculate the approximate win probability (Harrell's method)
  # ------------------------------
  # See: https://www.fharrell.com/post/wpo/
  capprox_harrell = exp(beta_1)^0.65 / (1 + exp(beta_1)^0.65)
  
  # ------------------------------
  # Calculate the difference in means under the linear model for the ordinal outcome
  # ------------------------------
  mean_tmt0 = (expit(params$alpha[1])) * 0 + 
    (expit(params$alpha[2]) - expit(params$alpha[1])) * 1 + 
    (expit(params$alpha[3]) - expit(params$alpha[2])) * 2 + 
    (expit(params$alpha[4]) - expit(params$alpha[3])) * 3 +
    (expit(params$alpha[5]) - expit(params$alpha[4])) * 4 + 
    (expit(params$alpha[6]) - expit(params$alpha[5])) * 5 +
    (expit(params$alpha[7]) - expit(params$alpha[6])) * 6 + 
    (expit(params$alpha[8]) - expit(params$alpha[7])) * 7 + 
    (expit(params$alpha[9]) - expit(params$alpha[8])) * 8 + 
    (expit(params$alpha[10]) - expit(params$alpha[9])) * 9 +
    (1 - expit(params$alpha[10])) * 10
  
  mean_tmt1 = (expit(params$alpha[1] - beta_1)) * 0 + 
    (expit(params$alpha[2] - beta_1) - expit(params$alpha[1] - beta_1)) * 1 + 
    (expit(params$alpha[3] - beta_1 - betaa6) - expit(params$alpha[2] - beta_1)) * 2 + 
    (expit(params$alpha[4] - beta_1 - betaa6) - expit(params$alpha[3] - beta_1 - betaa6)) * 3 +
    (expit(params$alpha[5] - beta_1 - betaa6) - expit(params$alpha[4] - beta_1 - betaa6)) * 4 + 
    (expit(params$alpha[6] - beta_1 - betaa6) - expit(params$alpha[5] - beta_1 - betaa6)) * 5 +
    (expit(params$alpha[7] - beta_1 - betaa6) - expit(params$alpha[6] - beta_1 - betaa6)) * 6 + 
    (expit(params$alpha[8] - beta_1 - betaa10) - expit(params$alpha[7] - beta_1 - betaa6)) * 7 + 
    (expit(params$alpha[9] - beta_1 - betaa10) - expit(params$alpha[8] - beta_1 - betaa10)) * 8 + 
    (expit(params$alpha[10] - beta_1 - betaa10) - expit(params$alpha[9] - beta_1 - betaa10)) * 9 +
    (1 - expit(params$alpha[10] - beta_1 - betaa10)) * 10
  
  # Set seed for reproducibility of simulations
  set.seed(3)
  
  # Initialize matrices to store simulation results for different methods:
  # Each row corresponds to one simulated dataset.
  POout = matrix(rep(0, n_sim * 4), nrow = n_sim, ncol = 4)          # For PO model results
  PRwinp = POout                                                      # (Not explicitly used in this snippet)
  PostTestwControl = POout                                             # For winP method results
  BinReg = POout                                                       # For binary logistic regression (cut-point at ybinhigh)
  LinReg = POout                                                       # For linear regression on ordinal outcomes
  nonPOout = matrix(rep(0, n_sim * 11), nrow = n_sim, ncol = 11)         # For nonPO model estimates (if applicable)
  
  # ------------------------------
  # Simulation Loop: Repeat for n_sim simulated datasets
  # ------------------------------
  for(j in 1:n_sim){
    
    # Initialize treatment indicator and outcome vector
    tmt = rep(9999, n)
    y = tmt
    df = data.frame(tmt, y)
    
    # Generate data for each individual
    for(i in 1:n){
      # Randomly assign treatment (1) or control (0) with equal probability
      df$tmt[i] = rbinom(1, size = 1, prob = 0.5)
      
      # Generate ordinal outcome (0 to 10) using the specified cumulative probabilities.
      # The probabilities are adjusted based on treatment assignment and, if applicable,
      # differential effects (betaa6, betaa10) for nonPO scenarios.
      df$y[i] = sample(0:10, 1, replace = TRUE, 
                       prob = c(expit(params$alpha[1] - beta_1 * df$tmt[i]), 
                                expit(params$alpha[2] - beta_1 * df$tmt[i]) - expit(params$alpha[1] - beta_1 * df$tmt[i]), 
                                expit(params$alpha[3] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]) - expit(params$alpha[2] - beta_1 * df$tmt[i]),
                                expit(params$alpha[4] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]) - expit(params$alpha[3] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]),
                                expit(params$alpha[5] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]) - expit(params$alpha[4] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]),
                                expit(params$alpha[6] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]) - expit(params$alpha[5] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]),
                                expit(params$alpha[7] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]) - expit(params$alpha[6] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]),
                                expit(params$alpha[8] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i]) - expit(params$alpha[7] - beta_1 * df$tmt[i] - betaa6 * df$tmt[i]),
                                expit(params$alpha[9] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i]) - expit(params$alpha[8] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i]),
                                expit(params$alpha[10] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i]) - expit(params$alpha[9] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i]),
                                1 - expit(params$alpha[10] - beta_1 * df$tmt[i] - betaa10 * df$tmt[i])
                       ))
    }
    
    # ------------------------------
    # Analysis of each simulated dataset
    # ------------------------------
    
    # Check frequency table of y to adjust model fitting for sparse categories
    tablecheck = table(df$y)
    
    # --- Proportional Odds (PO) Model ---
    # Fit an ordinal regression model using the 'orm' function from the rms package.
    fit <- orm(y ~ tmt, data = df)
    # Extract the treatment effect estimate corresponding to the highest observed outcome category
    betaPO = as.numeric(fit$coefficients)[dim(tablecheck)]
    sePO = sqrt(vcov(fit)[2, 2])
    # Determine whether the 95% confidence interval covers the true beta_1
    POoutcov = ifelse(betaPO - 1.96 * sePO <= beta_1 & betaPO + 1.96 * sePO >= beta_1, 1, 0)
    pvaluePOout = 2 * (1 - pnorm(abs(betaPO / sePO)))
    # Save PO model results: estimate, variance, coverage indicator, and significance indicator
    POout[j,] = c(betaPO, sePO^2, POoutcov, ifelse(pvaluePOout < 0.05, 1, 0))
    
    # --- Non-parametric Win Probability (winP) Method ---
    # Calculate win probabilities for the treatment and control groups.
    y = df$y
    tmt = df$tmt
    W1 = (rank(y)[tmt == 1] - rank(y[tmt == 1])) / (n - sum(tmt))
    winphatnonparam = sum(W1) / sum(tmt)
    W1bar = winphatnonparam
    W0 = (rank(y)[tmt == 0] - rank(y[tmt == 0])) / (sum(tmt))
    W0bar = sum(W0) / (n - sum(tmt))
    # Calculate the difference in average win probabilities between groups.
    betahat1 = winphatnonparam - (sum(W0) / (n - sum(tmt)))
    # Transform to the win probability scale for further comparison.
    WinPhat_lm_ = betahat1 / 2 + 0.5
    # Variance estimation for winP (using the method described by Zou, 2022)
    s1sq = sum((W1 - W1bar)^2) / (sum(tmt) - 1)
    s0sq = sum((W0 - W1bar)^2) / ((n - sum(tmt)) - 1)
    varhatwinphat = (s1sq / sum(tmt)) + (s0sq / (n - sum(tmt)))
    TT_norm = (WinPhat_lm_ - 0.5) / sqrt(varhatwinphat)
    p_TT_norm = 2 * (1 - pnorm(abs(TT_norm)))
    p_TT_normsig = ifelse(p_TT_norm < 0.05, 1, 0)
    cl_l_norm = WinPhat_lm_ - 1.96 * sqrt(varhatwinphat)
    cl_u_norm = WinPhat_lm_ + 1.96 * sqrt(varhatwinphat)
    cl_norm_cov = ifelse(cl_l_norm <= capprox_harrell & cl_u_norm >= capprox_harrell, 1, 0)
    # Save winP method results: estimate, variance, coverage indicator, and significance indicator
    PostTestwControl[j,] = c(WinPhat_lm_, varhatwinphat, cl_norm_cov, p_TT_normsig)
    
    # --- Binary Logistic Regression ---
    # Dichotomize the ordinal outcome at the specified cut-point ybinhigh.
    df$ybin = ifelse(df$y >= ybinhigh, 1, 0)
    fitbin = glm(ybin ~ tmt, family = "binomial", data = df)
    betabin = as.numeric(fitbin$coefficients[2])
    varbin = vcov(fitbin)[2, 2]
    # Check coverage against an adjusted treatment effect (beta_1 + betaa6)
    covbin = ifelse(betabin - 1.96 * sqrt(varbin) <= beta_1 + betaa6 & betabin + 1.96 * sqrt(varbin) >= beta_1 + betaa6, 1, 0)
    pbin = 2 * (1 - pnorm(abs(betabin / sqrt(varbin))))
    # Save binary regression results: estimate, variance, coverage, and significance indicator
    BinReg[j,] = c(betabin, varbin, covbin, ifelse(pbin < 0.05, 1, 0))
    
    # --- Linear Regression ---
    fitlin = lm(y ~ tmt, data = df)
    betalin = as.numeric(fitlin$coefficients[2])
    varlin = vcov(fitlin)[2, 2]
    # Coverage check based on the difference in means between treatment groups
    covlin = ifelse(betalin - 1.96 * sqrt(varlin) <= (mean_tmt1 - mean_tmt0) & betalin + 1.96 * sqrt(varlin) >= (mean_tmt1 - mean_tmt0), 1, 0)
    plin = 2 * (1 - pnorm(abs(betalin / sqrt(varlin))))
    # Save linear regression results: estimate, variance, coverage, and significance indicator
    LinReg[j,] = c(betalin, varlin, covlin, ifelse(plin < 0.05, 1, 0))
    
    # --- Non-Proportional Odds (nonPO) Model (Optional) ---
    # If runvglm is enabled, fit a cumulative logit model that allows non-constant treatment effects
    if (runvglm == 1) {
      nonpovglm <- vglm(y ~ tmt, family = cumulative(parallel = FALSE), data = df)
      # Check if the number of coefficients is as expected; extract coefficients for categories 2-11
      if (length(coef(nonpovglm)) == ((length(py) - 1) * 2)) {
        nonPOout[j, 1:10] = as.numeric(coef(nonpovglm)[11:20]) * (-1)
      }
      if (length(coef(nonpovglm)) < ((length(py) - 1) * 2)) {
        nonPOout[j, 11] = 1
      }
    } # end if runvglm
    
  } # end simulation loop
  
  # ------------------------------
  # Return Summary Metrics
  # ------------------------------
  # The function returns a vector containing:
  #   - Bias of the PO model estimate (mean difference from true beta_1)
  #   - Coverage and significance indicators for the PO model, winP method, binary regression, and linear regression
  #   - Bias test statistics (Z-tests) for each method and selected nonPO estimates
  return(c(
    mean(POout[,1]) - beta_1, mean(POout[,3]), mean(POout[,4]),
    mean(PostTestwControl[,1]) - capprox_harrell, mean(PostTestwControl[,3]), mean(PostTestwControl[,4]), 
    mean(BinReg[,1]) - (beta_1 + betaa6), mean(BinReg[,3]), mean(BinReg[,4]),
    mean(LinReg[,1]) - (mean_tmt1 - mean_tmt0), mean(LinReg[,3]), mean(LinReg[,4]),
    # Bias Z-test statistics: if the absolute value > 1.96, it may indicate a problem.
    (mean(POout[,1]) - beta_1) / (sd(POout[,1]) / sqrt(n_sim)),
    (mean(PostTestwControl[,1]) - capprox_harrell) / (sd(PostTestwControl[,1]) / sqrt(n_sim)),
    (mean(BinReg[,1]) - (beta_1 + betaa6)) / (sd(BinReg[,1]) / sqrt(n_sim)),
    (mean(LinReg[,1]) - (mean_tmt1 - mean_tmt0)) / (sd(LinReg[,1]) / sqrt(n_sim)),
    # Bias Z-tests for additional nonPO estimates across categories
    (mean(nonPOout[,1]) - beta_1) / (sd(nonPOout[,1]) / sqrt(n_sim)),  # for category 1
    (mean(nonPOout[,2]) - beta_1) / (sd(nonPOout[,2]) / sqrt(n_sim)),  # for category 2
    (mean(nonPOout[,3]) - (beta_1 + betaa6)) / (sd(nonPOout[,3]) / sqrt(n_sim)),  
    (mean(nonPOout[,4]) - (beta_1 + betaa6)) / (sd(nonPOout[,4]) / sqrt(n_sim)),  
    (mean(nonPOout[,5]) - (beta_1 + betaa6)) / (sd(nonPOout[,5]) / sqrt(n_sim)),  
    (mean(nonPOout[,6]) - (beta_1 + betaa6)) / (sd(nonPOout[,6]) / sqrt(n_sim)),  
    (mean(nonPOout[,7]) - (beta_1 + betaa6)) / (sd(nonPOout[,7]) / sqrt(n_sim)),  # for category 7
    (mean(nonPOout[,8]) - (beta_1 + betaa10)) / (sd(nonPOout[,8]) / sqrt(n_sim)),  
    (mean(nonPOout[,9]) - (beta_1 + betaa10)) / (sd(nonPOout[,9]) / sqrt(n_sim)),  
    (mean(nonPOout[,10]) - (beta_1 + betaa10)) / (sd(nonPOout[,10]) / sqrt(n_sim)), 
    sum(nonPOout[,11])
  ))
  
} # end of runsimfn function


# ------------------------------------------------------------------------------
# Simulation Specifications and Execution
#
# The following section sets up the simulation parameters based on the study design,
# generates simulated datasets for varying treatment effects, and computes performance
# measures for each method (PO model, winP, binary logistic regression, and linear regression).
# ------------------------------------------------------------------------------
# Specify treatment effect scenarios (odds ratios converted to log scale)
beta_1list = c(log(0.8), log(0.9), log(1.0)) #, log(1.1), log(1.2))

# Specify whether to use non-proportional odds and whether to attenuate the effect
nonPOlist = 0
attenuateeffectlist = 0

# Specify the probabilities for each ordinal outcome category (skewed distribution)
pylist = c(0.6, 0.1, 0.04, 0.04, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02, 0.06)

# Set the sample size and number of simulation replications
nlist = 1500     # Total sample size for each dataset
n_simlist = 1000  # Number of simulated datasets

# Specify the threshold for the binary outcome (for binary logistic regression)
ybinhighlist = 7

# Control flag for running the nonPO model (avoid errors in small datasets)
runvglm = 0
if (nlist < 1000){
  runvglm = 0  # Do not run the nonPO logit regression if the sample size is too small
}

# Initialize a matrix to store simulation summary results
tableout = matrix(rep(0, length(beta_1list) * 29), nrow = length(beta_1list), ncol = 29)

# Loop through each treatment effect scenario
for (p in 1:length(beta_1list)){
  
  print(Sys.time())
  tableout[p, 1] = nlist
  tableout[p, 2] = exp(beta_1list[p])  # Convert log-odds to odds ratio for reporting
  # Run simulation for current scenario and store summary metrics in tableout
  tableout[p, 3:29] = runsimfn(beta_1list[p], pylist, nlist, n_simlist, ybinhighlist, nonPOlist, attenuateeffectlist)
  print(p)
}

# Convert simulation results to a data frame for reporting
dfout = as.data.frame(tableout)
colnames(dfout) = c('n', 'beta1', 'poregbias', 'poregcov', 'poregp', 
                    'wpbias', 'wpcov', 'wpp', 
                    'bin7bias', 'bin7cov', 'bin7p', 
                    'linbias', 'lincov', 'linp',
                    'pobiasz', 'wpbiasz', 'bin7biasz', 'linbiasz',
                    'bin1biasz', 'bin2biasz', 'bin3biasz', 'bin4biasz', 
                    'bin5biasz', 'bin6biasz', 'bin7biasz', 'bin8biasz',
                    'bin9biasz', 'bin10biasz', 'sumnonpo')

# Output the results
dfout

# Play a sound notification upon completion (requires the 'beepr' package)
beep("fanfare")

# Copy the rounded simulation summary table to clipboard (for pasting into reports)
write.table(round(tableout, 3), "clipboard", sep = "\t")

# Extract and display selected performance measures for PO, winP, binary, and linear regressions
dfoutp = dplyr::select(dfout, c("n","beta1", "poregp", "wpp", "bin7bias", "bin7cov", "bin7p", "linbias", "lincov", "linp"))
write.table(round(dfoutp, 3), "clipboard", sep = "\t")
write.csv(dfoutp, "dfoutp - 1500.csv", row.names = FALSE)
cat("Simulation complete.\nCSV File with dfoutp was saved in working directory.")




