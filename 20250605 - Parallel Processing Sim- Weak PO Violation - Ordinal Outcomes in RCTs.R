# Clear environment & console
rm(list = ls());  cat("\014")

setwd("C://Users//aalja//OneDrive//01_Work//01_Projects//03_Analysis of ordinal data//01_A guide for ordinal data analysis//Code//R Code")

################################################################################
# SIMULATION OF ORDINAL OUTCOMES IN RCTS (GENERALIZED & PARALLELIZED)
#
# This updated script implements a two-step speed-up strategy:
# 1. Vectorized Data Generation: Efficiently creates simulated datasets.
# 2. Parallelized Outer Loop: Uses the 'future' framework for modern,
#    cross-platform parallel processing with reproducible seeds.
#
# Key Improvements:
# - Generalization: Works with any number of ordinal categories (pylist length).
# - Modern Parallelism: Uses `future` and `future.apply` for easy setup,
#   automatic variable export, and clean shutdowns.
# - Progress Reporting: Integrates with `progressr` for clear updates.
# - Robustness: Code is structured with helper functions to minimize
#   duplication and improve clarity.
# - Reproducibility: A single seed (14424) is used for all parallel tasks.
################################################################################
# Clear environment & console
rm(list = ls());  cat("\014")

# Clear environment & console
rm(list = ls());  cat("\014")

################################################################################
# SIMULATION OF ORDINAL OUTCOMES IN RCTS (GENERALIZED & PARALLELIZED)
#
# This updated script implements a two-step speed-up strategy:
# 1. Vectorized Data Generation: Efficiently creates simulated datasets.
# 2. Parallelized Outer Loop: Uses the 'future' framework for modern,
#    cross-platform parallel processing with reproducible seeds.
#
# Key Improvements:
# - Generalization: Works with any number of ordinal categories (pylist length).
# - Modern Parallelism: Uses `future` and `future.apply` for easy setup,
#   automatic variable export, and clean shutdowns.
# - Progress Reporting: Integrates with `progressr` for clear updates.
# - Robustness: Code is structured with helper functions to minimize
#   duplication and improve clarity.
# - Reproducibility: A single seed (14424) is used for all parallel tasks.
################################################################################

# ------------------------------------------------------------------------------
# Function: load_libraries
# ------------------------------------------------------------------------------
load_libraries <- function() {
  # Added future, future.apply, progressr for modern parallel processing
  packages <- c("rms", "MASS", "winprob", "dplyr", "VGAM", "beepr",
                "future", "future.apply", "progressr")
  
  cat("Loading required R packages...\n")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "winprob") {
        message(paste("Package 'winprob' is not installed.",
                      "Please install it from GitHub using:",
                      "install.packages('devtools')",
                      "devtools::install_github('proshano/winprob')", sep="\n"))
        stop("Required package 'winprob' not found.")
      }
      suppressMessages(install.packages(pkg, dependencies = TRUE))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  message("Libraries successfully loaded.")
}

# Run the library-loading function
load_libraries()

# ==============================================================================
# 1. CORE SIMULATION AND WRAPPER FUNCTIONS
# ==============================================================================

# helper to save with fallback
save_csv_with_fallback <- function(df, proposed_name) {
  tryCatch(
    {
      write.csv(df, proposed_name, row.names = FALSE)
      message("✔ Saved to ", proposed_name)
    },
    error = function(e) {
      # generate 10 random digits
      rand_tag <- paste0(sample(0:9, 10, replace = TRUE), collapse = "")
      fallback <- paste0("file_name_random_", rand_tag, ".csv")
      write.csv(df, fallback, row.names = FALSE)
      warning(
        "Failed to save to '", proposed_name, 
        "'; saved to fallback '", fallback, "' instead."
      )
    }
  )
}

# ------------------------------------------------------------------------------
# Helper Function: calculate_intercepts
# Description: Solves for the intercepts (α_k) for any number of categories.
# ------------------------------------------------------------------------------
calculate_intercepts <- function(py) {
  K <- length(py)
  if (K <= 1) return(numeric(0))
  
  alpha_vec <- numeric(K - 1)
  
  # Calculate first intercept for P(Y=0)
  alpha_vec[1] <- stats::qlogis(py[1])
  
  # Iteratively solve for remaining intercepts
  if (K > 2) {
    for (k in 2:(K - 1)) {
      cumulative_prob_prev <- stats::plogis(alpha_vec[k - 1])
      alpha_vec[k] <- stats::qlogis(py[k] + cumulative_prob_prev)
    }
  }
  return(alpha_vec)
}

# ------------------------------------------------------------------------------
# Helper Function: generate_npo_effects
# Description: Generates the non-proportional odds effect vector.
# Note: This function now correctly implements the PO violation specified in the prompt.
# ------------------------------------------------------------------------------
generate_npo_effects <- function(K, beta_1, nonPO, attenuateeffect) {
  if (!nonPO) return(rep(0, K - 1))
  
  npo_effects <- numeric(K - 1)
  
  if (K == 11) { # Logic for 11 categories (outcome 0-10)
    # Per prompt: PO violation is δ1=0; δ2–6=0.05; δ7–10=0.10
    effect_val1 <- 0.05
    effect_val2 <- 0.10
    
    if (attenuateeffect == 1) { # Attenuate
      if (beta_1 > 0) { effect_val1 <- -effect_val1; effect_val2 <- -effect_val2 }
    } else { # Strengthen
      if (beta_1 < 0) { effect_val1 <- -effect_val1; effect_val2 <- -effect_val2 }
    }
    
    # npo_effects vector is 1-based, corresponding to thresholds for P(Y<=0), P(Y<=1), ...
    # Prompt: δ_1=0 is the default, as the vector is initialized with zeros.
    # Prompt: δ_2 to δ_6 = 0.05 => Corresponds to R indices 2 through 6
    npo_effects[2:6] <- effect_val1
    # Prompt: δ_7 to δ_10 = 0.10 => Corresponds to R indices 7 through 10
    npo_effects[7:10] <- effect_val2
    
  } else {
    message("NOTE: Default non-PO effect logic is specific to 11 categories. Returning zeros for K != 11.")
  }
  
  return(npo_effects)
}


# ------------------------------------------------------------------------------
# Helper Function: get_pmf_from_lp
# Description: Calculates the probability mass function (PMF) from linear predictors
#              of cumulative probabilities.
# ------------------------------------------------------------------------------
get_pmf_from_lp <- function(lp_vec) {
  K_minus_1 <- length(lp_vec)
  K <- K_minus_1 + 1
  
  pmf <- numeric(K)
  if (K_minus_1 == 0) {
    pmf[1] <- 1
    return(pmf)
  }
  
  cum_probs <- stats::plogis(lp_vec)
  
  pmf[1] <- cum_probs[1]
  if (K > 2) {
    for (k in 2:K_minus_1) {
      pmf[k] <- cum_probs[k] - cum_probs[k - 1]
    }
  }
  pmf[K] <- 1 - cum_probs[K_minus_1]
  
  return(pmf)
}


# ────────────────────────────────────────────────────────────────
# 1.1.  Core function that runs ONE replicate
# ────────────────────────────────────────────────────────────────
run_one_simulation <- function(sim_index,
                               beta_1, K, n, ybinhigh,
                               alpha_vec, npo_effects_vec,
                               capprox_harrell, mean_diff_true,
                               runvglm, true_beta_bin,
                               true_grp_coeffs) { # Pass in true aggregate values
  
  tmt_vec <- rbinom(n, 1, 0.5)
  eff_beta_1_vec      <- beta_1 * tmt_vec
  eff_npo_effects_mat <- t(npo_effects_vec %o% tmt_vec)
  
  lp_mat   <- matrix(alpha_vec, nrow = n, ncol = K - 1, byrow = TRUE) -
    eff_beta_1_vec - eff_npo_effects_mat
  P_cum_mat <- plogis(lp_mat)
  
  y_vec <- rowSums(runif(n) > P_cum_mat)
  df    <- data.frame(tmt = tmt_vec, y = y_vec)
  
  PO_res <- WinP_res <- BinReg_res <- LinReg_res <- rep(NA_real_, 4)
  
  nonPO_grp_coeffs <- rep(NA_real_, 3)
  nonPO_grp_covs   <- rep(NA_real_, 3)
  nonPO_grp_pows   <- rep(NA_real_, 3)
  any_sig_flag     <- NA_real_ # NEW: Flag for if any group is significant
  nonPO_flags      <- c(1, NA, 0)   # (errFlag, LRT p<0.05, converged)
  
  n0 <- sum(df$tmt == 0); n1 <- n - n0
  if (n0 > 1 && n1 > 1 && var(df$y, na.rm=TRUE) > 0) {
    
    # ... (PO, WinP, Linear, Binary models remain unchanged) ...
    fitPO <- try(rms::orm(y ~ tmt, data = df), silent = TRUE)
    if (!inherits(fitPO, "try-error") && "tmt" %in% names(coef(fitPO))) {
      betaPO <- coef(fitPO)["tmt"]; sePO <- try(sqrt(vcov(fitPO)["tmt","tmt"]), silent=TRUE)
      if (!inherits(sePO, "try-error") && is.finite(sePO) && sePO > 0) {
        PO_res <- c(betaPO, sePO^2, ifelse(betaPO - 1.96*sePO <= beta_1 & betaPO + 1.96*sePO >= beta_1, 1, 0), ifelse(2*pnorm(abs(betaPO/sePO), lower.tail = FALSE) < .05, 1, 0))
      }
    }
    W1 <- (rank(df$y)[df$tmt == 1] - rank(df$y[df$tmt == 1])) / n0; W0 <- (rank(df$y)[df$tmt == 0] - rank(df$y[df$tmt == 0])) / n1
    WinPhat <- (mean(W1) - mean(W0))/2 + .5; varhatWP <- var(W1)/n1 + var(W0)/n0
    if (is.finite(varhatWP) && varhatWP > 0) {
      WinP_res <- c(WinPhat, varhatWP, ifelse(WinPhat - 1.96*sqrt(varhatWP) <= capprox_harrell & WinPhat + 1.96*sqrt(varhatWP) >= capprox_harrell, 1, 0), ifelse(2*pnorm(abs((WinPhat - 0.5)/sqrt(varhatWP)), lower.tail = FALSE) < .05, 1, 0))
    }
    fitLin <- try(lm(y ~ tmt, data = df), silent=TRUE)
    if (!inherits(fitLin, "try-error")) {
      betLin <- coef(fitLin)["tmt"]; varLin <- vcov(fitLin)[2,2]
      if(is.finite(betLin) && is.finite(varLin) && varLin > 0) {
        LinReg_res <- c(betLin, varLin, ifelse(betLin - 1.96*sqrt(varLin) <= mean_diff_true & betLin + 1.96*sqrt(varLin) >= mean_diff_true, 1, 0), ifelse(2*pnorm(abs(betLin/sqrt(varLin)), lower.tail = FALSE) < .05, 1, 0))
      }
    }
    df$ybin <- as.integer(df$y >= ybinhigh)
    if (length(unique(df$ybin))==2) {
      fitBin <- try(glm(ybin ~ tmt, binomial, data = df), silent=TRUE)
      if (!inherits(fitBin, "try-error")) {
        betBin <- coef(fitBin)["tmt"]; varBin <- vcov(fitBin)["tmt","tmt"]
        if(is.finite(betBin) && is.finite(varBin) && varBin > 0) {
          BinReg_res <- c(betBin, varBin, ifelse(betBin - 1.96*sqrt(varBin) <= true_beta_bin & betBin + 1.96*sqrt(varBin) >= true_beta_bin, 1, 0), ifelse(2*pnorm(abs(betBin/sqrt(varBin)), lower.tail = FALSE) < .05, 1, 0))
        }
      }
    }
    
    ## ---- VGAM non-PO (optional) ----
    if (runvglm == 1 && length(unique(df$y)) >= 3) {
      df$y <- factor(df$y, levels = 0:(K-1), ordered = TRUE)
      suppressWarnings({
        fitNPO  <- try(VGAM::vglm(y ~ tmt, family = VGAM::cumulative(parallel = FALSE, reverse = TRUE), data = df, trace = FALSE), silent = TRUE)
        fitNull <- try(VGAM::vglm(y ~ 1, family = VGAM::cumulative(parallel = TRUE, reverse = TRUE), data = df, trace = FALSE), silent = TRUE)
      })
      
      if (!inherits(fitNPO, "try-error") && !inherits(fitNull, "try-error")) {
        nonPO_flags[3] <- 1 # converged
        nonPO_flags[1] <- 0 # no error
        co <- VGAM::coef(fitNPO)
        tnames <- grep("^tmt", names(co), value = TRUE)
        vcv <- try(VGAM::vcov(fitNPO), silent = TRUE)
        
        if (length(tnames) == (K-1) && !inherits(vcv, "try-error")) {
          grp_indices <- list(g1 = 1, g2 = 2:6, g3 = 7:10)
          
          for (i in 1:3) {
            grp_idx <- grp_indices[[i]]; k <- length(grp_idx)
            avg_coeff <- mean(co[tnames[grp_idx]])
            nonPO_grp_coeffs[i] <- avg_coeff
            
            vcv_sub <- vcv[tnames[grp_idx], tnames[grp_idx], drop=FALSE]
            var_avg_coeff <- sum(vcv_sub) / (k^2)
            
            if (is.finite(var_avg_coeff) && var_avg_coeff > 0) {
              se_avg_coeff <- sqrt(var_avg_coeff)
              lower_ci <- avg_coeff - 1.96 * se_avg_coeff; upper_ci <- avg_coeff + 1.96 * se_avg_coeff
              nonPO_grp_covs[i] <- ifelse(lower_ci <= true_grp_coeffs[i] & upper_ci >= true_grp_coeffs[i], 1, 0)
              p_value <- 2 * pnorm(abs(avg_coeff / se_avg_coeff), lower.tail = FALSE)
              nonPO_grp_pows[i] <- ifelse(p_value < 0.05, 1, 0)
            }
          }
          # NEW: Check if any of the group power flags is 1
          any_sig_flag <- as.numeric(any(nonPO_grp_pows == 1, na.rm = TRUE))
        }
        
        ll_full <- suppressWarnings(as.numeric(VGAM::logLik(fitNPO))); ll_null <- suppressWarnings(as.numeric(VGAM::logLik(fitNull)))
        if (is.finite(ll_full) && is.finite(ll_null)) {
          df_diff <- length(co) - length(VGAM::coef(fitNull))
          if (df_diff > 0) {
            lrt_stat <- 2 * (ll_full - ll_null)
            if (is.finite(lrt_stat) && lrt_stat < 0) lrt_stat <- 0
            if (is.finite(lrt_stat)) {
              pval <- pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
              if (is.finite(pval)) nonPO_flags[2] <- as.numeric(pval < 0.05)
            }
          }
        }
      }
    }
  }
  
  # UPDATED: Return vector now contains the new "any significant" flag
  c(PO_res, WinP_res, BinReg_res, LinReg_res,
    nonPO_grp_coeffs, nonPO_grp_covs, nonPO_grp_pows,
    any_sig_flag, nonPO_flags)
}


# ──────────────────────────────────────────────────────────────
# 1.2.  Wrapper for ONE scenario – returns full + power rows
# ──────────────────────────────────────────────────────────────
run_one_scenario <- function(beta_1, nonPO, atten, sym_flag, n, n_simlist = 1000,
                             ybinhigh  = 7, runvglm = 1, seed = 14424) {
  
  pylist <- if (sym_flag == 1) {
    c(0.20,0.10,0.05,0.05,0.05,0.10,0.05,0.05,0.05,0.10,0.20)
  } else {
    c(0.60,0.10,0.04,0.04,0.04,0.04,0.02,0.02,0.02,0.02,0.06)
  }
  if (abs(sum(pylist)-1) > 1e-6) stop("pylist must sum to 1")
  K <- length(pylist)
  
  ident <- data.frame(beta1_log = beta_1, nonPO = nonPO, atten = atten,
                      symmetric = sym_flag, n = n,
                      sym_tag   = ifelse(sym_flag == 1, "Symmetric","Heavy0"),
                      po_tag    = ifelse(nonPO == 0, "POMet","POViolated"),
                      eff_tag   = ifelse(atten == 1, "attenuatedeffects", "strengthenedeffects"),
                      stringsAsFactors = FALSE)
  
  alpha_vec <- calculate_intercepts(pylist)
  npo_vec   <- generate_npo_effects(K, beta_1, nonPO, atten)
  capprox_harrell <- exp(beta_1)^0.65 / (1+exp(beta_1)^0.65)
  pmf_t0    <- get_pmf_from_lp(alpha_vec)
  pmf_t1    <- get_pmf_from_lp(alpha_vec - beta_1 - npo_vec)
  mean_diff_true <- sum((0:(K-1)) * (pmf_t1 - pmf_t0))
  true_beta_bin  <- beta_1 + npo_vec[ybinhigh]
  
  true_npo <- beta_1 + npo_vec
  grp_indices <- list(g1 = 1, g2 = 2:6, g3 = 7:10)
  true_grp_coeffs <- c(mean(true_npo[grp_indices$g1]),
                       mean(true_npo[grp_indices$g2]),
                       mean(true_npo[grp_indices$g3]))
  
  future.apply::future_lapply(
    1:n_simlist,
    future.seed = seed,
    FUN = run_one_simulation,
    beta_1 = beta_1, K = K, n = n, ybinhigh = ybinhigh,
    alpha_vec = alpha_vec, npo_effects_vec = npo_vec,
    capprox_harrell = capprox_harrell, mean_diff_true = mean_diff_true,
    runvglm = runvglm, true_beta_bin = true_beta_bin,
    true_grp_coeffs = true_grp_coeffs
  ) -> res_list
  
  raw <- do.call(rbind, res_list)
  
  #Standardised bias
  # z_fun <- function(est, truth) {
  #   valid_est <- est[!is.na(est)]; if (length(valid_est) < 2) return(NA)
  #   (mean(valid_est) - truth) / (sd(valid_est) / sqrt(length(valid_est)))
  # }
  
  # Relative bias  
  z_fun <- function(est, truth) {
    valid_est <- est[!is.na(est)]
    if (length(valid_est) < 2) return(NA)
    
    if (truth == 0) { # Use absolute bias when logOR = 0
      return(mean(valid_est) - truth)
    } else { # Use relative bias when logOR != 0
      return((mean(valid_est) - truth) / truth)
    }
  }
  
  # UPDATED: Indexing for the new "any significant" flag
  idx_po            <- 1:4; idx_wp <- 5:8; idx_bin <- 9:12; idx_lin <- 13:16
  idx_npo_est       <- 17:19 # 3 group coefficients
  idx_npo_cov       <- 20:22 # 3 group coverage flags
  idx_npo_pow       <- 23:25 # 3 group power flags
  idx_npo_any_pow   <- 26    # NEW: The "any significant" flag
  idx_flag          <- 27:29 # Shifted
  
  po  <- c(mean(raw[,1], na.rm=T) - beta_1, colMeans(raw[,c(3,4)], na.rm=T))
  wp  <- c(mean(raw[,5], na.rm=T) - capprox_harrell, colMeans(raw[,c(7,8)], na.rm=T))
  bin <- c(mean(raw[,9], na.rm=T) - true_beta_bin, colMeans(raw[,c(11,12)], na.rm=T))
  lin <- c(mean(raw[,13], na.rm=T)- mean_diff_true, colMeans(raw[,c(15,16)], na.rm=T))
  zb  <- c(z_fun(raw[,1], beta_1), z_fun(raw[,5], capprox_harrell),
           z_fun(raw[,9], true_beta_bin), z_fun(raw[,13], mean_diff_true))
  
  if (runvglm == 1) {
    npo_bias <- colMeans(raw[, idx_npo_est], na.rm = TRUE) - true_grp_coeffs
    znpo <- sapply(1:3, \(i) z_fun(raw[, idx_npo_est[i]], true_grp_coeffs[i]))
    npo_cov <- colMeans(raw[, idx_npo_cov], na.rm = TRUE)
    npo_pow_per_beta <- colMeans(raw[, idx_npo_pow], na.rm = TRUE)
    # NEW: Calculate power based on "any significant" flag
    npo_any_pow <- mean(raw[, idx_npo_any_pow], na.rm = TRUE)
    
    npo_err <- sum(raw[, idx_flag[1]], na.rm = TRUE); npo_conv <- mean(raw[, idx_flag[3]], na.rm = TRUE)
    good_lrt <- raw[, idx_flag[3]] == 1 & !is.na(raw[, idx_flag[2]])
    npo_lrt_pow <- mean(raw[good_lrt, idx_flag[2]], na.rm = TRUE)
  } else {
    npo_bias <- znpo <- npo_cov <- npo_pow_per_beta <- rep(NA, 3)
    npo_any_pow <- npo_err <- npo_conv <- npo_lrt_pow <- NA
  }
  
  results_list <- list(
    n = n, beta1_OR = exp(beta_1),
    poregbias = po[1], poregcov = po[2], poregp = po[3],
    wpbias = wp[1], wpcov = wp[2], wpp = wp[3],
    bin7bias = bin[1], bin7cov = bin[2], bin7p = bin[3],
    linbias = lin[1], lincov = lin[2], linp = lin[3],
    pobiasz = zb[1], wpbiasz = zb[2], bin7biasz = zb[3], linbiasz = zb[4],
    
    nonPO_d0_bias = npo_bias[1], nonPO_d15_bias = npo_bias[2], nonPO_d30_bias = npo_bias[3],
    nonPO_d0_Zbias = znpo[1], nonPO_d15_Zbias = znpo[2], nonPO_d30_Zbias = znpo[3],
    nonPO_d0_cov = npo_cov[1], nonPO_d15_cov = npo_cov[2], nonPO_d30_cov = npo_cov[3],
    nonPO_d0_pow = npo_pow_per_beta[1], nonPO_d15_pow = npo_pow_per_beta[2], nonPO_d30_pow = npo_pow_per_beta[3],
    
    sumnonpoErrFlag = npo_err,
    nonPO_pow = npo_any_pow, # NEW Power Column
    nonpoLRTpower = npo_lrt_pow,
    nonpoConvProp = npo_conv
  )
  
  dfout <- as.data.frame(results_list)
  dfout$runvglm <- runvglm
  dfout$n_simulations <- n_simlist
  
  # UPDATED: Add the new power column to the power summary data frame
  dfpow <- dfout[, c("n","beta1_OR","runvglm","n_simulations",
                     "poregp","wpp","bin7p","linp",
                     "nonpoLRTpower", "nonPO_pow"), drop = FALSE]
  names(dfpow) <- c("TotalSampleSize","TrueOR","RunVGAMflag", "NumSimulations",
                    "PowerPO","PowerWinP", "PowerBinaryReg","PowerLinearReg",
                    "PowerNonPOLRT", "PowerNonPOAnySig")
  
  list(full  = cbind(ident, dfout),
       power = cbind(ident, dfpow))
}


# ==============================================================================
# 2. GRID DRIVER
# ==============================================================================
run_simulation_grid <- function(
    beta_1list          = c(log(1.0)),
    nonPOlist           = c(1),
    attenuateeffectlist = c(1),
    symmetric           = c(0,1),
    nlist               = c(1500),
    n_simlist           = 1000,
    ybinhigh            = 7,
    runvglm             = 1,
    seed                = 14424) {
  
  grid <- expand.grid(beta_1 = beta_1list,
                      nonPO  = nonPOlist,
                      atten  = attenuateeffectlist,
                      sym    = symmetric,
                      n      = nlist,
                      KEEP.OUT.ATTRS = FALSE)
  
  future::plan(future::multisession,
               workers = max(1, future::availableCores() - 2))
  
  progressr::handlers("txtprogressbar")
  
  full_rows  <- vector("list", nrow(grid))
  power_rows <- vector("list", nrow(grid))
  
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid))
    
    for (i in seq_len(nrow(grid))) {
      pars <- grid[i, ]
      cat(sprintf("\nScenario %d/%d : logOR=%.3f  nonPO=%d  atten=%d  sym=%d  n=%d\n",
                  i, nrow(grid), pars$beta_1, pars$nonPO,
                  pars$atten, pars$sym, pars$n))
      
      res <- run_one_scenario(beta_1      = pars$beta_1,
                              nonPO       = pars$nonPO,
                              atten       = pars$atten,
                              sym_flag    = pars$sym,
                              n           = pars$n,
                              n_simlist   = n_simlist,
                              ybinhigh    = ybinhigh,
                              runvglm     = runvglm,
                              seed        = seed)
      full_rows[[i]]  <- res$full
      power_rows[[i]] <- res$power
      p() # update progress
    }
  })
  
  summary_df <- do.call(rbind, full_rows)
  power_df   <- do.call(rbind, power_rows)
  
  save_csv_with_fallback(summary_df, "simulation_summary_grid_weak_po_violation.csv")
  save_csv_with_fallback(power_df,   "simulation_power_grid_weak_po_violation.csv")
  
  beepr::beep("fanfare")
  invisible(list(summary = summary_df, power = power_df))
}



# ==============================================================================
# 3. RUN THE SIMULATION
# ==============================================================================
# Example run with a comprehensive grid of parameters.
# NOTE: This will take a significant amount of time to complete.
# For a quick test, reduce n_simlist to a small number (e.g., 10 or 100).
run_simulation_grid(beta_1list          = c(log(0.5),log(0.6),log(0.7),log(0.8),log(0.9),log(1.0)),
                    nonPOlist           = c(0, 1),
                    attenuateeffectlist = c(0, 1),
                    symmetric           = c(0, 1),
                    nlist               = c(300,600,1200,2400),
                    n_simlist           = 1000, # Using 100 for a quick test run
                    ybinhigh            = 7,
                    runvglm             = 1,
                    seed                = 14424)