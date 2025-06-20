# Simulation Study for Ordinal Outcomes in RCTs under the PO and Non-PO Models

This repository contains R scripts for simulating ordinal outcome data in randomized controlled trials (RCTs) under both Proportional Odds (PO) and Non-Proportional Odds (Non-PO) assumptions. The simulation adheres to the ADEMP framework, as described by Morris et al. (2019).

---

## Overview

We evaluated the performance of parametric and non-parametric methods for analyzing ordinal outcomes in individual-level randomized controlled trials. Our simulations cover scenarios with both adherence to and violation of the PO assumption.

---

## Repository Contents

- `Simulated Data Generation for Ordinal Outcomes in RCTs.R`: Main R script containing the functions and simulation workflow.
- `requirements.txt`: Lists R version and required packages with specific versions.

---

## Required Packages

Install necessary packages using:

```R
install.packages("devtools")
devtools::install_github("proshano/winprob")

install.packages(c("rje", "rms", "MASS", "dplyr", "lmtest", "sandwich", "VGAM", "beepr", "future", "future.apply", "progressr"))
```

---

## Simulation Configuration

### Adjusting Parameters:

- **Ordinal Outcome Distribution**:
  - Modify `pylist` vector in the script to reflect desired outcome probabilities.

- **Model Type (PO vs. Non-PO)**:
  - Set `nonPOlist`:
    - `0`: PO assumption holds
    - `1`: Non-PO scenario

  - If `nonPOlist` = `1`, set `attenuateeffectlist`:
    - `0`: PO violation strengthens treatment effect of the common odds ratio
    - `1`: PO violation attenuates treatment effect of the common odds ratio

- **Sample Size**:
  - Adjust `nlist` to desired total sample size (e.g., 300, 600, 1200, 2400).

---

## Running the Simulation

Execute the main script:

```R
source("20250609 - Parallel Processing Sim- Strong PO Violation - Ordinal Outcomes in RCTs.R")
```

Results will be saved as a CSV file in your working directory.

---

## Results Interpretation

Output CSV provides:
- Bias, coverage probability, and p-values for:
  - Proportional odds logistic regression
  - Non-parametric win probability (winP)
  - Binary logistic regression
  - Linear regression
  - Partial proportional odds logistic regresssion

---

## Citation

If you use this code in your research, please cite our manuscript:

> TBD [Currently under review.)

---

For issues or questions, please open an issue on GitHub.

