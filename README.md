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

install.packages(c("rje", "rms", "MASS", "dplyr", "lmtest", "sandwich", "VGAM", "beepr"))
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
    - `0`: Strengthen treatment effect
    - `1`: Attenuate treatment effect

- **Sample Size**:
  - Adjust `nlist` to desired total sample size (e.g., 300, 1000, 1500, 2000).

---

## Running the Simulation

Execute the main script:

```R
source("20250402 - Simulated Data Generation for Ordinal Outcomes in RCTs.R")
```

Results will be saved as a CSV file (`dfoutp - <sample-size>.csv`) in your working directory.

---

## Results Interpretation

Output CSV provides:
- Bias, coverage probability, and p-values for:
  - Proportional Odds regression
  - Non-parametric win probability (winP)
  - Binary logistic regression
  - Linear regression

---

## Citation

If you use this code in your research, please cite our manuscript:

> TBD [Currently under review.)

---

For issues or questions, please open an issue on GitHub.

