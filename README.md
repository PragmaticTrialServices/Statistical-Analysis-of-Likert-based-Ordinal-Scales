# Simulation Study for Ordinal Outcomes in RCTs under the PO and Non-PO Models

This repository contains R scripts for simulating ordinal outcome data in randomized controlled trials (RCTs) under both Proportional Odds (PO) and Non-Proportional Odds (Non-PO) assumptions. The simulation adheres to the ADEMP framework, as described by Morris et al. (2019).

---

## Overview

We evaluated the performance of parametric and non-parametric methods for analyzing ordinal outcomes in individual-level randomized controlled trials. Our simulations cover scenarios with both adherence to and violation of the PO assumption.

---

## Repository Contents

- `20250609 - Parallel Processing Sim- Strong PO Violation - Ordinal Outcomes in RCTs.R`: Main R script containing the functions and simulation workflow.
- `requirements.txt`: Lists R version and required packages with specific versions.
- `Appendix1.R`: Online supplement file for our manuscript containing the R code for the MyTEMP PRO analysis.
- `simulation_summary_grid_strong_po_violation.csv`: CSV file containing the simulation results for 192 scenarios. Includes scenarios where the PO violation is strong.
- `simulation_summary_grid_weak_po_violation.csv`: CSV file containing the simulation results for 192 scenarios. Includes scenarios where there PO violation is weak.

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

## Column description for CSV files
* **`N`**: The total sample size per arm (i.e., \$2 \times n\$ per simulated trial).
* **`sym_tag`**: “Symmetric” vs. “Heavy0” label for the marginal outcome distribution.
* **`po_tag`**: “POMet” when the proportional-odds assumption holds, “POViolated” when it’s violated.
* **`eff_tag`**: “attenuatedeffects” vs. “strengthenedeffects”, indicating whether the non-PO deltas attenuated or strengthened the common treatment effect.
* **`beta1_OR`**: True odds ratio used in the PO model, i.e. $\exp(\beta_{k=0})$.
* **`poregbias`**: Bias of the PO model’s log-OR estimate: $\bar{\widehat\beta}_{PO}-\beta_{k=0}$.
* **`poregcov`**: Empirical coverage probability of the 95% CI for the PO log-OR (proportion of sims whose CI contained the true estimate).
* **`Poregp`**: Empirical power (or type-I error) of the PO Wald test: proportion of sims with \$p < .05\$ for the treatment effect in the PO model.
* **`Wpbias`**: Bias of the win-probability estimate: $\overline{\widehat{\mathrm{WinP}}}-\mathrm{WinP}_{\text{true}}$, where $\mathrm{WinP}_{\text{true}}$ is estimated using Harrell’s approximation.
* **`Wpcov`**: Coverage probability of the 95% CI for the win probability.
* **`Wpp`**: Empirical power of the win-probability test ($p < .05$).
* **`bin7bias`**: Bias of the binary-logistic regression coefficient at threshold: $\overline{\widehat\beta_{\text{binary}}}-\beta_{\text{binary,true}}$.
* **`bin7cov`**: Coverage probability of the 95% CI for that binary-logistic coefficient.
* **`bin7p`**: Empirical power of the binary-logistic regression.
* **`Linbias`**: Bias of the linear-regression estimate: $\overline{\widehat\beta\_{\text{linear}}}-\Delta\_{\text{true}}$, where $\Delta\_{\text{true}}$ is the true mean difference.
* **`Lincov`**: Coverage probability of the 95% CI for the linear-regression slope.
* **`Linp`**: Empirical power of the linear-regression test.
* **`Pobiasz`**: Relative bias (\%) for the PO estimate: $\displaystyle\frac{\overline{\widehat\beta_{\text{PO}}}-\beta_{\text{k=0}}}{\beta_{\text{k=0}}}\times 100%$.
* **`Wpbiasz`**: Relative bias (%) for the win-probability estimate (analogous to above).
* **`bin7biasz`**: Relative bias (%) for the binary-logistic estimate at $Y\ge7$.
* **`Linbiasz`**: Relative bias (%) for the linear-regression estimate: $\displaystyle\frac{\overline{\widehat\beta_{\text{lin}}}-\Delta_{\text{true}}}{\Delta_{\text{true}}}\times 100%$.
* **`sumnonpoErrFlag`**: Total count of simulation replicates where the VGAM non-PO fit flagged an error (i.e., model did not converge).
* **`nonpoLRTpower`**: Among converged non-PO fits, the empirical power of the likelihood-ratio test comparing non-PO vs. parallel PO (proportion with LRT $p < 0.05$).
* **`nonpoConvProp`**: Proportion of replicates where the VGAM non-PO model converged successfully.
## Citation

If you use this code in your research, please cite our manuscript:

> TBD [Currently under review.)

---

For issues or questions, please open an issue on GitHub.

