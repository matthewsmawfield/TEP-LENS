# The Temporal Equivalence Principle: A Geometric Route-Closure Test in Multiply-Imaged Supernovae

This repository contains the analysis pipeline and manuscript for **Paper 14** of the Temporal Equivalence Principle (TEP) series.

## Overview
The Temporal Equivalence Principle (TEP) predicts that light rays traversing different gravitational potential depths experience differential temporal scaling ($\alpha \neq 0$). This work applies the test to SN Refsdal, showing that an expansion-mode coupling ($\alpha = -0.05$) resolves the "Low $H_0$" bias in lensed supernovae. The key observed comparison uses disjoint data: seven blind pre-reappearance GR predictions of $\Delta t_{\rm SX,S1}$ versus the later Kelly et al. (2023) measurement. The strongest single observed result is a Wilcoxon signed-rank test: all 7 non-zero model residuals are positive ($p = 0.0078$, 2.4$\sigma$), matching the TEP expansion prediction. Additional observed strands are directionally consistent: weighted residual $\mathcal{R}_{\rm obs}=+14.6\pm11.6$ d (within $0.12\sigma$ of the TEP prediction), delay-magnification correlation ($r=0.932$, SX-driven), and 45% wRMS reduction after TEP correction. The extended TDCOSMO-2025 and SN Encore dataset independently corroborates this systematic effect. Robustness analyses added here show stability to model dependence and to 10-30% flux-proxy perturbations, while hierarchical Bayesian model comparison remains inconclusive at current uncertainties (Bayes factors $\approx 1$). The evidence is therefore strong in directional terms but not yet model-selection decisive: consistent support for TEP, with model-selection-level discrimination limited by present lens-model errors. Recognizing that the various evidence strands for SN Refsdal are structurally correlated, invalid meta-analytic combinations are avoided. The most defensible approach is to select the single most robust non-parametric test (Wilcoxon signed-rank, all 7 non-zero residuals positive) as the conservative headline significance for the system ($p=0.0078$, $z=2.4\sigma$). Directional-odds metrics are treated as interpretive complements. A cross-system consistency check shows that the parameter-free TEP correction ($\alpha=-0.05$) resolves the "Low $H_0$" bias in SN Refsdal, SN Encore, and SN H0pe (TD-only), shifting them from $\sim 61-66$ to $\sim 63-69$ km s$^{-1}$ Mpc$^{-1}$, improving agreement with Planck.

## Repository Structure
- `data/`: Raw light curves, lens models, and compiled catalogs.
- `scripts/steps/`: Sequential analysis pipeline.
- `scripts/simulations/`: Theoretical TEP lensing predictors.
- `scripts/figures/`: Publication-quality plotting scripts.
- `results/`: Pipeline outputs and figures.
- `site/`: Source code for the interactive HTML manuscript.
- `manuscripts/14manuscript-tep-lens.md`: The compiled Markdown manuscript.

## Status
Work in progress.

## Latest Robustness Upgrade (Steps 11-14)
- **Step 11:** model-dependence robustness (effective sample size, leave-one-out stress tests, exact sign-flip and Wilcoxon tests).
- **Step 12:** microlensing nuisance Monte Carlo (10-30% flux-proxy perturbations propagated to residuals and fit preference metrics).
- **Step 13:** hierarchical Bayesian GR-vs-TEP comparison with baseline and H0pe-2025-informed prior sensitivity scenarios.
- **Step 14:** public H0LiCOW/TDCOSMO chain ingestion bridge for immediate external-data integration (without waiting for new observations).

Current interpretation remains conservative: directional consistency with TEP is robust across stress tests, while decisive model-selection evidence is still limited by lens-model uncertainties.
