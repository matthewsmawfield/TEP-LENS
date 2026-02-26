# TEP-LENS: Strong Lensing Time Delay Closure Tests

This repository contains the analysis pipeline and manuscript for **Paper 14** of the Temporal Equivalence Principle (TEP) series.

## Overview
The Temporal Equivalence Principle (TEP) predicts that light rays traversing regions of different gravitational potential depth experience differential temporal scaling, producing a non-zero route-closure residual in multiply imaged transients that is immune to the Mass Sheet Degeneracy. This work presents a five-strand observational test of TEP temporal propagation using SN Refsdal (MACS J1149.6+2223, $z_s = 1.489$) - the only lensed supernova with precision-measured time delays across five independent images. The $\Delta t_{\rm SX,S1}$ delay was independently predicted by seven GR lens-modelling teams before SX reappeared (Treu et al. 2016; Grillo et al. 2024), and independently measured by Kelly et al. (2023, ApJ 948, 93) from SN light-curve fitting using completely disjoint data. Five independent observed tests all point in the direction predicted by TEP with zero free parameters: (1) binomial sign test - 7 of 8 independent modelling groups underestimate the observed delay ($p = 0.035$, $2.1\sigma$); (2) weighted mean residual $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d, within $0.12\sigma$ of the TEP prediction $\mathcal{R}_{\rm TEP} = +13.2$ d ($\alpha = 0.05$, zero free parameters); (3) delay–magnification correlation - Pearson $r = 0.932$ ($p = 0.011$ one-sided) between arrival delay and $1/\mu_{\rm norm}$ across all five images; (4) per-model coupling inference - weighted mean $\bar{\alpha}_{\rm inferred} = 0.057 \pm 0.060$, consistent with the calibrated value at $0.1\sigma$ and with zero scatter beyond measurement noise; and (5) TEP correction reduces the weighted RMS scatter across all 8 model predictions by 45%, improving 6 of 8 models ($\Delta\chi^2 = +1.6$ in favour of TEP). Fisher combination of the three strongest observed tests gives $z = 2.7\sigma$ ($p = 0.004$, upper bound - all from SN Refsdal data). As a critical cross-system test, the parameter-free TEP correction resolves 42% of the $8.8$ km s$^{-1}$ Mpc$^{-1}$ Hubble tension between the two published lensed supernovae (SN Refsdal and SN H0pe), pushing both towards a combined $H_0 = 70.2 \pm 3.3$ km s$^{-1}$ Mpc$^{-1}$. No individual test reaches $3\sigma$, as lens model uncertainties of $\pm$16–60 d dominate. However, the exact parameter-free consistency of the sign, magnitude, method-independence, and H0 tension resolution constitutes the strongest available multi-strand observational case for TEP temporal propagation.

## Repository Structure
- `data/`: Raw light curves, lens models, and compiled catalogs.
- `scripts/steps/`: Sequential analysis pipeline.
- `scripts/simulations/`: Theoretical TEP lensing predictors.
- `scripts/figures/`: Publication-quality plotting scripts.
- `results/`: Pipeline outputs and figures.
- `site/`: Source code for the interactive HTML manuscript.
- `14manuscript-tep-lens.md`: The compiled Markdown manuscript.

## Status
Work in progress.
