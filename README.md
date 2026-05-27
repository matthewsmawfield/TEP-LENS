# Temporal Equivalence Principle: A Blind-Prediction Residual Test in Multiply-Imaged Supernovae

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Lisboa)  
**First published:** 02 March 2026  
**Status:** Preprint (Open for Collaboration)  
**Website:** [https://matthewsmawfield.github.io/TEP-LENS/](https://matthewsmawfield.github.io/TEP-LENS/)  
**Paper Series:** TEP Series: Paper 19 (Strong Lensing Time Delays)

## Abstract

SN Refsdal provides a rare blind-prediction test of potential-dependent temporal propagation: seven GR lens-modelling teams predicted the SX reappearance delay before the image was observed, while the later Kelly et al. (2023) measurement provides an independent comparator. This single-system probe reduces to a single-contrast measurement: signal-energy partitioning (Step 35) shows that 99.9% of the predicted proxy-model signal resides in the long-baseline S4--SX contrast, with an effective dimensionality $D_{\rm eff} \approx 2.0$; the inner Einstein cross serves as a *null region* providing no probative leverage. The underlying algebraic loop identity is insensitive to a uniform Mass Sheet Degeneracy -- a uniform convergence sheet rescales all delays symmetrically and cannot generate a differential residual -- but the empirical blind-prediction residual remains limited by GR lens-model precision and correlated modelling systematics.
The pre-specified primary non-parametric directional test is a Wilcoxon signed-rank test on the six non-zero residuals among the seven blind models, all of which are positive ($p = 0.016$, approximately $2.2\sigma$), matching the log-magnification proxy prediction for a negative temporal-shear coupling. The supplementary all-eight-model Wilcoxon (including one post-blind precision update) gives $p = 0.0078$ ($2.4\sigma$). Additional strands are directionally consistent: a 46% reduction in weighted RMS after proxy-model correction, and a binomial sign test of 6/7 positive among blind models ($p = 0.0625$). Hierarchical Bayesian model comparison remains inconclusive (BF $= 1.06$ fixed-$\alpha$ baseline / $0.997$ h0pe-informed, BF $= 0.615$ free-$\alpha$ baseline / $0.492$ h0pe-informed with a GR-centred prior), indicating that present lens-model uncertainties dominate formal model-selection metrics. The measured coupling $\alpha_{\rm lens} = -0.055 \pm 0.044$ is calibrated from the same SN Refsdal data, so the magnitude agreement with the proxy-model prediction is definitional rather than an independent confirmation; the probative content lies in the sign consistency across modelling groups that use independent codes but share the same lens and image constraints. In the absence of a solved TEP lensing transfer function for $\kappa(\boldsymbol{\theta})$ or $\psi(\boldsymbol{\theta})$, this test should be read as a falsifiable phenomenological screen, not a fundamental coupling measurement. Internal consistency checks, including a numerical $H_0$ exercise confined to SN Refsdal, are discussed in Section 4.8. Additional observed strands are directionally consistent: weighted residual $\mathcal{R}_{\rm obs}=+14.6\pm11.6$ d (within $0.01\sigma$ of the TEP prediction), delay-magnification correlation ($r=0.932$, SX-driven), and 46% wRMS reduction after TEP correction. The extended TDCOSMO-2025 and SN Encore dataset shows a predicted systematic scaling trend consistent with TEP at the measured coupling, though these systems do not permit an independent closure test. Robustness analyses show stability to model dependence and to 10--30% flux-proxy perturbations, while hierarchical Bayesian model comparison remains inconclusive at current uncertainties (non-decisive Bayes factors: BF ≈ 1.00 for the fixed-alpha test, BF ≈ 0.49 for the free-alpha model with a GR-centred prior). The evidence is therefore strong in directional terms but not yet model-selection decisive: consistent support for TEP, with model-selection-level discrimination limited by present lens-model errors. Recognizing that the various evidence strands for SN Refsdal are structurally correlated, invalid meta-analytic combinations are avoided. The Wilcoxon signed-rank test was pre-registered as the primary non-parametric directional test because it equal-weights independent modelling groups, eliminating inverse-variance downweighting from heterogeneous uncertainties. All 7 non-zero residuals are positive ($p=0.0078$, $z=2.4\sigma$). Directional-odds metrics are treated as interpretive complements. A cross-system consistency check shows that applying the measured TEP coupling ($\alpha_{\rm lens} \approx -0.055$) provides a candidate resolution for the "Low $H_0$" bias in SN Refsdal (shifting from $66.6$ to $69.3$ km s$^{-1}$ Mpc$^{-1}$), though this is not an independent confirmation because $\alpha_{\rm lens}$ was calibrated on the same data.

## Key Results

- **Wilcoxon signed-rank (primary, blind 7):** 6/6 non-zero residuals positive, $p = 0.016$ ($\approx 2.2\sigma$)
- **Wilcoxon signed-rank (supplementary, all 8):** 7/7 non-zero positive, $p = 0.0078$ ($2.4\sigma$)
- **Weighted blind-prediction residual:** $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d ($p = 0.10$, $1.26\sigma$; definitional, not independent)
- **Single-contrast dominance:** 99.9% signal energy in S4--SX, $D_{\rm eff} \approx 2.0$
- **Cross-system $H_0$:** internal consistency check only; not independent evidence
- **Robustness:** model dependence, microlensing MC (10–30%), hierarchical Bayes, external H0LiCOW/TDCOSMO chains

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102062](https://doi.org/10.5281/zenodo.19102062) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454863](https://doi.org/10.5281/zenodo.19454863) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |
| **Paper 18** | [TEP-HC](https://github.com/matthewsmawfield/TEP-HC) | EFT Mapping and Acoustic Peak Constraints via hi_class | — |
| **Paper 19** | **TEP-LENS** (This repo) | Blind-Prediction Residual Test in Multiply-Imaged Supernovae | — |

## Directory Structure

```text
TEP-LENS/
├── data/
│   ├── raw/                 # SN Refsdal, H0pe, TDCOSMO catalogs
│   ├── interim/             # Pipeline intermediates
│   └── cosmograil/          # CosmoGRAIL inputs (when used)
├── logs/                    # Step execution logs
├── manuscripts/             # Generated markdown (from site build)
├── results/                 # Figures and JSON outputs
├── scripts/
│   ├── steps/               # Numbered analysis pipeline
│   │   └── run_all_steps.py
│   └── utils/               # Shared utilities
├── site/
│   └── components/          # HTML source of truth for manuscript
├── README.md
├── CITATION.cff
├── VERSION.json
├── version.txt
├── zenodo.txt
└── requirements.txt
```

## Installation

```bash
git clone https://github.com/matthewsmawfield/TEP-LENS.git
cd TEP-LENS
pip install -r requirements.txt
```

## Reproduction Pipeline

```bash
# Full pipeline (steps 00–20, 21 total)
python scripts/steps/run_all_steps.py

# Build manuscript from HTML components (static site + markdown)
cd site && npm ci && npm run build
# Output: 19-TEP-LENS-v0.1-Lisboa.md (repo root and manuscripts/)

# Generate PDF (requires playwright: pip install playwright && playwright install chromium)
python scripts/generate_site_pdf.py --quality high --wait-time 5
# Output: site/public/docs/19-TEP-LENS-v0.1-Lisboa.pdf and repo root copy

# Deploy static site
./deploy.sh
```

## Citation

```bibtex
@article{tep_lens_paper,
  title={Temporal Equivalence Principle: A Blind-Prediction Residual Test in Multiply-Imaged Supernovae},
  author={Smawfield, Matthew Lukin},
  year={2026},
  note={Preprint v0.1 (Lisboa)},
  url={https://github.com/matthewsmawfield/TEP-LENS}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons licenses to encourage replication. Feedback and collaboration are warmly invited.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
