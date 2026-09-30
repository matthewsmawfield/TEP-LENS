# Temporal Equivalence Principle: A Blind-Prediction Residual Test in Multiply-Imaged Supernovae

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20572720.svg)](https://doi.org/10.5281/zenodo.20572720)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Lisboa)  
**First published:** 10 June 2026 · **Last updated:** 30 September 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.20572720](https://doi.org/10.5281/zenodo.20572720)  
**Website:** [https://mlsmawfield.com/tep/lens](https://mlsmawfield.com/tep/lens)  
**Paper Series:** TEP Series: Paper 19 (Strong Lensing Time Delays)

## Abstract





The Temporal Equivalence Principle modifies gravitational lensing not through the photon path but through mass inference. The conformal transport of photons produces exactly zero additional delay (null-curve invariance); the scalar field's backreaction on the photon path, evaluated from the common action (Step 61), is foreclosed at ~ 2 × 10³ below the required phantom mass; and the disformal contribution, evaluated on the cluster path with the GW170817 constraint applied to the coupling B_0 rather than to the local cone tilt (Step 62), falls a factor ≈ 5 short of the observed residual under a genuine deprojection of the GLAFIC convergence map — B_0 req ≈ 420 against the bound B_0 ≲ 78 at the best-fit NFW halo scale, reaching the envelope only for coherent line-of-sight extents ≳ 1.8 Mpc well beyond standard cluster concentrations (Step 63). The surviving channel is the lensing–dynamical mass inference: stellar kinematics can break the lens-model mass-sheet degeneracy, making the dynamical–lensing mass difference — which measures the phantom mass — a concrete, sign-testable consequence of the coupled reconstruction, and the image-plane calculation below is a forward force-response realization and target scale for that mechanism. SN Refsdal offers an unusually high-leverage case. Seven pre-reappearance lens-model variants spanning five modelling families published predictions for the long-baseline SX reappearance delay before the image was observed; Kelly et al. (2023) later measured that same delay independently from SN light-curve fitting. Six of seven strictly blind original predictions and all seven delay-blind revised predictions yield positive residuals (observed delay longer than predicted). Under a uniform time-delay rescaling, the +30.1 d residual corresponds to the required target scale shift of approximately 8.7%. This is the measured response that the coupled reconstruction must reproduce, with external-convergence and mass-model uncertainty propagated in the amplitude fit. The correlation-aware family-sign-flip test gives p = 0.031 on the delay-blind ensemble. The evidence establishes a single-system, sign-only directional result; the decisive amplitude test is the preregistered geometric forecast for the next long-baseline multiply-imaged supernova (SN 2025wny). The remaining theoretical task is narrow: the action-level scalar backreaction and screening law are already specified and tested in their native benchmarks. The open closure is their lens-domain projection through one solved lens geometry, with the GR halo contribution kept separate and the absolute delay amplitude fixed without fitting the Refsdal residual. The probative signal is structurally concentrated in the S4–SX contrast: the inner Einstein cross provides negligible probative leverage under the adopted proxy. A transport audit establishes that the archived ψ map measures the ordinary Fermat/Shapiro potential already included in the GR lens predictions; because tilde g_μν=A²g_μν preserves null curves exactly, it cannot be counted again as a static conformal residual. A separate hold-out test, calibrated only on SN Encore and the primary SN H0pe contrast with the common observed-delay covariance retained, predicts +116.8 ± 74.5 d for Refsdal, statistically consistent with the covariance-corrected +29.6 ± 10.2 d at 1.16σ but too imprecise and centrally too large to close the amplitude channel. The blind sign result survives; with the direct backreaction kernel now evaluated and foreclosed (Step 61) and the disformal transport bounded ~ 5 ×  short under the measured deprojection (Step 63), a quantitative amplitude prediction requires the surviving channel — the phantom mass carried by the mass-sheet/dynamics inference — to be derived from the action before application to the next long-baseline system.


## Key Results


- **Family-sign-flip (headline, delay-blind 7):** correlation-aware exact test, $p = 0.031$ ($\approx 1.86\sigma$)
- **Wilcoxon signed-rank (benchmark, delay-blind 7):** 7/7 non-zero residuals positive, $p = 0.0078$ ($\approx 2.4\sigma$)
- **Wilcoxon signed-rank (supplementary, all 8):** 8/8 non-zero positive, $p = 0.0039$ ($2.8\sigma$)
- **Weighted blind-prediction residual:** $\mathcal{R}_{\rm obs} = +30.1 \pm 8.9$ d ($p = 0.0006$, $3.39\sigma$; definitional, not independent)
- **Single-contrast dominance:** 99.9% signal energy in S4--SX, $D_{\rm eff} \approx 2.0$
- **Cross-system $H_0$:** internal consistency check only; not independent evidence
- **Robustness:** model dependence, microlensing MC (10–30%), hierarchical Bayes, external H0LiCOW/TDCOSMO chains

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |
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
# Full pipeline (36 registered steps: 00-20 plus extended 30-42 diagnostics)
python scripts/steps/run_all_steps.py

# Build manuscript from HTML components (static site + markdown)
cd site && npm ci && npm run build
# Output: 19-TEP-LENS-v0.2-Lisboa.md (repo root and manuscripts/)

# Generate PDF (requires playwright: pip install playwright && playwright install chromium)
python scripts/generate_site_pdf.py --quality high --wait-time 5
# Output: site/public/docs/19-TEP-LENS-v0.2-Lisboa.pdf and repo root copy

# Deploy static site
./deploy.sh
```

## Known Data Caveats

- **GLAFIC v3 map normalization offset (§3.3.1):** The archived GLAFIC v3 lensing-potential and convergence maps (used in Steps 44–54) exhibit a ~4 ×  normalization offset relative to the Kelly et al. (2023) tabulated parameters. Diagnostic checks show this is not a simple source-redshift rescaling ($D_{\rm ls}/D_{\rm s}$ and $\Sigma_{\rm crit}$ ratios do not account for the discrepancy). The most probable origin is either a localized pixel-to-arcsecond unit-translation error in the raw archived map files or an undocumented mass-sheet projection effect inside the GLAFIC inversion coordinate system. Neither invalidates the tabulated parameters, but researchers re-running the 3D Abel deprojection or Jacobian transfer-kernel reconstruction (Step 54) should treat the map-derived absolute amplitudes as exploratory and rely on the Kelly+2023 tabulated values for quantitative comparison.

## Citation

```bibtex
@article{tep_lens_paper,
  title={Temporal Equivalence Principle: A Blind-Prediction Residual Test in Multiply-Imaged Supernovae},
  author={Smawfield, Matthew Lukin},
  year={2026},
  note={Preprint v0.2 (Lisboa)},
  url={https://github.com/matthewsmawfield/TEP-LENS}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons licenses to encourage replication. Feedback and collaboration are warmly invited.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)