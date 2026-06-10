# Temporal Equivalence Principle: A Blind-Prediction Residual Test in Multiply-Imaged Supernovae

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)


## Abstract

The Temporal Equivalence Principle posits that the effective rate of temporal propagation scales with the depth of the local gravitational potential. Strong gravitational lensing provides a geometric test: a single source imaged through multiple sightlines accumulates differential temporal shear along each path, producing a blind-prediction residual between the observed delay and the delay predicted by a standard GR lens model. SN Refsdal, currently unique among multiply-imaged supernovae in having five resolved images and precision-measured delays, offers an unusually high-leverage case. Seven lens-model variants spanning five modelling families published blind predictions for the long-baseline SX reappearance delay before the image was observed; Kelly et al. (2023) later measured that same delay independently from SN light-curve fitting. All seven delay-blind model variants yield positive residuals (observed delay longer than predicted), matching the sign expected for a negative temporal-shear coupling. The directional evidence does not depend on amplitude calibration.

The probative signal is structurally concentrated in the S4–SX contrast: the inner Einstein cross provides negligible probative leverage under the adopted proxy, while the 376-day SX baseline amplifies the differential temporal-shear signature. Direct potential-map and 3D geodesic reconstructions preserve the residual sign but predict sub-day amplitudes, indicating that the corrected-compilation response scale points to a magnification-sensitive amplification mechanism. The operational log-magnification response captures the sign and the corrected-compilation order of magnitude, yet the amplitude match remains phenomenological: deriving the corresponding transfer kernel from the scalar-field action remains the central open theoretical task. The decisive next step is a prospective amplitude test on a future long-baseline multiply-imaged supernova.

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Lisboa)  
**First published:** 10 June 2026  
**Status:** In Development  
**Website:** [https://mlsmawfield.com/tep/lens](https://mlsmawfield.com/tep/lens)  
**Paper Series:** TEP Series: Paper 19 (Strong Lensing Time Delays)

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
# Full pipeline (36 registered steps: 00-20 plus extended 30-42 diagnostics)
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

## Known Data Caveats

- **GLAFIC v3 map normalization offset (§3.3.1):** The archived GLAFIC v3 lensing-potential and convergence maps (used in Steps 44–54) exhibit a ~4× normalization offset relative to the Kelly et al. (2023) tabulated parameters. Diagnostic checks show this is not a simple source-redshift rescaling ($D_{\rm ls}/D_{\rm s}$ and $\Sigma_{\rm crit}$ ratios do not account for the discrepancy). The most probable origin is either a localized pixel-to-arcsecond unit-translation error in the raw archived map files or an undocumented mass-sheet projection effect inside the GLAFIC inversion coordinate system. Neither invalidates the tabulated parameters, but researchers re-running the 3D Abel deprojection or Jacobian transfer-kernel reconstruction (Step 54) should treat the map-derived absolute amplitudes as exploratory and rely on the Kelly+2023 tabulated values for quantitative comparison.

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