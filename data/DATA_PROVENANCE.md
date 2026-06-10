# TEP-LENS Data Provenance

This document provides complete provenance information for all external data used in the TEP-LENS analysis pipeline.

## 1. Multiply-Imaged Supernova and Time-Delay Data

**Primary system:** SN Refsdal (MACS J1149.6+2223)  
**Sources:** Kelly et al. (2023, ApJ 948, 93) — measured pairwise delays (Table 15, "Combined" method) and magnification proxies (total flux ratios); original modeling papers — pre-reappearance blind GR model predictions (per-team references in scripts/steps/step_07_observed_vs_predicted.py); Grillo et al. (2024, ApJ 971, 49) — post-blind precision update  
**Ingestion:** `scripts/steps/step_01_fetch_snh0pe_data.py` compiles values from published literature into `data/raw/sn_lensing/lensed_sn_catalog.json`

**Secondary system:** SN H0pe (PLCK G165.7+67.0)  
**Sources:** Pierel et al. (2024, ApJ 967, 50); Frye et al. (2024, ApJ 961, 171); Grayling et al. (2025, arXiv:2510.11719)  
**Ingestion:** Same catalog file via step_01

**Prediction target:** SN 2025wny (z=2.011)  
**Source:** Johansson et al. (2025, ApJ 995, L17)  
**Ingestion:** Same catalog file via step_01

**Additional target:** SN Encore  
**Source:** Pierel et al. (2026, ApJ, arXiv:2509.12301)  
**Ingestion:** Same catalog file via step_01

---

## 2. TDCOSMO Quad-Lens Quasar Time Delays

**Source:** TDCOSMO Collaboration (2020–2025) — published time-delay measurements for quad-lens quasar systems  
**Systems:** HE0435-1223, WFI2033-4723, RXJ1131-1231, PG1115+080, and others  
**Data format:** `.rdb`, `.dat`, and `.csv` files in `data/cosmograil/`  
**Ingestion:** `scripts/steps/step_05_tdcosmo_shear.py` reads published delay and flux-ratio values from `TDCOSMO_QUADS` dictionary (hardcoded from literature)

---

## 3. H0LiCOW Distance Chains

**Source:** H0LiCOW Collaboration — posterior MCMC chains for time-delay cosmography  
**Ingestion:** `scripts/steps/step_14_external_chain_ingestion.py` discovers and ingests chain files from `data/interim/external/h0licow_distance_chains/`

---

## 4. TDCOSMO 2025 Public Data

**Source:** TDCOSMO Collaboration 2025 release  
**Ingestion:** `scripts/steps/step_19_tdcosmo2025_ingestion.py` processes data from `data/interim/external/tdcosmo2025_public/`

---

## 5. Lensed SN Registry

**Source:** arXiv preprint registry and literature survey  
**File:** `data/interim/external/arxiv_lensed_sn_registry.json`  
**Ingestion:** `scripts/steps/step_00_fetch_literature_and_cross_paper_data.py`

---

## Verification Commands

```bash
# Run full pipeline with data regeneration
python scripts/steps/run_all_steps.py

# Check catalog contents
python -c "import json; d=json.load(open('data/raw/sn_lensing/lensed_sn_catalog.json')); print(list(d.keys()))"
```

---

## Data Quality Flags

| Dataset | Source | Method | Verified |
|---------|--------|--------|----------|
| SN Refsdal delays | Kelly+2023 ApJ 948, 93 | Manual transcription from Table 15 | ✓ |
| SN Refsdal flux ratios | Kelly+2023 ApJ 948, 93 | Manual transcription | ✓ |
| Blind model predictions | Treu+2016 ApJ 817, 60; Kelly+2023 Sci Supp | Manual transcription from Tables 2, S4 | ✓ |
| Grillo+2024 update | Grillo+2024 ApJ 971, 49 | Manual transcription | ✓ |
| SN H0pe delays | Pierel+2024 ApJ 967, 50; Grayling+2025 arXiv:2510.11719 | Manual transcription | ✓ |
| TDCOSMO quasar delays | TDCOSMO Collaboration | Hardcoded from published values | ✓ |
| TDCOSMO 2025 | TDCOSMO Collaboration 2025 | Automated ingestion | ✓ |
| H0LiCOW chains | H0LiCOW Collaboration | Automated discovery | ✓ |

---

*Last updated: 2026-05-25*
