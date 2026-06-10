# TEP-LENS Target Monitoring Tracker

Last updated: 2026-05-28

## Active Systems with Blind-Prediction Residual Tests

| System | Images | z_src | z_lens | Blind Models | Time Delays | TEP Evidence | Status |
|--------|--------|-------|--------|--------------|-------------|--------------|--------|
| SN Refsdal | 5 (S1,S2,S3,S4,SX) | 1.49 | 0.54 | 7 (+1 post-blind) | Measured (Kelly+2023) | **z ≈ +2.4, p = 0.0078** (blind Wilcoxon) | Primary evidence strand |
| SN Encore | 2 (1a,1b) | 1.95 | 0.338 | 8 | Measured (Pierel+2026) | Consistency check only; predicted shift ~0.5 d swamped by scatter | Secondary check |
| SN H0pe | 3 (A,B,C) | 1.783 | 0.351 | 7 | Measured (Pierel+2024) | Consistency check only; predicted shift ~1-2 d swamped by scatter | Secondary check |

## High-Priority Future Targets

### SN 2025wny (SN Winny)
- **Host:** PS1J0716+3821 (z_lens = 0.3754)
- **Source:** z = 2.011
- **Images:** 4 (A, B, C, D) — Einstein cross pattern
- **Magnification:** mu_A ~ 20-50 (extremely high contrast)
- **Predicted time delays** (post-hoc Witt-Wynne model, arXiv:2605.11090):
  - A trails C by ~20 d
  - A trails B by ~30 d
  - A trails D by ~175 d (long baseline!)
- **Why it matters for TEP:**
  - 4 images → closed 3-image loops possible (A-B-C, A-B-D, A-C-D, B-C-D)
  - High magnification contrast → large predicted TEP shifts (comparable to Refsdal S4-SX)
  - Long delay baseline (~175 d) → high SNR
- **What is needed:**
  1. Time delay measurements from ongoing photometric monitoring
  2. BLIND lens model predictions (published BEFORE delays are measured)
- **Ongoing programs:**
  - HST PID 17611 (PI: Goobar)
  - JWST PID 5564 (PI: Goobar)
  - Ground-based: Maidanak, Lulin, COLIBRI, Wendelstein
- **Literature:**
  - Discovery: Johansson et al. 2025, ApJ 995, L17
  - HOLISMOKES characterization: Taubenberger et al. 2025 (arXiv:2510.21694)
  - LSST protocol: arXiv:2605.11090
  - Host lens search: Cañameras et al. 2020 (HOLISMOKES II)
- **TEP Priority:** **HIGH** — Most promising system after Refsdal

## Lower-Priority / Non-Viable Targets

### SN Eos
- **Host:** MACS J1931.8-2635 (z_lens = 0.35)
- **Source:** z = 5.133 (farthest spectroscopically confirmed SN)
- **Images:** 2 confirmed
- **Why not viable for TEP:** Only 2 images → no closed 3-image loop possible
- **Literature:**
  - Discovery: Coulter et al. 2026 (arXiv)
  - VENUS lens model: arXiv:2602.14074
- **TEP Priority:** LOW

## What Makes a System Viable for TEP?

1. **Minimum 3 images** (preferably 4+) to form closed algebraic loops
2. **High magnification contrast** between images (mu_max/mu_min > 5)
3. **Long delay baselines** (>50 d) for high SNR
4. **Blind lens model predictions** — models published BEFORE time-delay measurements
5. **Time delay measurements** from light-curve photometry

## Monitoring Actions

- [ ] Check for SN 2025wny time-delay measurement papers (quarterly)
- [ ] Check for blind lens-model predictions for PS1J0716+3821 (check Cañameras+2020, any HOLISMOKES follow-up)
- [ ] Monitor arXiv for new multiply-imaged SN discoveries with 4+ images
- [ ] Monitor JWST/HST proposal databases for new lensed SN follow-up programs
- [ ] Set up Google Scholar alerts for: "multiply-imaged supernova" + "time delay" + "lens model"
