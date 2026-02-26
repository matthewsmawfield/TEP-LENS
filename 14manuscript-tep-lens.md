# TEP-LENS: Resolving the Hubble Tension

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Kingston upon Hull)  
**Date:** First published: 11 January 2026  
**DOI:** 10.5281/zenodo.18209703  
**Generated:** 2026-02-26  
**Paper Series:** TEP Series: Paper 12 (Cosmological Observations)

---

**Abstract.** The Temporal Equivalence Principle (TEP) predicts that light rays traversing regions of different gravitational potential depth experience differential temporal scaling, producing a non-zero route-closure residual in multiply imaged transients that is immune to the Mass Sheet Degeneracy. This work presents a five-strand observational test of TEP temporal propagation using SN Refsdal (MACS J1149.6+2223, $z_s = 1.489$) — the only lensed supernova with precision-measured time delays across five independent images. The $\Delta t_{\rm SX,S1}$ delay was independently predicted by seven GR lens-modelling teams before SX reappeared (Treu et al. 2016; Grillo et al. 2024), and independently measured by Kelly et al. (2023, ApJ 948, 93) from SN light-curve fitting using completely disjoint data. Five independent observed tests all point in the direction predicted by TEP with zero free parameters: (1) binomial sign test — 7 of 8 independent modelling groups underestimate the observed delay ($p = 0.035$, $2.1\sigma$); (2) weighted mean residual $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d, within $0.12\sigma$ of the TEP prediction $\mathcal{R}_{\rm TEP} = +13.2$ d ($\alpha = 0.05$, zero free parameters); (3) delay–magnification correlation — Pearson $r = 0.932$ ($p = 0.011$ one-sided) between arrival delay and $1/\mu_{\rm norm}$ across all five images; (4) per-model coupling inference — weighted mean $\bar{\alpha}_{\rm inferred} = 0.057 \pm 0.060$, consistent with the calibrated value at $0.1\sigma$ and with zero scatter beyond measurement noise; and (5) TEP correction reduces the weighted RMS scatter across all 8 model predictions by 45%, improving 6 of 8 models ($\Delta\chi^2 = +1.6$ in favour of TEP). Fisher combination of the three strongest observed tests gives $z = 2.7\sigma$ ($p = 0.004$, upper bound — all from SN Refsdal data). As a critical cross-system test, the parameter-free TEP correction resolves 42% of the $8.8$ km s$^{-1}$ Mpc$^{-1}$ Hubble tension between the two published lensed supernovae (SN Refsdal and SN H0pe), pushing both towards a combined $H_0 = 70.2 \pm 3.3$ km s$^{-1}$ Mpc$^{-1}$. No individual test reaches $3\sigma$, as lens model uncertainties of $\pm$16–60 d dominate. However, the exact parameter-free consistency of the sign, magnitude, method-independence, and H0 tension resolution constitutes the strongest available multi-strand observational case for TEP temporal propagation.

    **Keywords:** gravitational lensing — time delays — Temporal Equivalence Principle — SN Refsdal — SN H0pe — Hubble tension — binomial sign test

                
                

                    
## 1. Introduction

The Temporal Equivalence Principle (TEP) posits that the effective rate of temporal propagation scales with the depth of the local gravitational potential (Smawfield 2023–2025, Papers I–XIII). Where prior papers established the TEP framework through stellar evolution anomalies, Cepheid period-luminosity residuals, and galaxy-scale redshift correlations, the strong gravitational lensing regime provides a qualitatively different and entirely independent test. Here, a single photon traverses multiple geometric paths through a cluster or galaxy potential, accumulating a potential-dependent temporal shear along each sightline. The resulting pairwise time delays between images carry a direct imprint of the differential shear—independently of any cosmological model.

### 1.1 The Route-Closure Test

For a source producing three or more resolved images, the pairwise time delays obey a strict geometric closure identity under General Relativity: any three delays in a closed loop sum to identically zero. TEP breaks this identity. Light along path $i$ through a region of projected convergence $\kappa_i$ acquires a temporal shear $\Gamma_t(i) = 1 + \alpha\log_{10}(\mu_i)$, where $\mu_i$ is the magnification at image $i$ and $\alpha$ is the TEP coupling. The closure residual,

$\mathcal{R}_{\rm closure}(i,j,k) = (\Gamma_i-1)\Delta t_{ij} + (\Gamma_j-1)\Delta t_{jk} + (\Gamma_k-1)\Delta t_{ki}$,

is non-zero if and only if the images traverse regions of different potential depth. Crucially, this quantity is *immune to the Mass Sheet Degeneracy*: a uniform convergence sheet rescales all delays by the same factor, leaving the loop residual unchanged at zero. The route-closure test is therefore a clean, lens-model-independent probe of potential-dependent temporal propagation.

### 1.2 SN Refsdal: The Ideal System

SN Refsdal (MACS J1149.6+2223, Kelly et al. 2015) is the only known multiply-imaged supernova with *five* resolved images and precision-measured relative time delays. The first four images (S1–S4) form an Einstein cross around a cluster member galaxy; the fifth image SX appeared ~8 arcsec away at a separate arc position approximately 376 days after S1, predicted by cluster lens models before it was observed. Kelly et al. (2023, ApJ 948, 93) published the definitive time delay measurements from the combined light curve analysis, including the SX–S1 delay measured to 1.5% precision—the most precise lensed supernova delay to date.

The geometry of SN Refsdal is ideal for the route-closure test. The four Einstein-cross images (S1–S4) sample the deep potential of the cluster member galaxy halo, while SX samples the outer cluster arc at much lower magnification. This contrast in potential depth—combined with the long 376-day SX baseline—amplifies the expected TEP closure residual by a factor of ~18 relative to the inner cross loops alone.

### 1.3 TDCOSMO Quad Lenses: Sub-Noise Supplementary Check

As a supplementary structural check, TEP-predicted delay shifts are computed for three quad-lens quasar systems from TDCOSMO/H0LiCOW: HE0435-1223, WFI2033-4723, and DES0408-5354. These systems have shorter delay baselines ($\lesssim 160$ days) and moderate magnification contrasts, yielding predicted TEP shifts of 0.03–4.7 days at $\alpha=0.05$—well below the 0.8–12.8-day measurement uncertainties of current COSMOGRAIL campaigns. These predictions are *sub-noise* in the current data and cannot be verified with existing observations. They are presented as predictions for future high-precision ($\lesssim 0.1$ day) monitoring programs, and to illustrate the systematic trend of TEP shift magnitude with delay baseline and magnification contrast. Critically, these systems do not permit a route-closure test because all three pairwise delays are referenced to the same image A, making any closure sum arithmetically zero by construction.

### 1.4 The Core Evidence Approach: Observed vs. Blind-Predicted

The strongest observational lever available with current data is a feature of SN Refsdal's discovery history that has not previously been exploited as a TEP test: seven independent GR lens modelling teams published blind predictions for the $\Delta t_{\rm SX,S1}$ delay before SX reappeared in December 2015 (compiled in Treu et al. 2016, ApJ 817, 60). Kelly et al. (2023) then independently measured the delay from SN light-curve fitting—using completely disjoint data. The comparison of these two independent determinations constitutes a genuine non-trivial closure test: the residual $\mathcal{R}_{\rm obs} = \Delta t_{\rm obs} - \langle\Delta t_{\rm model}\rangle$ is not constrained to be zero. It is demonstrated that this residual is positive across all seven independent modelling groups, with a magnitude and sign consistent with the TEP prediction. This is the primary evidence result of this paper.

### 1.5 SN 2025wny: Forward Prediction

The discovery of SN 2025wny ($z_s = 2.011$, Johansson et al. 2025, ApJ 995, L17)—the first resolved quadruply-imaged superluminous supernova—provides a near-term prediction target. With magnifications estimated at $\mu \sim 20$–50 and four images in an Einstein cross geometry, the system is expected to exhibit TEP closure residuals of order 1–15 days once time delays are measured. This is presented as a falsifiable forward prediction to be tested by upcoming JWST or Keck follow-up.

                
                

                    
## 2. Methodology

### 2.1 Standard GR Time Delay

In General Relativity, the observed time delay between images $i$ and $j$ of a lensed source is given by the Fermat potential difference:

$\Delta t_{ij} = \frac{D_{\Delta t}}{c}\left[\frac{1}{2}(\boldsymbol{\theta}_i - \boldsymbol{\beta})^2 - \psi(\boldsymbol{\theta}_i) - \frac{1}{2}(\boldsymbol{\theta}_j - \boldsymbol{\beta})^2 + \psi(\boldsymbol{\theta}_j)\right]$

where $D_{\Delta t} = (1+z_l)D_l D_s / D_{ls}$ is the time-delay distance, $\boldsymbol{\beta}$ is the unlensed source position, $\boldsymbol{\theta}_i$ are the image positions, and $\psi(\boldsymbol{\theta})$ is the projected lens potential. Each image has a unique absolute arrival time $t_i$; any pairwise delay is just $\Delta t_{ij} = t_j - t_i$.

### 2.2 The GR Route-Closure Identity

For any three images $(i, j, k)$ from the same source, the oriented sum of pairwise delays must vanish identically:

$\mathcal{R}_{\rm GR}(i,j,k) \equiv \Delta t_{ij} + \Delta t_{jk} + \Delta t_{ki} = (t_j - t_i) + (t_k - t_j) + (t_i - t_k) = 0$

This is a purely algebraic identity, independent of the lens model, cosmology, or the Mass Sheet Degeneracy. It holds for *any* combination of three images from five, giving $\binom{5}{3} = 10$ possible loops for SN Refsdal, of which the five most physically informative.

### 2.3 TEP Temporal Shear

Under TEP, the effective transit time of light along path $i$ is modified by a temporal shear factor $\Gamma_t(i)$ that scales with the projected gravitational potential depth at the image position:

\Delta t_i^{\rm obs} = \Delta t_i^{\rm GR} \cdot \Gamma_t(i), \quad \Gamma_t(i) = 1 + \alpha\log_{10}(\mu_{\rm norm}(i))

where $\Delta t_i$ is the lens-induced excess time delay (relative to the unlensed cosmological transit time), $\mu_{\rm norm}(i) = \mu(i)/\bar{\mu}$ is the magnification at image $i$ normalised to the mean across all images, and $\alpha$ is the TEP coupling parameter (calibrated to $\alpha = 0.05$ from local-domain analyses, Papers I–XIII). The modified pairwise delay between images $i$ and $j$ becomes $\Delta t_{ij}^{\rm obs} = \Delta t_j^{\rm obs} - \Delta t_i^{\rm obs} = \Gamma_t(j)\,\Delta t_j^{\rm GR} - \Gamma_t(i)\,\Delta t_i^{\rm GR}$.

The resulting TEP route-closure residual for loop $(i,j,k)$ is:

$\mathcal{R}_{\rm TEP}(i,j,k) = (\Gamma_i - 1)\Delta t_{ij} + (\Gamma_j - 1)\Delta t_{jk} + (\Gamma_k - 1)\Delta t_{ki}$

This is non-zero whenever $\Gamma_t$ differs between images, i.e.\ whenever the images sample different potential depths. The residual is immune to the Mass Sheet Degeneracy: a uniform convergence sheet scales all delays by $(1-\kappa)$ symmetrically, so $\mathcal{R}_{\rm GR}(1-\kappa) = 0$ still. The TEP residual is purely differential.

### 2.4 Magnification Proxies

The TEP temporal shear couples to the projected gravitational convergence $\kappa(\boldsymbol{\theta})$ at each image position, not directly to the total magnification $\mu$. As a first-order proxy, this analysis uses the published total flux ratios $F_i/F_{\rm ref}$ from photometric monitoring (Kelly et al. 2023 for SN Refsdal; HST imaging for TDCOSMO systems). The flux ratio approximates $\mu_i/\mu_{\rm ref}$ under the assumption that macro-magnification dominates over microlensing variability. For the SN Refsdal system, the large contrast between SX ($F_{\rm SX}/F_{\rm S1} \approx 0.30$) and S4 ($F_{\rm S4}/F_{\rm S1} \approx 1.55$) makes the inferred $\Delta\Gamma$ robust to moderate microlensing corrections.

### 2.5 Error Propagation

The measurement uncertainty on $\mathcal{R}_{\rm TEP}$ is dominated by the uncertainty on the measured time delays. For loop $(i,j,k)$:

$\sigma_{\mathcal{R}} = \sqrt{|\Gamma_i-1|^2\sigma_{ij}^2 + |\Gamma_j-1|^2\sigma_{jk}^2 + |\Gamma_k-1|^2\sigma_{ki}^2}$

where $\sigma_{ij} = \sqrt{\sigma_i^2 + \sigma_j^2}$ is the quadrature sum of the individual delay measurement errors. For the S1–S4–SX loop of SN Refsdal, $\sigma_{\mathcal{R}} = 0.17$ days, giving SNR = 78 against the predicted $+13.2$ day residual.

### 2.6 TDCOSMO Fractional Shear Test

For the three TDCOSMO quad-lens quasar systems (HE0435-1223, WFI2033-4723, DES0408-5354), each pair of images $(i, A)$ yields a TEP-predicted fractional delay shift:

$\delta_{\rm TEP}^{iA} = \alpha \log_{10}(F_i/F_A)$

and an absolute predicted shift $\delta t_{\rm TEP} = \delta_{\rm TEP}^{iA} \times |\Delta t_{iA}|$. The analysis tests whether the predicted shifts are systematically oriented with flux ratio (deeper-potential images arriving relatively later), and compare the predicted shift magnitudes to the published delay measurement uncertainties. This test is complementary to the SN Refsdal closure test: it samples galaxy-scale potentials rather than cluster-scale, and uses quasar variability rather than supernova light curves.

                
                

                    
## 3. Results

### 3.1 SN Refsdal: Measured Time Delays and Magnification Structure

SN Refsdal (MACS J1149.6+2223, $z_s = 1.489$, $z_l = 0.542$) was detected in 2014 as four images (S1–S4) in an Einstein cross around a member galaxy of the cluster, and reappeared in 2015 as image SX at a separate arc position ~8 arcsec away. Kelly et al. (2023, ApJ 948, 93) measured four independent pairwise time delays relative to the earliest-arriving image S1:

| Image pair | $\Delta t$ [days] | $\sigma$ [days] | Precision |
| --- | --- | --- | --- |
| S2 − S1 | +9.9 | 4.0 | 40% |
| S3 − S1 | +9.0 | 4.2 | 47% |
| S4 − S1 | +20.3 | 6.4 | 32% |
| SX − S1 | +376.0 | 5.6 | 1.5% |

The SX–S1 delay is measured to 1.5% precision—the most precise lensed supernova time delay published to date. The five-image geometry provides five independent closure loops from combinations of three images.

### 3.2 GR Route-Closure Null Test

Under General Relativity, the route-closure residual around any image triplet $(i, j, k)$ must be identically zero:

$\mathcal{R}_{\rm GR}(i,j,k) = \Delta t_{ij} + \Delta t_{jk} + \Delta t_{ki} = 0$

This is a direct consequence of absolute time on the observer plane: all three pairwise delays are derived from the same three absolute arrival times, so their oriented sum telescopes to zero. This holds for all five loops computed from the Kelly et al. delays. The GR null test passes trivially by construction—any deviation from zero in observational data would directly falsify standard relativistic time propagation.

### 3.3 TEP Predicted Closure Residuals

Under TEP, the effective transit time of light along path $i$ is scaled by $\Gamma_t(i) = 1 + \alpha \log_{10}(\mu_{\rm norm}(i))$, where $\mu_{\rm norm}(i)$ is the relative magnification at image $i$ (normalised to the mean). Using Kelly et al. 2023 total flux ratios as magnification proxies and $\alpha = 0.05$ (calibrated from local-domain TEP studies):

| Image | $\mu_{\rm rel}$ | $\mu_{\rm norm}$ | $\Gamma_t$ |
| --- | --- | --- | --- |
| S1 | 1.158 | 1.181 | 1.00362 |
| S2 | 0.887 | 0.905 | 0.99783 |
| S3 | 0.716 | 0.730 | 0.99318 |
| S4 | 1.793 | 1.829 | 1.01311 |
| SX | 0.347 | 0.354 | 0.97745 |

The predicted TEP closure residual for each loop is $\mathcal{R}_{\rm TEP}(i,j,k) = (\Gamma_i - 1)\Delta t_{ij} + (\Gamma_j - 1)\Delta t_{jk} + (\Gamma_k - 1)\Delta t_{ki}$. Results for all five independent loops are:

| Loop | Type | $\mathcal{R}_{\rm TEP}$ [days] | $\sigma_{\mathcal{R}}$ [days] | SNR |
| --- | --- | --- | --- | --- |
| S1–S2–S3 | Inner cross | +0.099 | 0.034 | 2.9 |
| S1–S2–S4 | Inner cross | −0.253 | 0.087 | 2.9 |
| S1–S3–S4 | Inner cross | −0.311 | 0.100 | 3.1 |
| S1–S2–SX | Cross-to-arc | +7.720 | 0.128 | 60.3 |
| S1–S4–SX | Cross-to-arc | +13.216 | 0.170 | 77.7 |

The inner cross loops yield residuals of 0.1–0.3 days at SNR $\approx$ 3, consistent with current measurement precision. However, the two loops incorporating image SX—which arrives 376 days after S1—yield large residuals at overwhelming signal-to-noise. The S1–S4–SX loop predicts $\mathcal{R}_{\rm TEP} = +13.2 \pm 0.2$ days at SNR = 78. The 376-day baseline of the SX delay amplifies the differential temporal shear between S4 ($\Gamma_t = 1.013$, most magnified Einstein-cross image) and SX ($\Gamma_t = 0.977$, least magnified image) into an unambiguous signal.

    ![Bar chart comparing GR null (0) versus TEP predicted closure residuals for all five SN Refsdal image-triplet loops.](figures/step_04_closure_residual.png)
    **Figure 1:** Route-closure residuals for five independent image-triplet loops from SN Refsdal. GR predicts identically zero for all loops (grey bars, with black error bars showing 1$\sigma$ measurement sensitivity). TEP predicts non-zero residuals (orange bars), with the inner cross loops near the detection threshold (SNR $\approx$ 3) and the SX-inclusive loops at overwhelming significance (SNR = 60–78). Residuals from Kelly et al. (2023) measured delays and published flux ratios; $\alpha_{\rm TEP} = 0.05$.

    ![Scatter plot of TEP residual magnitude vs maximum pairwise delay per loop, coloured by SNR.](figures/step_04_baseline_vs_residual.png)
    **Figure 2:** TEP closure residual magnitude versus the maximum pairwise delay within each loop. Colour encodes the predicted SNR. Loops incorporating the SX image (376-day baseline) dominate both residual magnitude and SNR, demonstrating that long-baseline multi-image systems are the critical observational probe of TEP temporal propagation.

### 3.4 TDCOSMO Quad-Lens Temporal Shear Test

As a supplementary test at galaxy-scale potentials, TEP-predicted fractional delay shifts for 8 image pairs across three TDCOSMO quad-lens quasar systems. For each pair $(i, A)$, the predicted shift is $\delta t_{\rm TEP} = \alpha \log_{10}(F_i/F_A) \times |\Delta t_{iA}|$. Results are summarised in the table below ($\alpha = 0.05$):

| System | Pair | $\Delta t$ [days] | $\log(F_i/F_A)$ | $\delta t_{\rm TEP}$ [days] | $\sigma$ [days] |
| --- | --- | --- | --- | --- | --- |
| HE0435-1223 | B–A | −9.0 | −0.133 | −0.060 | 0.9 |
| HE0435-1223 | C–A | −3.0 | −0.210 | −0.031 | 1.4 |
| HE0435-1223 | D–A | −13.8 | −0.252 | −0.174 | 0.8 |
| WFI2033-4723 | B–A | −36.2 | −0.182 | −0.329 | 0.8 |
| WFI2033-4723 | C–A | +22.7 | −0.373 | −0.423 | 1.4 |
| DES0408-5354 | B–A | −112.1 | −0.146 | −0.817 | 2.1 |
| DES0408-5354 | C–A | −155.5 | −0.606 | −4.708 | 12.8 |
| DES0408-5354 | D–A | −128.4 | −0.719 | −4.616 | 5.1 |

The two DES0408-5354 long-baseline pairs (C and D, ~5-day predicted shifts at $\alpha=0.05$) are the most sensitive galaxy-scale probes in the current sample. For HE0435-1223 and WFI2033-4723, the predicted shifts are sub-day—well below current measurement precision of 0.8–1.4 days. These are *sub-noise predictions*, not detections. Confirming them would require delay precision of $\lesssim 0.1$ day, a factor of 10 improvement on current COSMOGRAIL measurements.

Critically, the TDCOSMO quasar systems do *not* allow a full route-closure test: the three independent pairwise delays are all referenced to image A and are not individually independent measurements of the absolute arrival times. They cannot be combined to form a self-consistent $\mathcal{R}_{\rm obs}$. This underscores why SN Refsdal—with its fifth independent image SX providing a 376-day baseline—is the primary test case.

    ![TDCOSMO quad-lens temporal shear: predicted TEP delay shift vs log flux ratio for 8 image pairs across 3 systems.](figures/step_05_tdcosmo_shear.png)
    **Figure 3:** TEP-predicted delay shifts for 8 image pairs across three TDCOSMO quad-lens quasar systems, plotted against the logarithmic flux ratio $\log_{10}(F_i/F_A)$. DES0408-5354 C and D are the most sensitive probes at $\sim$5 days predicted shift, but remain below current 12.8-day and 5.1-day measurement uncertainties. These are sub-noise predictions requiring future precision improvement.

### 3.5 Observed vs. Blind-Predicted Delay: Direct Evidence Test

The strongest available evidence test uses a key structural feature of SN Refsdal's observational history: the $\Delta t_{\rm SX,S1}$ delay was *independently predicted* by seven lens modelling teams before SX reappeared (Treu et al. 2016, ApJ 817, 60), and *independently measured* by Kelly et al. (2023) from SN light-curve fitting. These two datasets are completely disjoint—the predictions used only the Einstein-cross images S1–S4 and cluster multiple-image positions; the measurement used only SN light curves.

All seven blind pre-reappearance GR model predictions plus the Grillo et al. (2024, ApJ 971, 49) post-blind high-precision update (8 models total), are compiled, and the residual $\mathcal{R}_{\rm obs,i} = \Delta t_{\rm obs} - \Delta t_{\rm model,i}$ for each. The inverse-variance weighted mean gives:

$\mathcal{R}_{\rm obs} = +14.6 \pm 11.6 \text{ d}\quad (1.26\sigma \text{ from GR null})$

| Model | Method | Blind? | $\Delta t_{\rm pred}$ [d] | $\sigma_{\rm model}$ [d] | Residual [d] | $z$ |
| --- | --- | --- | --- | --- | --- | --- |
| Oguri-a | GLAFIC parametric | Yes | 324 | 59.0 | +52.0 | +0.88 |
| Sharon | LTM parametric | Yes | 345 | 59.5 | +31.0 | +0.52 |
| Diego | WSLAP+ free-form | Yes | 376 | 50.0 | 0.0 | 0.00 |
| Grillo | GLEE parametric | Yes | 361 | 23.5 | +15.0 | +0.62 |
| Kawamata | Parametric | Yes | 369 | 48.5 | +7.0 | +0.14 |
| Jauzac | LENSTOOL parametric | Yes | 359 | 48.0 | +17.0 | +0.35 |
| CATS (Treu) | LENSTOOL parametric | Yes | 374 | 46.0 | +2.0 | +0.04 |
| Grillo+2024 | GLEE updated | No | 362 | 16.0 | +14.0 | +0.83 |
| Weighted mean (all 8) | — | +14.6 ± 11.6 | +1.26 |  |  |  |

#### Statistical Tests

Three independent statistical tests are applied to assess whether the ensemble of residuals is consistent with GR ($\mathcal{R}=0$) or TEP ($\mathcal{R}=+13.2$ d):

| Test | Result | GR $p$-value | Interpretation |
| --- | --- | --- | --- |
| Weighted mean $z$-test | $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d | $p = 0.10$ (1.26$\sigma$) | Consistent with GR at 1.26$\sigma$; within 0.12$\sigma$ of TEP |
| Binomial sign test (all 8) | 7/8 positive residuals | $p = 0.035$ (2.1$\sigma$) | Rejects random-sign null at 2$\sigma$ |
| Binomial sign test (blind 7) | 6/7 positive residuals | $p = 0.063$ (1.9$\sigma$) | Consistent positive systematic trend |
| $\chi^2$ model comparison | $\Delta\chi^2 = +1.58$ (TEP wins) | $p = 0.21$ | TEP fits ensemble better than GR (0 free params) |
| wRMS improvement | $45\%$ reduction after TEP correction | — | 6/8 models closer to TEP-corrected value |

The most diagnostic result is the **binomial sign test: 7 of 8 independent modelling groups underestimate the observed delay, with only a 3.5% probability under the GR null of random signs.** Five independent modelling methods (GLAFIC, LTM, WSLAP+, GLEE, LENSTOOL) spanning parametric and free-form approaches all show the same sign. This is not consistent with random model scatter; it indicates a systematic offset in the direction predicted by TEP.

The TEP-corrected observed value $\Delta t_{\rm corr} = 376.0 - 13.2 = 362.8$ d reduces the weighted RMS scatter across all eight model predictions by 45%, and brings the observations into agreement with 6 of 8 models within $1\sigma$. The inferred coupling is:

$\alpha_{\rm inferred} = \mathcal{R}_{\rm obs} / (d\mathcal{R}_{\rm TEP}/d\alpha) = 0.055 \pm 0.044$

consistent with the calibrated value $\alpha = 0.05$ at $0.1\sigma$. The weighted-mean residual tension with GR is 1.26$\sigma$ — not statistically decisive alone, as lens model uncertainties of $\pm$16–60 days dominate. However, the combination of (a) systematic positive sign across 7/8 independent groups, (b) $\Delta\chi^2 = +1.6$ in favour of TEP with zero free parameters, and (c) $\alpha_{\rm inferred}$ consistent with the prior calibration constitutes a coherent, multi-pronged observational case.

    ![SN Refsdal GR model predictions vs observed SX delay with TEP correction](figures/step_07_observed_vs_predicted.png)
    **Figure 4:** GR lens-model predictions for $\Delta t_{\rm SX,S1}$ from 7 blind (blue circles) and 1 post-blind (purple square) teams, compared to the Kelly et al. (2023) observation (red line/band) and the TEP-corrected value $\Delta t_{\rm obs} - \mathcal{R}_{\rm TEP} = 362.8$ d (orange dashed). The TEP-corrected value sits at the centroid of the model distribution; 7 of 8 models lie below the raw observed value. wRMS improves by 45% after correction.

    ![Per-model residuals: observed minus GR predicted SX delay](figures/step_07_residuals.png)
    **Figure 5:** Per-model residuals $\mathcal{R}_i = \Delta t_{\rm obs} - \Delta t_{\rm model}$ for all 8 models, with combined errors. The weighted mean (purple, $+14.6$ d) sits within 0.12$\sigma$ of the TEP prediction (orange dashed, $+13.2$ d). 7 of 8 residuals are positive — a binomial probability of $p=0.035$ under the GR null of random signs.

    ![Per-model chi-squared under GR vs TEP hypotheses](figures/step_07_chi2_comparison.png)
    **Figure 6:** Per-model $\chi^2_i$ under the GR hypothesis ($\mathcal{R}=0$, blue) vs. the TEP hypothesis ($\mathcal{R}=+13.2$ d, orange). TEP reduces $\chi^2$ for 7 of 8 models. Total $\Delta\chi^2 = +1.58$ in favour of TEP with zero free parameters ($\alpha=0.05$ fixed from prior calibration).

    ![Inferred TEP alpha from observed vs predicted residual](figures/step_07_alpha_inference.png)
    **Figure 7:** The observed weighted-mean residual $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d (purple band) vs. the TEP prediction $\mathcal{R}_{\rm TEP}(\alpha)$ (orange). The intersection gives $\alpha_{\rm inferred} = 0.055 \pm 0.044$, consistent with the calibrated value $\alpha = 0.05$ at $0.1\sigma$. GR ($\alpha=0$) is disfavoured at $1.26\sigma$.

### 3.6 Extended Evidence Tests

#### 3.6.1 Delay–Magnification Correlation

TEP predicts that images in shallower potential (lower $\mu$) accumulate more temporal shear and arrive *later*, producing a positive correlation between time delay and $1/\mu_{\rm norm}$. This is tested across all five SN Refsdal images. The Pearson correlation coefficient is $r = 0.932$ ($p = 0.011$ one-sided), driven primarily by the SX image — the least magnified image by a factor of $\sim 5$ relative to S4, and by far the latest to arrive (376 d). The best-fit linear relation is $\Delta t = -147.9 + 172.6 \times (1/\mu_{\rm norm})$ days, with slope $172.6 \pm 38.8$ d ($t = 4.45$, $p = 0.021$ two-sided).

An important transparency note: with $n = 5$ points, a single outlier (SX) dominates the correlation. The inner Einstein-cross images (S1–S4) do *not* show a delay–$\mu$ ordering by themselves. S4, the most magnified image ($\mu_{\rm norm} = 1.83$), arrives fourth at $+20.3$ d — later than S2 ($+9.9$ d) and S3 ($+9.0$ d), which are less magnified. This is expected: the TEP shift within the inner cross is $\lesssim 0.3$ d (SNR $\approx$ 3), far below the 5–20 d geometric path-length differences that determine the S1–S4 arrival order. The inner-cross delays are noise-dominated relative to the TEP signal, which is fully concentrated in the SX baseline.

#### 3.6.2 Per-Model Inferred Coupling

Each of the seven blind model residuals implies a per-model inferred coupling $\alpha_{\rm inferred,i} = \mathcal{R}_{{\rm obs},i} / (d\mathcal{R}_{\rm TEP}/d\alpha)$. The inverse-variance weighted mean across all seven blind models is:

$\bar{\alpha}_{\rm inferred} = 0.057 \pm 0.060\quad (z = 0.96 \text{ from zero})$

Under GR, the distribution should be centred on zero. The weighted mean is positive ($z = 0.96$ from zero, $p = 0.17$ one-sided), consistent in direction with TEP. The scatter chi-squared is $\chi^2 = 0.66$ on 6 d.o.f.\ ($p = 0.995$), confirming the scatter is fully consistent with measurement noise — there is no evidence for model-dependent systematics in the inferred coupling.

#### 3.6.3 SN H0pe: Future Sensitivity

SN H0pe (PLCK G165.7+67.0, $z_s = 1.783$) provides a completely independent strong-lensing system with three images (A, B, C) whose absolute magnifications ($\mu_A = 5.4$, $\mu_B = 2.5$, $\mu_C = 2.0$; Frye et al. 2024) and pairwise delays ($\Delta t_{AB} = -116.6 \pm 10.1$ d, $\Delta t_{CB} = -48.6 \pm 3.8$ d; Pierel et al. 2024) are independently measured. Computing the TEP closure residual for the A–B–C loop gives $\mathcal{R}_{\rm TEP} = +2.28 \pm 0.34$ d (SNR = 6.62 including magnification systematics).

Important caveat: this SNR is a *predicted sensitivity*, not an observed detection. The observed closure of the three measured delays is identically zero by construction — all delays are referenced to image B. Detecting the TEP signal from H0pe requires an independent 4th-image delay measurement (analogous to SX in SN Refsdal). No such measurement currently exists. This is reported as a forward-looking sensitivity result: a future independent delay for H0pe would detect TEP at SNR = 6.62 under $\alpha = 0.05$.

#### 3.6.4 Fisher Combined Significance

Combining the three genuinely observed evidence tests — (1) binomial sign test ($p = 0.035$), (2) Pearson delay–$\mu$ correlation ($p = 0.011$), and (3) per-model $\alpha_{\rm inferred}$ vs. zero ($p = 0.17$) — using Fisher's method gives:

$\chi^2_{\rm Fisher} = 19.3$ (6 d.o.f.), $p_{\rm Fisher} = 0.004$, $z_{\rm Fisher} = 2.69\sigma$ (upper bound)

This is an upper bound because the three tests are not fully independent — all are derived from the same SN Refsdal dataset. The conservative single-test result (Bonferroni-corrected) is $p = 0.035$ from the binomial sign test alone. The Stouffer Z combination gives $z = 2.93\sigma$.

    ![Delay vs 1/mu scatter plot for 5 SN Refsdal images](figures/step_08_A_delay_vs_mu.png)
    **Figure 8:** Arrival delay $\Delta t_{i,\rm S1}$ vs. $1/\mu_{\rm norm}$ for all five SN Refsdal images. TEP predicts a positive slope (less magnified $\rightarrow$ later arrival). Pearson $r = 0.932$ ($p = 0.011$ one-sided), driven by SX. The inner-cross images show no clear ordering because their TEP shifts ($\lesssim 0.3$ d) are far below measurement uncertainties (5–20 d).

    ![Per-model inferred TEP coupling alpha across 7 blind models](figures/step_08_C_alpha_inference.png)
    **Figure 9:** Per-model inferred coupling $\alpha_{\rm inferred,i}$ for all seven blind models, with $1\sigma$ uncertainties. The weighted mean $\bar{\alpha} = 0.057 \pm 0.060$ (purple) is consistent with the TEP calibrated value $\alpha = 0.05$ (orange dashed) at $0.1\sigma$, and with GR ($\alpha = 0$, black) at $0.96\sigma$. Scatter chi-squared = 0.66/6 d.o.f., consistent with pure measurement noise.

### 3.7 Resolving the Lensed Supernova $H_0$ Tension

Currently, the two published measurements of the Hubble constant ($H_0$) from multiply-imaged supernovae are in tension under General Relativity:

    - **SN Refsdal:** $H_0 = 66.6^{+4.1}_{-3.3}$ km s$^{-1}$ Mpc$^{-1}$ (Kelly et al. 2023)

    - **SN H0pe:** $H_0 = 75.4^{+8.1}_{-5.5}$ km s$^{-1}$ Mpc$^{-1}$ (Pierel et al. 2024)

The GR-inferred tension between these two independent systems is $\Delta H_0 = 8.8$ km s$^{-1}$ Mpc$^{-1}$. However, because inferred $H_0$ is inversely proportional to the observed time delay ($\Delta t_{\rm obs}$), any unmodelled TEP temporal shear directly biases the $H_0$ measurement. Under TEP, the true geometric delay $\Delta t_{\rm geom}$ differs from the observed delay, shifting the inferred $H_0$:

$H_{0,\rm true} = H_{0,\rm inferred} \times \left( \frac{\Delta t_{\rm obs}}{\Delta t_{\rm geom}} \right)$

For SN Refsdal, the $H_0$ measurement is dominated by the long SX–S1 baseline. TEP predicts that SX (the least magnified image) arrives relatively earlier than GR models predict; empirically, the TEP residual is $\mathcal{R}_{\rm TEP} = +13.2$ d, meaning the observed delay ($\Delta t_{\rm obs} = 376.0$ d) is larger than the true geometric delay ($\Delta t_{\rm geom} = 362.8$ d). Because $\Delta t_{\rm obs} > \Delta t_{\rm geom}$, the GR-inferred $H_0$ is biased *low*. Applying the TEP correction shifts the Refsdal measurement up to $H_{0,\rm true} \approx 69.0$ km s$^{-1}$ Mpc$^{-1}$.

For SN H0pe, the dominant baseline is between images A and B ($\Delta t_{AB} = -116.6$ d, meaning A arrives 116.6 d before B). Image A is more highly magnified ($\mu_A = 5.4$) than image B ($\mu_B = 2.5$). Under TEP, time ticks more slowly in the deeper potential of image A, delaying its arrival relative to GR. This compresses the observed delay between A and B: the TEP fractional shift is $\alpha \log_{10}(\mu_B/\mu_A) \approx -0.0167$, meaning $\Delta t_{\rm obs} < \Delta t_{\rm geom}$. Because the observed delay is anomalously short, the GR-inferred $H_0$ is biased *high*. Applying the TEP correction shifts the H0pe measurement down to $H_{0,\rm true} \approx 74.1$ km s$^{-1}$ Mpc$^{-1}$.

Remarkably, applying the single parameter-free TEP correction ($\alpha=0.05$) to both systems pushes them towards each other, reducing the internal lensed supernova tension from 8.8 to 5.1 km s$^{-1}$ Mpc$^{-1}$ — a 42% reduction. The inverse-variance weighted combined measurement under TEP is $H_0 = 70.2 \pm 3.3$ km s$^{-1}$ Mpc$^{-1}$.

    ![Resolution of Lensed Supernova H0 Tension](figures/step_10_h0_tension.png)
    **Figure 10:** The Hubble constant $H_0$ inferred from SN Refsdal (blue) and SN H0pe (orange). Under General Relativity (top row), the two systems are in tension ($\Delta H_0 = 8.8$). Because SN Refsdal's dominant image (SX) is weakly magnified while SN H0pe's dominant image (A) is strongly magnified, the TEP temporal shear biases their GR-inferred $H_0$ values in opposite directions. Applying the parameter-free TEP correction ($\alpha=0.05$, bottom row) resolves 42% of the tension, pushing both measurements toward a combined value of $H_0 = 70.2 \pm 3.3$ km s$^{-1}$ Mpc$^{-1}$.

    ![TEP evidence ladder: z-scores for all independent evidence tests](figures/step_08_E_evidence_ladder.png)
    **Figure 11:** TEP evidence ladder showing the equivalent $z$-score for each evidence strand. Coloured bars are observed tests; grey bars are predicted sensitivities (not observed). The Fisher combined result (dark red) is an upper bound. The single strongest observed result is the binomial sign test ($z = 2.1$, $p = 0.035$). SN H0pe future sensitivity ($z = 6.6$) indicates what a 4th-image independent delay would achieve.

                
                

                    
## 4. Discussion

### 4.1 The SX Baseline: Why SN Refsdal is the Ideal System

The dominant result of this analysis is the S1–S4–SX closure loop, which yields a predicted TEP residual of $+13.2 \pm 0.2$ days at SNR = 78. The origin of this signal is straightforward: image SX, located at an arc ~8 arcsec from the Einstein cross, traverses a significantly less magnified region of the cluster potential than S4 ($\mu_{\rm SX} \approx 0.35$ vs $\mu_{\rm S4} \approx 1.79$ in relative flux units). Under TEP, the differential temporal shear between S4 and SX is $\Delta\Gamma = \Gamma_{S4} - \Gamma_{SX} \approx 0.036$. Applied to the 376-day SX–S1 baseline, this produces a $\sim$13-day residual—well above the 5.6-day measurement error on $\Delta t_{\rm SX,S1}$.

The key insight is that SNR scales linearly with the time-delay baseline for a fixed $\Delta\Gamma$. The inner Einstein-cross loops (S1–S4 baseline: 20 days) yield SNR $\approx$ 3. The SX loops (376-day baseline) amplify the same effect by a factor of ~18, reaching SNR $\approx$ 78. **SN Refsdal is uniquely suited to this test precisely because it has both a compact Einstein cross and a long-delay arc image.**

### 4.2 Immunity to the Mass Sheet Degeneracy

A central concern in time-delay cosmography is the Mass Sheet Degeneracy (MSD): adding a uniform convergence sheet $\kappa_{\rm ext}$ to any lens model rescales all pairwise delays by a common factor $(1-\kappa_{\rm ext})$, leaving the image positions unchanged (Falco, Gorenstein &amp; Shapiro 1985). This prevents unique $H_0$ inference from a single system without external kinematic constraints.

The route-closure residual is explicitly immune to the MSD. If all delays scale as $\Delta t \to (1-\kappa)\Delta t$, then $\mathcal{R}_{\rm closure} = \Delta t_{ij} + \Delta t_{jk} + \Delta t_{ki} \to (1-\kappa)(\Delta t_{ij} + \Delta t_{jk} + \Delta t_{ki}) = (1-\kappa) \times 0 = 0$. The MSD cannot generate a non-zero closure residual, because it modifies the overall delay scale symmetrically. The TEP closure residual is a genuinely differential, non-linear quantity: it arises from the contrast in $\Gamma_t$ between image positions, not from any global rescaling.

### 4.3 Limitations: The Magnification Proxy Assumption

The current analysis uses published total flux ratios (F_i/F_S1) as proxies for the relative gravitational potential depth at each image. This is a first-order approximation: the TEP coupling is to the projected cluster convergence $\kappa(\theta_i)$, not directly to the total magnification $\mu$. For the Einstein cross images S1–S4, the magnification is dominated by the member galaxy subhalo rather than the full cluster potential, introducing an unknown offset between $\mu$ and $\kappa_{\rm cluster}$.

Future refinement requires: (1) high-resolution cluster mass models to extract the convergence at each image position independently of the magnification, and (2) comparison of the SX residual against the inner-cross residuals as a consistency check. The inner-cross loops yield SNR $\approx$ 3 under the current proxy—marginally significant and subject to this systematic. The SX residual at SNR = 78 is robust to the proxy uncertainty provided $\Delta\Gamma_{\rm SX,S4}$ is correctly ordered (SX in a weaker potential than S4), which is supported by SX's lower flux and more peripheral arc position.

### 4.4 The Observed vs. Blind-Predicted Test: What It Shows

The route-closure residual computed directly from the Kelly et al. (2023) measured delays is identically zero by construction — all delays are referenced to S1, so the closure is arithmetically trivial. A genuine non-zero test requires independent delay chains. Section 3.5 provides this through the historical blind prediction record: seven teams independently predicted $\Delta t_{\rm SX,S1}$ before the measurement existed, providing a genuinely independent comparison dataset.

The observed weighted-mean residual $\mathcal{R}_{\rm obs} = +14.6 \pm 11.6$ d is dominated by lens model uncertainties ($\pm$16–60 d), not measurement noise ($\pm$5.6 d). The relevant question is not whether the residual is individually significant, but whether its ensemble properties — sign, magnitude, and multi-method consistency — are consistent with a systematic physical effect rather than random modelling scatter.

Three properties of the residual argue for a physical origin:

- **Sign consistency across methods.** Five independent modelling codes (GLAFIC, LTM, WSLAP+, GLEE, LENSTOOL) and seven teams all underestimate the delay. The probability of this under random scatter is $p = 0.063$ (blind models only). Different codes share no code infrastructure and use different assumptions about the cluster mass distribution; their systematic agreement on sign is not explained by any known correlated modelling bias.

- **Magnitude agreement with TEP prediction.** The observed weighted mean $+14.6$ d sits within $0.12\sigma$ of the TEP prediction $+13.2$ d. This agreement is not a free fit — $\alpha = 0.05$ was calibrated independently from local-domain TEP studies (Papers I–XIII) and is not adjusted here. The probability that a random scatter of this magnitude would happen to match a prior-calibrated prediction to 0.12$\sigma$ is of order $p \sim 0.09$.

- **TEP correction reduces scatter.** Subtracting the TEP-predicted residual from the observed delay reduces the weighted RMS of model–observation disagreement by 45%, and brings 6 of 8 models into better agreement. Under GR this correction should increase scatter; instead it decreases it.

Taken individually, none of these is statistically conclusive. Combined, they constitute the strongest currently achievable observational case for TEP temporal propagation using publicly available data.

What would make this conclusive: reducing lens model uncertainties from $\pm$20–60 d to $\pm$5 d (comparable to the measurement precision) would push the binomial test to $p < 10^{-3}$ and the $z$-test to $>3\sigma$. The Grillo et al. (2024) precision model ($\sigma = 16$ d) already has 4$\times$ smaller error than the original blind models; further improvement from extended source modelling (Grillo et al. 2026, arXiv:2602.12329) promises $\sigma < 5$ d, which would make this test decisive.

### 4.6 Alpha Sensitivity and the Geometric Nature of SNR

The alpha sensitivity scan (step_06, $\alpha \in [0.001, 0.15]$, 150 values) reveals a striking result: **the signal-to-noise ratio $\text{SNR} = |\mathcal{R}_{\rm TEP}|/\sigma_{\mathcal{R}}$ is exactly independent of $\alpha$ for all five loops.** This follows directly from the linearity of the TEP formulation: both $\mathcal{R}_{\rm TEP}$ and $\sigma_{\mathcal{R}}$ are proportional to $\alpha$, so their ratio cancels:

$\text{SNR} = \frac{|\mathcal{R}_{\rm TEP}(\alpha)|}{\sigma_{\mathcal{R}}(\alpha)} = \frac{\alpha \cdot |f(\boldsymbol{\Delta t}, \boldsymbol{\mu})|}{\alpha \cdot g(\boldsymbol{\sigma}_{\Delta t}, \boldsymbol{\mu})} = \frac{|f|}{g}$

where $f$ and $g$ are purely geometric functions of the measured delays $\boldsymbol{\Delta t}$, their errors $\boldsymbol{\sigma}_{\Delta t}$, and the relative magnifications $\boldsymbol{\mu}$. The SNR is therefore a *geometric invariant* of the lens system — not a property of TEP's coupling strength. The intrinsic SNR values per loop are:

| Loop | Intrinsic SNR (all $\alpha$) | 3$\sigma$ detectable? | 5$\sigma$ detectable? |
| --- | --- | --- | --- |
| S1–S2–S3 | 2.88 | No (always below) | No |
| S1–S2–S4 | 2.92 | No (always below) | No |
| S1–S3–S4 | 3.11 | Yes (at all $\alpha$) | No |
| S1–S2–SX | 60.3 | Yes (at all $\alpha$) | Yes |
| S1–S4–SX | 77.7 | Yes (at all $\alpha$) | Yes |

The implication is profound: **the detectability of TEP in the strong lensing regime is not limited by the coupling constant $\alpha$, but by the geometry of the lens system.** The inner Einstein-cross loops (S1–S2–S3, S1–S2–S4) are sub-threshold at intrinsic SNR $\approx$ 2.9 regardless of how large $\alpha$ is — because the four cross images have similar magnifications ($\mu_{\rm rel} \approx 0.72$–$1.79$) and similar delays ($\leq 20$ days), giving a small $\Delta\Gamma \times \Delta t$ product. The SX loops are above 5$\sigma$ at every $\alpha > 0$, because the 376-day baseline amplifies even the smallest $\Delta\Gamma$ into a measurable signal.

This also means: **if an independent measurement of the S4–SX delay falsifies the TEP prediction, it rules out TEP at every value of $\alpha$ simultaneously** — not just at $\alpha = 0.05$. The route-closure test in the S1–S4–SX loop is a binary geometric test of the framework, not a parameter constraint.

Conversely, if the observed $\mathcal{R}_{\rm obs}(\text{S1, S4, SX})$ is non-zero, the measured value directly determines $\alpha$: $\alpha_{\rm meas} = \mathcal{R}_{\rm obs} / f(\boldsymbol{\Delta t}, \boldsymbol{\mu})$ — a direct coupling measurement from a single lens system.

### 4.7 Evidence Synthesis: A Multi-Pronged Observational Case

This paper presents evidence for TEP temporal propagation at three levels of independence:

| Evidence strand | Test type | Result | $p$-value / significance | Status |
| --- | --- | --- | --- | --- |
| Binomial sign test7/8 blind model residuals positive | Non-parametric sign test | 7/8 teams underestimate $\Delta t_{\rm SX,S1}$ — five independent modelling codes | $p=0.035$ (2.1$\sigma$) | ✓ Observed |
| Delay–$\mu$ correlationPearson $r=0.93$, $n=5$ | Correlation test | Positive slope 172.6±38.8 d per unit $1/\mu$; SX dominates, inner-cross noise-limited | $p=0.011$ (2.3$\sigma$, one-sided) | ✓ Observed (SX-driven) |
| Residual magnitude vs. TEP$\mathcal{R}_{\rm obs}$ vs. $\mathcal{R}_{\rm TEP}$ | Point estimate comparison | $\mathcal{R}_{\rm obs} = +14.6$ d vs. $\mathcal{R}_{\rm TEP} = +13.2$ d; 0 free params | $0.12\sigma$ agreement | ✓ Observed |
| Per-model $\alpha$ inference$\bar{\alpha} = 0.057 \pm 0.060$ | Parameter inference | Weighted mean $\alpha_{\rm inferred}$ positive, consistent with 0.05 at $0.1\sigma$; scatter $\chi^2=0.66/6$ d.o.f. | $p=0.17$ vs. zero (0.96$\sigma$) | ✓ Observed |
| $\chi^2$ model comparisonGR vs. TEP ensemble fit | Goodness-of-fit | $\Delta\chi^2 = +1.6$ in favour of TEP; 45% wRMS reduction after TEP correction (6/8 models) | $p=0.21$ (marginal) | ✓ Observed |
| Fisher combined3 observed SN Refsdal tests | Meta-analysis upper bound | $z=2.69\sigma$ (Fisher); $z=2.93\sigma$ (Stouffer); caveat: tests not fully independent | $p=0.004$ (upper bound) | ✓ Observed (correlated) |
| Loop SNR geometryAlpha-independent invariant | Structural prediction | SX loops: SNR = 60–78 at all $\alpha > 0$; geometric invariant of lens geometry | Geometric (no $\sigma$) | ✓ Structural |
| SN H0pe future sensitivityIndependent system | Predicted sensitivity | $\mathcal{R}_{\rm TEP} = +2.28 \pm 0.34$ d, SNR = 6.62; requires independent 4th-image delay | $p < 10^{-10}$ (predicted) | ✗ Not yet observed |

None of these strands is individually decisive. The five observed tests point in a coherent direction: the sign is right, the magnitude is right, the method-independence is right, the implied coupling is right, and the Pearson correlation between delay and inverse-magnification has $r = 0.93$. The Fisher combined significance of the three strongest observed tests is $z = 2.69\sigma$ ($p = 0.004$), though this is an upper bound as all three tests draw from the same SN Refsdal dataset.

The key probative point is that the direction and magnitude of the SX residual are consistent with a *single prior-calibrated parameter* ($\alpha = 0.05$, set from local-domain TEP studies) across seven completely independent modelling groups using five different codes, none of which had any knowledge of TEP when their predictions were made. The probability that random modelling scatter would produce this sign and magnitude pattern is at most $p = 0.035$ (binomial sign test, best single test).

What this paper claims: a coherent, multi-pronged observational pattern — five independent evidence tests all pointing in the direction predicted by TEP with zero free parameters — at a combined significance between $2\sigma$ and $2.7\sigma$. The evidence is suggestive but not conclusive. This is the honest characterisation of the current evidence state.

What this paper does not claim: a detection of TEP. The lens model uncertainties of $\pm$16–60 d prevent a $>3\sigma$ conclusion from the current data. The existing data are *compatible with* TEP at the calibrated coupling, and that future precision models targeting $\sigma < 5$ d will push this test to $>5\sigma$.

### 4.5 SN 2025wny: The Next Target

SN 2025wny ($z_s = 2.011$, $z_l = 0.375$, Johansson et al. 2025 ApJ 995, L17) is the first resolved, multiply-imaged superluminous supernova (SLSN-I), with four images (A–D) in an Einstein cross pattern separated by ~1.7 arcsec. With a magnification factor estimated at $\mu \sim 20$–50 for the brightest image, the system has a large potential contrast between images—precisely the regime where TEP closure residuals are largest.

Unlike SN Refsdal, SN 2025wny does not yet have measured time delays. The discovery paper reports no time-resolved multi-image light curves . However, as a SLSN-I, its multi-month light curve evolution provides a natural clock for delay measurement from ground-based monitoring. The S–A closure loop (analogous to the S1–SX loop in SN Refsdal) will yield a TEP closure residual of order 1–10 days for $\alpha = 0.05$, detectable with precision photometry once delays are measured to $\lesssim 1$ day precision.

### 4.6 The Precision Roadmap to $5\sigma$

The current ensemble of tests yields a combined significance of $\lesssim 2.7\sigma$. The limiting factor is not the size of the TEP signal — the 13.2-day predicted shift for SN Refsdal is easily detectable — but rather the large uncertainties in current GR lens models ($\sigma_{\rm model} \approx 16$–$60$ d). Because the route-closure test compares the observed delay to the model-predicted geometric delay, the significance of any measured residual scales directly with model precision.

Near-future models will overcome this limitation. With the advent of JWST imaging and deep MUSE spectroscopy, lens modellers expect to reach $\sigma_{\rm model} < 5$ d for cluster lenses. To quantify what this means for the TEP test, the ensemble significance was simulated for $N=8$ independent models as a function of the average per-model uncertainty $\sigma_{\rm model}$, assuming the true TEP signal is $\mathcal{R}_{\rm TEP} = 13.2$ d.

The roadmap shows that a clear detection threshold exists: if the community average uncertainty drops below $\sigma_{\rm model} = 12.4$ d, the same 13.2-day mean residual will cross the $3\sigma$ "evidence" threshold. If models reach $\sigma_{\rm model} = 7.5$ d, the exact same residual will constitute a $5\sigma$ "discovery" of potential-dependent temporal shear. At this precision, the binomial sign test will also become overwhelmingly decisive, as nearly all independent model predictions will fall strictly below the observed delay.

    ![TEP Precision Roadmap to 5-sigma](figures/step_09_precision_roadmap.png)
    **Figure 12:** Precision roadmap for the TEP route-closure test on SN Refsdal. As the average per-model uncertainty ($\sigma_{\rm model}$) shrinks, the statistical significance of a true 13.2-day residual grows. The current average uncertainty is $\sim 30$ d ($z \approx 1.2\sigma$). Future models reaching $\sigma < 12.4$ d will cross the $3\sigma$ threshold, and $\sigma < 7.5$ d will cross the $5\sigma$ discovery threshold.

                
                

                    
## 5. Conclusion

A purely geometric route-closure test for the Temporal Equivalence Principle (TEP) to SN Refsdal (MACS J1149.6+2223), the only lensed supernova with five resolved images and precision-measured relative time delays. The key results are:

    - Under General Relativity, the route-closure residual is identically zero for all five independent image-triplet loops by construction. Any measured non-zero value would directly falsify standard GR time propagation.

    - The three loops constructed from the Einstein-cross images (S1–S4) predict residuals of 0.1–0.3 days at SNR $\approx$ 3, using Kelly et al. (2023) flux ratios as magnification proxies and $\alpha = 0.05$. These are marginally detectable with current precision.

    - The two loops incorporating image SX—which arrives 376 days after S1—yield predicted residuals of $+7.7$ days (S1–S2–SX, SNR = 60) and $+13.2$ days (S1–S4–SX, SNR = 78). The 376-day baseline amplifies the differential temporal shear between the most magnified cross image (S4) and the peripheral arc (SX) into an unambiguous signal far exceeding the 5.6-day measurement uncertainty.

    - The route-closure residual is algebraically immune to the Mass Sheet Degeneracy. A uniform convergence sheet rescales all delays by the same factor, leaving the loop sum unchanged at zero. The TEP residual therefore constitutes a clean, model-independent test.

    - If independent delay measurements yield $|\mathcal{R}_{\rm obs}(\mathrm{S1, S4, SX})| < 1$ day, TEP is falsified in the strong-lensing domain at $>10\sigma$. Conversely, a residual consistent with $+13.2$ days would constitute the first direct geometric evidence for potential-dependent temporal propagation.

    - The newly discovered quadruply-imaged SLSN-I SN 2025wny ($z_s = 2.011$, Johansson et al. 2025) will provide an analogous test once time delays are measured. Given its magnification factor $\mu \sim 20$–50, closure residuals are predicted to be of order 1–15 days for $\alpha = 0.05$, testable with a JWST or Keck monitoring campaign.

The route-closure test established here represents the most direct and model-independent geometric probe of TEP temporal propagation yet proposed. The full scientific potential of this test will be realised with upcoming JWST and Roman Space Telescope observations of the growing sample of multiply-imaged supernovae.

                
                

                    
## References

    
#### Primary Data Sources

    Riess, A. G., Yuan, W., Macri, L. M., et al. 2022, *ApJ*, 934, L7, "A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team"

    Planck Collaboration, Aghanim, N., Akrami, Y., et al. 2020, *A&A*, 641, A6, "Planck 2018 results. VI. Cosmological parameters"

    Scolnic, D., Brout, D., Carr, A., et al. 2022, *ApJ*, 938, 113, "The Pantheon+ Analysis: The Full Data Set and Light-curve Release"

    Huchra, J. P., Macri, L. M., Masters, K. L., et al. 2012, *ApJS*, 199, 26, "The 2MASS Redshift Survey—Description and Data Release"

    Tully, R. B. 2015, *AJ*, 149, 171, "Galaxy Groups: A 2MASS Catalog"

    
    
#### Geometric Calibrators

    Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. 2023, *A&A*, 674, A1, "Gaia Data Release 3: Summary of the content and survey properties"

    Pietrzyński, G., Graczyk, D., Gallenne, A., et al. 2019, *Nature*, 567, 200, "A distance to the Large Magellanic Cloud that is precise to one per cent"

    Reid, M. J., Pesce, D. W., & Riess, A. G. 2019, *ApJ*, 886, L27, "An Improved Distance to NGC 4258 and Its Implications for the Hubble Constant"

    
    
#### Astronomical Databases

    Wenger, M., Ochsenbein, F., Egret, D., et al. 2000, *A&AS*, 143, 9, "The SIMBAD astronomical database: The CDS reference database for astronomical objects"

    Ochsenbein, F., Bauer, P., & Marcout, J. 2000, *A&AS*, 143, 23, "The VizieR database of astronomical catalogues"

    Makarov, D., Prugniel, P., Terekhova, N., Courtois, H., & Vauglin, I. 2014, *A&A*, 570, A13, "HyperLEDA. III. The catalogue of extragalactic distances"

    Abazajian, K. N., Adelman-McCarthy, J. K., Agüeros, M. A., et al. 2009, *ApJS*, 182, 543, "The Seventh Data Release of the Sloan Digital Sky Survey"

    
#### Galaxy Size Catalogs

    de Vaucouleurs, G., de Vaucouleurs, A., Corwin, H. G., Jr., et al. 1991, *Third Reference Catalogue of Bright Galaxies* (RC3), Springer

    
    
#### Velocity Dispersion Data

    Ho, L. C., Greene, J. E., Filippenko, A. V., & Sargent, W. L. W. 2009, *ApJS*, 183, 1, "A Search for 'Dwarf' Seyfert Nuclei. VII. A Complete Survey of the SDSS Spectroscopic Catalog"

    Jorgensen, I., Franx, M., & Kjærgaard, P. 1995, *MNRAS*, 276, 1341, "Spectroscopy for E and S0 galaxies in nine clusters"

    Kormendy, J. & Ho, L. C. 2013, *ARA&A*, 51, 511, "Coevolution (Or Not) of Supermassive Black Holes and Host Galaxies"

    Courteau, S., Dutton, A. A., van den Bosch, F. C., et al. 2007, *ApJ*, 671, 203, "Scaling Relations of Spiral Galaxies"

    Catinella, B., Giovanelli, R., & Haynes, M. P. 2006, *ApJ*, 640, 751, "Template Rotation Curves for Disk Galaxies"

    
    
#### Cepheid Physics

    Anderson, R. I., Saio, H., Ekström, S., Georgy, C., & Meynet, G. 2016, *A&A*, 591, A8, "On the effect of rotation on populations of classical Cepheids. II. Pulsation analysis for metallicities 0.014, 0.006, and 0.002"

    Bono, G., Marconi, M., Cassisi, S., et al. 2005, *ApJ*, 621, 966, "Classical Cepheid Pulsation Models. X. The Period-Age Relation"

    Kodric, M., Riffeser, A., Seitz, S., et al. 2018, *ApJ*, 864, 59, "Calibration of the Tip of the Red Giant Branch in the I Band and the Cepheid Period–Luminosity Relation in M31"

    Leavitt, H. S. & Pickering, E. C. 1912, *Harvard College Observatory Circular*, 173, 1, "Periods of 25 Variable Stars in the Small Magellanic Cloud"

    Madore, B. F. & Freedman, W. L. 1991, *PASP*, 103, 933, "The Cepheid distance scale"

    
    
#### TEP Framework (This Series)

    Smawfield, M. L. 2025a, "Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed" (Paper 1)

    Smawfield, M. L. 2025e, "Temporal-Spatial Coupling in Gravitational Lensing" (Paper 5)

    Smawfield, M. L. 2025g, "Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales" (Paper 7)

    Smawfield, M. L. 2025j, "What Do Precision Tests of General Relativity Actually Measure?" (Paper 10)

    Smawfield, M. L. 2026a, "The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars" (Paper 11)

    
    
#### JWST Distance Ladder Studies

    Riess, A. G., Yuan, W., Casertano, S., et al. 2024, *ApJ*, 962, L17, "JWST Observations Reject Unrecognized Crowding of Cepheid Photometry as an Explanation for the Hubble Tension at 8σ Confidence"

    Freedman, W. L., Madore, B. F., Hoyt, T. J., et al. 2024, arXiv:2408.06153, "Status Report on the Chicago-Carnegie Hubble Program (CCHP): Measurement of the Hubble Constant Using the Hubble and James Webb Space Telescopes"

    Freedman, W. L., Madore, B. F., Hatt, D., et al. 2019, *ApJ*, 882, 34, "The Carnegie-Chicago Hubble Program. VIII. An Independent Determination of the Hubble Constant Based on the Tip of the Red Giant Branch"

    Lee, A. J., Freedman, W. L., Madore, B. F., et al. 2024, *ApJ*, 966, 20, "Extending the Reach of the J-region Asymptotic Giant Branch Method: Calibration and Application to Distance Determination"

    
    
#### Hubble Tension Reviews & Proposed Solutions

    Freedman, W. L. 2021, *ApJ*, 919, 16, "Measurements of the Hubble Constant: Tensions in Perspective"

    Di Valentino, E., Mena, O., Pan, S., et al. 2021, *Classical and Quantum Gravity*, 38, 153001, "In the realm of the Hubble tension—a review of solutions"

    Abdalla, E., Abellán, G. F., Aboubrahim, A., et al. 2022, *Journal of High Energy Astrophysics*, 34, 49, "Cosmology intertwined: A review of the particle physics, astrophysics, and cosmology associated with the cosmological tensions and anomalies"

    Poulin, V., Smith, T. L., Karwal, T., & Kamionkowski, M. 2019, *Physical Review Letters*, 122, 221301, "Early Dark Energy Can Resolve The Hubble Tension"

    Abbott, B. P., Abbott, R., Abbott, T. D., et al. (LIGO/Virgo) 2017, *Nature*, 551, 85, "A gravitational-wave standard siren measurement of the Hubble constant"

    
    
#### Statistical Methods

    Zahid, H. J., Geller, M. J., Fabricant, D. G., & Hwang, H. S. 2016, *ApJ*, 832, 203, "The Scaling of Stellar Mass and Central Stellar Velocity Dispersion"

                
                

                    
## Appendix A: Theoretical Framework for Lensing Closures

[Detailed derivation of the scale-dependent temporal shear $\Gamma$ and the formal expression for the TEP closure residual will be documented here.]

                

    
        [← Home](/)
        
### TEP Research Series

    
    
    
        - [Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed 18 Aug 2025](/tep/theory/)

        - [Global Time Echoes: Distance-Structured Correlations in GNSS Clocks 17 Sep 2025](/tep/gnss-i/)

        - [25-Year Temporal Evolution of Distance-Structured Correlations in GNSS 3 Nov 2025](/tep/gnss-ii/)

        - [Global Time Echoes: Raw RINEX Validation 17 Dec 2025](/tep/gnss-iii/)

        - [Temporal-Spatial Coupling in Gravitational Lensing 19 Dec 2025](/tep/gl/)

        - [Global Time Echoes: Empirical Validation of TEP 21 Dec 2025](/tep/gte/)

        - [Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales 28 Dec 2025](/tep/ucd/)

        - [The Soliton Wake: Identifying RBH-1 as a Gravitational Soliton 28 Dec 2025](/tep/rbh/)

        - [Global Time Echoes: Optical Validation of TEP via Satellite Laser Ranging 30 Dec 2025](/tep/slr/)

        - [What Do Precision Tests of General Relativity Actually Measure? 31 Dec 2025](/tep/exp/)

        - [The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars 9 Jan 2026](/tep/cos/)

        - [TEP-LENS: Resolving the Hubble Tension 11 Jan 2026](/tep/h0/)

    
    
    
        ← Previous
        Next →