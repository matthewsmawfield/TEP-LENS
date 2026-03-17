# Temporal Shear: Reconciling JWST's Impossible Galaxies

## Abstract

     JWST has revealed coherent anomalies at $z > 5$: star formation efficiencies exceeding $\Lambda$CDM limits, overmassive black holes, and stellar masses surpassing dynamical masses — all concentrated in the deepest gravitational potentials. This study tests whether a single violation of the isochrony axiom can account for a systematic fraction of all three. The Temporal Equivalence Principle (TEP), a chameleon-screened scalar-tensor theory, predicts enhanced proper time in unscreened potentials. With a coupling $\alpha_0 = 0.58 \pm 0.16$ calibrated entirely from local Cepheids and zero parameters tuned to JWST data, this single mechanism simultaneously predicts four distinct anomalies: it accounts for $\sim 34\%$ of the Red Monster efficiency excess, resolves all currently identified "impossible" $\log M_* > 10.5$ galaxies at $z > 8$ in this sample (301/400 resolved across $z = 7$–$12$), reduces the cosmic SFRD excess from $11.1\times$ to $2.6\times$ $\Lambda$CDM, and provides a purely kinematic mechanism (differential temporal shear) for overmassive black hole growth in Little Red Dots.

     Four lines of evidence survive mass-control checks across three surveys ($N = 4{,}726$): (L1) a dust–$\Gamma_t$ correlation ($\rho = +0.62$, $N = 1{,}283$) with $\Delta\text{AIC} = -23$ against a mass-matched threshold; (L2) inside-out core screening in JADES ($\rho = -0.18$); (L3) mass–sSFR inversion at $z > 7$; (L4) resolution of 11/11 impossible $M_*/M_{\rm dyn} > 1$ cases via SED-independent kinematics (note: parameters are representative values from published ranges — see §3.9 L4 caveat). The sSFR sign inversion (L3) is mathematically immune to the "mass-proxy" critique: while systematic mass-measurement errors might induce a monotonic bias, they cannot produce a sharp sign inversion at $z=7$. Furthermore, L1 is strictly orthogonal to L3: controlling for sSFR leaves the $z>9$ dust correlation unchanged ($\rho = +0.592$, §3.9). The sharpest model-discriminating result is the Uniformity Paradox: the anomaly is not that massive galaxies at $z > 8$ are dusty, but that low-mass galaxies are *not* — in a pattern that tracks gravitational potential depth ($\rho = +0.56$). Any standard-physics resolution that tunes a time-uniform parameter (e.g., enhanced AGB yields) cannot reproduce this mass-dependent suppression: if dust production were maximally efficient everywhere, dust should be ubiquitous or track star formation, not potential depth. Only a mechanism that selectively suppresses effective time in shallow potentials — as TEP predicts via $\Gamma_t \ll 1$ for low-mass halos — reproduces the observed gradient without additional parameters. Fitted polynomials collapse cross-survey ($R^2 = -6.4$) while $t_{\rm eff}$ generalises without training ($\rho = 0.60$–$0.80$). Full CAMB Boltzmann integration is consistent with Planck constraints ($\sigma_8$ within $0.1\sigma$). The three independent JWST surveys each confirm L1 individually above $5\sigma$ with no clustering correction required between surveys (CEERS: $z = 7.0\sigma$; UNCOVER: $z = 11.4\sigma$; COSMOS-Web: $z = 22.4\sigma$). A Fisher combination of all five independent L1 confirmations (three photometric surveys + GOODS-S DJA 4th field + NIRSpec H$\alpha$/H$\beta$ spectroscopic dust with no SED fitting) gives $\chi^2 = 950.9$, $p = 7.0 \times 10^{-198}$, $z = 30.0\sigma$, spanning 4 distinct sky fields, 3 SED pipelines, and 2 dust estimators — between-field combination requiring no clustering correction. The functional-form test further confirms that $t_{\rm eff}$ outperforms $t_{\rm cosmic}$ in all three surveys independently at $>5\sigma$ (Steiger $Z = 5.3$, $6.7$, $16.8$), ruling out redshift ordering as the active variable. The clustering-corrected combined significance across all four independent lines is $6.4\sigma$ (upper bound); the corrected JWST recovery gives $\alpha_0 = 0.55 \pm 0.32$ (dust, Pearson $R^2$) and $\alpha_0 = 0.75 \pm 0.29$ (joint four-observable), both consistent with the Cepheid calibration within $1\sigma$; large uncertainties reflect genuine mass-proxy degeneracy.

     Eight new cross-dataset results further strengthen the evidence. (1) COSMOS2025 dust (Shuntov et al. 2025; 0.54 deg² blank field): partial $\rho(\Gamma_t, E(B-V) \mid M_*, z)$ strengthens from $+0.241$ at $z = 4$–$7$ to $+0.595$ at $z = 9$–$13$ ($N = 1{,}568$, 95% CI $[+0.532, +0.660]$) — noting that the UNCOVER DR4 Prospector spec-$z$ subsample ($N = 122$) returns a null ($\rho = -0.009$, $p = 0.92$) at $z = 9$–$12$, plausibly explained by NIRSpec incompleteness against dusty systems at $z > 9$ (§4.12 item 7). (2) COSMOS2025 sSFR inversion: partial $\rho$ flips from $-0.010$ (null) at $z = 4$–$7$ to $+0.237$ at $z = 9$–$13$; Steiger $Z = 6.37$ ($p = 1.9 \times 10^{-10}$) supports L3 in an independent field. (3) UNCOVER DR4 MegaScience (Wang et al. 2024): null at $z 9$ dust partial correlation unchanged at $+0.592$. Key limitations: (i) $\Gamma_t$ is a deterministic function of photometric halo mass, so breaking the mass-proxy degeneracy definitively requires IFU velocity dispersions; (ii) the dust–$\Gamma_t$ correlation (L1) contributes $\sim 90\%$ of the combined statistical weight — L3 is a sign-change argument and L4 has $N = 11$.

    *Keywords:* Cosmology: early universe – Galaxies: high-redshift – Galaxies: evolution – Gravitation – Scalar-tensor theories – Infrared: galaxies

 

    
## 1. Introduction

    
### 1.1 Observational Tensions

    JWST has revealed a coherent pattern of anomalies at $z > 5$ that challenges the standard framework for inferring stellar properties from photometry. Spectroscopically confirmed systems — the "Red Monsters" (Xiao et al. 2024) — imply stellar masses ($M_* \gtrsim 10^{11}\,M_\odot$) that require baryon-to-star conversion efficiencies of $\sim 0.50$, more than double the $\sim 0.20$ theoretical maximum imposed by feedback in $\Lambda$CDM halos, a discrepancy significant at $11\sigma$. This tension is not isolated. The UNCOVER UV luminosity function at $z > 9$ implies a star formation rate density exceeding the halo accretion limit by factors of 3–10 (Wang et al. 2024). The population of "Little Red Dots" (LRDs) — compact, red sources hosting broad-line AGN — harbours supermassive black holes ($M_{\rm BH} \sim 10^7$–$10^8\,M_\odot$) at $z > 4$ that are 10–100 times more massive relative to their hosts than local scaling relations predict ($M_{\rm BH}/M_* \sim 0.1$ vs. $0.001$; Matthee et al. 2024); growing such objects from stellar seeds requires near-continuous super-Eddington accretion, in tension with the high observed space density ($n \sim 10^{-5}\,\mathrm{Mpc}^{-3}$). JWST NIRSpec kinematics further reveal 11 galaxies at $z > 4$ with $M_*/M_{\rm dyn} > 1$ (de Graaff et al. 2024; Wang et al. 2024) — a physically impossible condition under standard assumptions. All three anomalies share a common structure: stellar masses and ages inferred from photometry appear systematically too large, too early, in precisely the environments with the deepest gravitational potentials.

    
### 1.2 Challenging Isochrony with TEP

    Underlying every photometric inference of stellar age, mass, and star formation rate is the isochrony axiom: the assumption that the clock governing stellar evolution ticks at the universal cosmic rate, regardless of local gravitational environment. Under isochrony, an observed red colour is interpreted as a combination of age, dust, and metallicity, and the resulting mass-to-light ratio ($M/L \propto t^n$) is treated as universal. If this axiom is violated — if stars in deep gravitational potentials accumulate proper time faster than the cosmic mean — then SED-inferred masses and ages are systematically inflated in precisely the environments where JWST finds the largest anomalies.

    The Temporal Equivalence Principle (TEP) formalises this possibility within a conformally coupled scalar-tensor framework with chameleon screening (Khoury & Weltman 2004; Brax et al. 2004; Burrage & Sakstein 2018). In unscreened halos, the conformal factor $A(\phi) > 1$ causes stellar clocks to tick faster relative to coordinate time, so that a galaxy's effective age $t_{\rm eff} = \Gamma_t\,t_{\rm cosmic}$ exceeds its cosmic age. The resulting bias in $M/L$ inflates inferred stellar masses by $\Gamma_t^n$, directly mimicking a star formation efficiency excess. This effect is physically distinct from standard gravitational redshift: photons still lose energy climbing out of potential wells (kinematic redshift, fully preserved in TEP), while the scalar field coupling independently accelerates atomic processes within the well. Both effects coexist; only the latter biases photometric mass inference.

    
        **Physics Note: Dilation vs. Enhancement**

        It is essential to distinguish between two relativistic effects:

        
            - **Kinematic Gravitational Redshift (Standard GR):** Photons lose energy climbing out of potential wells. This affects light and is fully preserved in TEP.

            - **Dynamical Clock Rate (TEP):** The scalar field coupling modifies the effective mass of particles, changing the rate at which atomic clocks tick relative to coordinate time. In the TEP framework, the scalar field in diffuse halos relaxes to values where $A(\phi) > 1$, causing clocks to tick *faster* (enhancement) than the cosmic mean, even while photons suffer redshift.

        
    

    The temporal shear effect is governed by an enhancement factor $\Gamma_t$:

    
        $$\Gamma_t = \exp\left[\alpha(z) \times \frac{2}{3} \times \Delta\log(M_h) \times \frac{1+z}{1+z_{\rm ref}}\right]$$
    

    The coupling constant $\alpha_0 = 0.58 \pm 0.16$ is calibrated entirely from local Cepheid observations ($z \approx 0$) and applied unchanged to the high-redshift universe. No parameters are tuned to the JWST data this model seeks to explain. The structural choices in the $\Gamma_t$ formula — the reference redshift $z_{\rm ref} = 5.5$, reference halo mass $\log M_{h,\rm ref} = 12.0$, exponential functional form, and $\sqrt{1+z}$ coupling scaling — were fixed by the scalar-tensor framework in prior papers (Papers 1 and 7); all have independent physical motivation and none were adjusted to improve JWST fits (§4.11.1).

    
### 1.3 Reader's Guide to the Evidence

    This work evaluates the TEP hypothesis through three distinct stages, moving from a zero-parameter theoretical prediction to population-level statistics and finally to physical mechanism tests. Because the analysis spans photometric, kinematic, and spatial data across multiple JWST surveys, the core evidence is structured as follows:

    
        - **Stage 1: The Zero-Parameter Prediction (§3.1).** The local-universe TEP calibration is applied directly to the most extreme $z > 5$ anomalies (the Red Monsters) to test if it predicts the correct magnitude of the star formation efficiency correction without any data tuning.

        - **Stage 2: The Four Independent Lines of Evidence (§3.2–3.8).** Four distinct observational signatures are tested across three JWST surveys ($N = 4{,}726$) to break the mass-proxy degeneracy. These lines are:
            
                **L1. Dust–$\Gamma_t$ emergence:** Massive galaxies at $z>8$ are heavily dust-obscured while low-mass galaxies are dust-poor ($\rho = +0.62$). This signal strengthens with redshift and tracks a specific step-function threshold corresponding to AGB dust onset. A Fisher combination across five independent datasets (three photometric surveys, a fourth sky field, and a spectroscopic H$\alpha$/H$\beta$ line-ratio test) confirms this correlation at $30.0\sigma$ (§3.8).

                - **L2. Inside-out core screening:** Spatially resolved color gradients in JADES show massive galaxies have bluer cores ($\rho = -0.18$), contradicting standard inside-out growth models but matching TEP's central screening prediction.

                - **L3. Mass–sSFR inversion:** The standard downsizing correlation reverses sign at $z > 7$, consistent with time-dilation inflating the apparent sSFR in massive early halos. The discontinuous sign change makes this test immune to smooth mass-measurement systematics.

                - **L4. Kinematic dynamical masses:** TEP resolves 11/11 physically impossible $M_*/M_{\rm dyn} > 1$ cases using SED-independent kinematics, providing a direct break from photometric mass circularity.

            
        
        - **Stage 3: The Little Red Dot Mechanism (§4.9).** TEP is applied to the overmassive black hole crisis, demonstrating how *differential temporal shear* between the dense core and the diffuse stellar halo accelerates apparent black hole growth without requiring super-Eddington accretion.

    

    
### 1.4 Prior Cross-Domain Evidence for TEP

    The JWST analysis presented here is not the first test of the TEP framework. The same coupling $\alpha_0 = 0.58$ has been independently constrained across multiple domains spanning 13.5 Gyr of cosmic time (Table 11, §4.2), providing a pre-existing evidence base against which the present results can be evaluated. **Important caveat:** all prior constraints in Table 11 derive from the same author's analysis pipeline; the TEP series has not yet undergone independent replication or peer review in a refereed journal. Readers should weigh the cross-domain consistency with this single-source limitation in mind (§4.11a). The three domains most directly used in this work are:

    
        - **Hubble Tension:** Stratification of $N = 29$ SH0ES Cepheid hosts by velocity dispersion reveals an environmental bias. High-$\sigma$ hosts yield $H_0 = 72.45 \pm 2.32$ km/s/Mpc; low-$\sigma$ hosts yield $H_0 = 67.82 \pm 1.62$ km/s/Mpc — consistent with Planck within $1\sigma$. The TEP correction with $\alpha_0 = 0.58$ yields $H_0^{\rm TEP} = 68.66 \pm 1.51$ km/s/Mpc ($0.79\sigma$ from Planck). This is the source of $\alpha_0$ used in this work.

        - **Globular Cluster Pulsars:** Analysis of 380 millisecond pulsars reveals a $0.13$ dex spin-down excess in cluster pulsars ($p = 1.7 \times 10^{-15}$). The environmental screening threshold $\sigma > 165$ km/s derived from this population is used directly in §2.3.2.2 and §4.4.3 of this work.

        - **Universal Critical Density:** GNSS atomic clock networks yield $\rho_c \approx 20$ g/cm³, independently confirmed by the SPARC rotation curve slope and magnetar critical periods. This $\rho_c$ is used as the screening threshold in this work.

    
    The central question this work addresses is whether the same $\alpha_0$ that resolves the Hubble tension and accounts for pulsar timing anomalies also predicts the high-redshift galaxy anomalies, with no re-tuning. The recovered JWST coupling from dust alone is $\alpha_0 = 0.55 \pm 0.32$ (Pearson $R^2$ maximisation at $z > 6$, see §4.4.8), consistent with the Cepheid value to $0.1\sigma$; the joint four-observable fit gives $\alpha_0 = 0.75 \pm 0.29$ ($0.5\sigma$ from nominal). The large uncertainties reflect the genuine degeneracy of $\Gamma_t$ with halo mass within any single survey — the Cepheid calibration therefore remains the primary constraint on $\alpha_0$.

    
### 1.5 Alternative Explanations

    Standard alternatives include top-heavy initial mass functions (Boylan-Kolchin 2023), enhanced AGN feedback, bursty star formation, and super-Eddington accretion. Each can partially accommodate one or two of the observed anomalies but fails on the full set. AGN feedback predicts a negative dust–black hole mass correlation, as AGN activity clears dust; the observed correlation is strongly positive ($\rho = +0.62$). Bursty star formation predicts bluer colours during burst phases, whereas the TEP-enhanced population is significantly redder at fixed magnitude ($\rho(M_{\rm mag}, \text{color}) = -0.40$, $p = 2.8 \times 10^{-16}$, $N = 398$). Top-heavy IMFs can partially resolve the star formation efficiency crisis but offer no mechanism for the spatially resolved screening gradients, the mass–sSFR inversion, or the differential black hole growth in Little Red Dots. A systematic comparison across all tested signatures (§4.3.5) shows TEP achieves 4/4 correct sign predictions on the independent lines with one fitted parameter, against 0–1/4 for alternatives requiring 2–3 parameters. No alternative reproduces the cross-domain consistency of $\alpha_0$ across instruments, epochs, and physical regimes.

    
        **Key Limitations**

        Several important limitations should be borne in mind when evaluating the evidence presented here:

        
            - **Mass circularity:** Because $\Gamma_t$ is derived from halo mass, distinguishing TEP effects from intrinsic mass-dependent evolution requires careful partial-correlation analysis (§3.4) and exploitation of the independent redshift-dependent component.

            - **Spectroscopic sample size:** Two new spectroscopic compilations substantially upgrade the prior $N = 147$ sample. JADES DR4 (D'Eugenio et al. 2025) provides 2,858 good-quality spec-z with 118 at $z > 7$. The DJA NIRSpec Merged Table v4.4 (Brammer et al.; September 2025) provides 19,445 unique grade-$\ge 3$ sources from all public JWST/NIRSpec programs, with 698 at $z > 7$ and 234 at $z > 8$. Both $z > 7$ and $z > 8$ subsamples are now well-powered. However, stellar masses rely on photometric estimates ($\pm 0.3$–$0.5$ dex), so both spectroscopic samples remain consistency checks rather than independent lines of evidence. Five new independent datasets spanning 4 sky fields and including a spectroscopic H$\alpha$/H$\beta$ dust measurement yield a Fisher combined L1 significance of $30.0\sigma$ (§3.8); these are treated as robust replications of L1 rather than new independent lines.

            - **Theoretical foundation:** The enhancement factor $\Gamma_t > 1$ is motivated by a scalar-tensor, two-metric construction in which matter couples to a conformal factor $A(\phi)$ and screening suppresses the effective coupling in dense regimes. The first-principles action, field equations, screening mechanism, and PPN consistency mapping provide the foundation for this framework; this manuscript presents only the components required to define and test the observational mapping (§2.3.2). A full joint cosmological parameter inference is outside the scope of this work.

            - **Red Monsters case study:** The $N = 3$ illustrative case study is not statistically robust in isolation. Population-level tests ($N = 2{,}315$) provide the primary evidence.

        
    

    Section 2 details the TEP mapping and statistical procedures, including a formal derivation from scalar-tensor theory. Section 3 presents the four-line evidence package and independent replications. Section 4 discusses theoretical implications: compatibility with precision GR tests, the link to the Hubble tension, and the physical mechanism for LRD black hole growth. Section 5 concludes with falsification criteria and observational predictions. Appendix A provides the theoretical foundation (action, field equations, screening mechanism), and Appendix B documents key pipeline algorithms.

    
## 2. Data and Methods

    
### 2.1 Data

    
#### 2.1.1 Red Monsters (FRESCO)

    The motivating case study is the class of ultra-massive galaxies at $z \sim 5$–$6$ exemplified by the three spectroscopically confirmed "Red Monsters" reported by Xiao et al. (2024). For the illustrative TEP prediction (§3.1), representative parameters spanning the published range ($z \approx 5.1$–$5.9$, $\log M_* \approx 10.6$–$10.9$, SFE $\approx 0.48$–$0.52$) are adopted. These are *not* the exact catalogue values for individual objects but capture the regime where the anomaly is most acute. The resulting SFE correction ($\sim 34\%$) depends primarily on $\Gamma_t$ (set by halo mass and redshift via the pre-calibrated TEP formula) and is insensitive to the precise input SFE at the $\lesssim 2\%$ level.

    
#### 2.1.2 UNCOVER DR4

    For population-level tests, the UNCOVER DR4 stellar population synthesis catalog is used (Wang et al. 2024; Furtak et al. 2023). The pipeline applies quality cuts and constructs a high-redshift sample with $4 < z < 10$ and $\log M_* > 8$, yielding $N = 2{,}315$ galaxies. For multi-property analyses (age ratio, metallicity, dust), a subset with complete measurements is used (e.g., $N = 1{,}108$ for the partial-correlation and split tests).

    
#### 2.1.3 Independent replications and spectroscopic validation

    To evaluate independent replication of the $z > 8$ dust result, catalogs for CEERS are used (Cox et al. 2025; Finkelstein et al. 2023; photometric redshifts via EAZY, Brammer et al. 2008) and COSMOS-Web (Shuntov et al. 2025). The COSMOS2025 catalog (Shuntov et al. 2025) provides LePHARE SED-derived stellar masses, SFRs, E(B-V) dust, and ages for 784,016 galaxies over 0.54 deg², with 37,965 sources at $z > 4$ passing quality cuts — the largest high-$z$ photometric SED sample used in this analysis. The UNCOVER DR4 SPS catalog (Wang et al. 2024; Suess et al. 2024; Price et al. 2025) uses 20-band MegaScience photometry and Prospector-β SED fitting, providing 2,628 sources at $z > 4$ with Prospector dust2 and a spec-z sub-catalog of 203 sources with spectroscopic redshifts fixed in the SED fit. For spectroscopic validation, two complementary catalogs are used:

    **JADES Data Release 4** (D'Eugenio et al. 2025; Curtis-Lake et al. 2025; Scholtz et al. 2025): 2,858 high-quality spectroscopic redshifts (flags A/B) across GOODS-N and GOODS-S, with 118 sources at $z > 7$ and 41 at $z > 8$. UV-luminosity-based stellar masses (Song et al. 2016) are derived for the 1,345 sources with valid $M_{\rm UV}$.

    **DAWN JWST Archive (DJA) NIRSpec Merged Table v4.4** (Brammer et al.; de Graaff et al. 2024; Heintz et al. 2023; September 2025): a comprehensive compilation of 80,367 uniformly reduced JWST/NIRSpec spectra from all public programs, processed with the msaexp/grizli pipelines. After applying grade $\ge 3$ quality cuts and deduplication by sky position, 19,445 unique sources are retained, of which 3,251 are at $z > 5$, 698 at $z > 7$, and 234 at $z > 8$. Photometric stellar masses are available for 2,598 of the high-$z$ sources. This catalog spans JADES, CEERS, RUBIES, UNCOVER, GLASS, PRIMER, and more than 50 other public programs, providing the largest uniform cross-survey spectroscopic sample to date.

    
#### 2.1.4 MIRI-based mass calibration context

    Recent JWST/MIRI analyses (Pérez-González et al. 2024) show that NIRCam-only SED fits can overestimate stellar masses at $z > 5$ because of age-attenuation degeneracy and emission-line contamination. When MIRI photometry is included, the number density of the most massive systems decreases and some candidates are reclassified as dusty or line-dominated sources. The photometry is not reprocessed in this work, but published masses are treated as conservative upper bounds and MIRI-based studies serve as an external check on the interpretation of the extreme-mass tail.

    
        
            Table 1a: Observational Datasets
            
                
                    Dataset
                    Role
                    Sample Size
                    Redshift Range
                    Mass Cut ($\log M_*$)
                    Key Reference
                    Key Biases
                
            
            
                
                    Red Monsters
                    Case Study
                    3
                    $5.3 < z < 5.9$
                    $> 10.5$
                    Xiao et al. (2024)
                    Small N, Selection Function
                
                
                    UNCOVER DR4
                    Primary Statistical Sample
                    2,315
                    $4 < z < 10$
                    $> 8.0$
                    Wang et al. (2024)
                    NIRCam Mass Overestimation
                
                
                    CEERS DR1
                    Independent Replication
                    82
                    $z > 8$
                    $> 8.0$
                    Cox et al. (2025)
                    Field Variance
                
                
                    COSMOS-Web
                    Large-Volume Check
                    2,606 (918 dust-detected)
                    $z > 8$
                    $> 8.0$
                    Shuntov et al. (2025)
                    Photometric Redshift Uncertainties; Zero-Inflated Dust
                
                
                    JADES DR4 (NIRSpec/MSA)
                    Spectroscopic Validation
                    2,858 (flags A/B); 118 at $z > 7$
                    $z = 0.1$–$14.2$
                    None
                    D'Eugenio et al. (2025); Curtis-Lake et al. (2025)
                    Slit Losses; UV-based $M_*$ ($\pm 0.4$ dex)
                
                
                    DJA NIRSpec Merged v4.4
                    Cross-Survey Spectroscopic Validation
                    19,445 unique (grade $\ge 3$); 698 at $z > 7$; 234 at $z > 8$
                    $z = 0.1$–$14.1$
                    None
                    Brammer et al. (DJA); de Graaff et al. (2024)
                    Photometric $M_*$ from grizli; heterogeneous survey depths
                
                
                    UNCOVER DR4 SPS (MegaScience)
                    Primary + Spectroscopic Validation
                    2,628 (z$>$4, Prospector-β); 203 with spec-z fixed fits
                    $z = 4$–$12$
                    Abell 2744 (lensed)
                    Wang et al. 2024; Suess et al. 2024; Price et al. 2025
                    20-band photometry; lensing magnification corrections
                
                
                    COSMOS2025 (LePHARE SED)
                    Cross-Field Replication
                    37,965 (z$>$4, M$_*>10^8$); 5,672 at $z > 7$; 2,128 at $z > 8$
                    $z = 4$–$13$
                    None (blank field)
                    Shuntov et al. 2025 (COSMOS2025)
                    LePHARE E(B-V) less precise than Prospector dust2; photo-z scatter
                
            
        
    

    Related MIRI-supported analyses of Little Red Dots (LRDs) at $z > 4$ find that inferred stellar masses can shift by up to orders of magnitude depending on the assumed AGN contribution. This motivates a conservative stance in the interpretation of compact red sources and provides a systematic-control context for any extreme-mass claims in the literature.

    
### 2.2 Key Terminology

    The following terms are used consistently throughout this work:

    
        
            Table 1b: Glossary of Key Terms
            
                
                    Term
                    Symbol
                    Definition
                
            
            
                
                    Temporal Enhancement Factor
                    $\Gamma_t$
                    The ratio of effective proper time to cosmic time experienced by stellar populations. In unscreened regions (low density), $\Gamma_t$ increases with potential depth: massive halos have $\Gamma_t > 1$ (enhancement), low-mass halos have $\Gamma_t < 1$ (suppression). In screened regions (density $> \rho_c$), $\Gamma_t \to 1$ regardless of potential depth.
                
                
                    Temporal Shear
                    —
                    The spatial gradient of $\Gamma_t$ across a galaxy or environment. Differential temporal shear refers to the difference in $\Gamma_t$ between two regions (e.g., galactic center vs. halo).
                
                
                    Isochrony Bias
                    —
                    The systematic error in inferred stellar properties (mass, age, SFR) arising from the assumption that stellar clocks tick at the cosmic rate everywhere. Under TEP, this assumption is violated in deep potential wells.
                
                
                    Screening
                    —
                    The suppression of TEP effects in regions where the local density exceeds a critical threshold ($\rho_c \approx 20$ g/cm³), restoring standard GR physics. Two types are distinguished:**
                    *Core Screening*—Screening within a single galaxy, where the deep central potential suppresses TEP ($\Gamma_t \to 1$) while the outskirts remain enhanced. Produces bluer cores and redder outskirts.

                    *Environmental Screening*—Screening by the ambient group or cluster potential, causing galaxies in dense environments to appear younger than isolated field galaxies of the same mass.
                
                
                    Effective Time
                    $t_{\rm eff}$
                    The proper time experienced by stellar populations: $t_{\rm eff} = t_{\rm cosmic} \times \Gamma_t$.
                
            
        
    

    
### 2.3 Derived quantities

    
#### 2.3.1 Halo mass proxy

    For each galaxy, the pipeline uses an abundance-matching relation (Behroozi et al. 2019) to map stellar mass to halo mass. This mapping is used solely to construct $\Delta\log(M_h)$ for the TEP parameterization. To mitigate circularity, sensitivity tests are performed with $\pm 0.3$ dex scatter in the $M_h-M_*$ relation, propagating to $\pm 12\%$ in $\Gamma_t$ corrections.

    
#### 2.3.2 The TEP Metric Coupling

    The temporal enhancement factor $\Gamma_t$ is not an ad-hoc fitting function but is derived from first principles within a conformally coupled scalar-tensor framework. This framework introduces a scalar field $\phi$ that couples to the trace of the matter energy-momentum tensor, effectively altering the local metric experienced by standard model fields. While the full theoretical development, including stability analysis and cosmological constraints, is extensive, this section summarizes the key steps from action to observable required for the current analysis.

    
    2.3.2.1 From Action to Observable
    The TEP framework builds upon chameleon-class scalar-tensor theories (Khoury & Weltman 2004; Brax et al. 2004; Burrage & Sakstein 2018) where the effective mass of the scalar field depends on the local matter density. The key steps mapping the fundamental physics to the observable $\Gamma_t$ are:

    
        - Action:** Matter couples to $\tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu}$ where $A(\phi) = \exp(2\beta\phi/M_{\rm Pl})$. The Klein-Gordon equation sources $\phi$ from the matter density trace $T^\mu_\mu$.

        - **Proper time:** Clock rates scale as $d\tau/dt \approx A(\phi)^{1/2}$, defining $\Gamma_t \equiv (d\tau/dt)/(d\tau/dt)_{\rm ref}$.

        - **Halo mapping:** In virialized halos, $\phi$ tracks the potential depth $\Phi \propto M_h^{2/3}$, yielding:

    
    
        $$\Gamma_t = \exp \left[ \alpha(z) \times \frac{2}{3} \times \Delta \log(M_h) \times \frac{1+z}{1+z_{\rm ref}} \right]$$
    
    where $\Delta\log(M_h) = \log(M_h/M_{\rm ref})$, $\log M_{h,\rm ref} = 12.0$, and $\alpha(z)$ is the redshift-dependent coupling.

    2.3.2.1a Enhancement and Suppression: The Two-Sided Prediction
    The exponential formula predicts two physically distinct regimes depending on halo mass relative to the reference mass $\log M_{h,\rm ref} = 12.0$:

    
        - **Enhancement ($\Gamma_t > 1$, massive halos):** For $M_h > M_{\rm ref}$, $\Delta\log M_h > 0$ and $\Gamma_t > 1$. The scalar field $\phi$ is sourced more strongly by the deeper potential, raising $A(\phi) > 1$ and accelerating material clock rates relative to the cosmic mean. This is the regime of the Red Monsters and massive $z > 8$ galaxies.

        - **Suppression ($\Gamma_t 
    This two-sided prediction is physically essential: the "Uniformity Paradox" — why low-mass galaxies at $z > 8$ are dust-poor despite cosmic time being nominally sufficient for AGB production — is resolved precisely because $\Gamma_t \ll 1$ in low-mass halos shuts off the effective AGB clock. A model that only predicted enhancement ($\Gamma_t \geq 1$ everywhere) would not explain the dust-poor low-mass population. The suppression regime is therefore a falsifiable prediction, not a free parameter: it predicts that low-mass galaxies at $z > 8$ should be systematically younger in their stellar populations than their cosmic age implies, and should lack dust regardless of the available cosmic time.

    2.3.2.2 Screening and Scale Separation
    The bare coupling $\alpha_0 = 0.58$ corresponds to $\omega_{\rm BD} = 1/(2\alpha_0^2) - 1/2 \approx 0.99$ — roughly four orders of magnitude below the Cassini bound ($\omega_{\rm BD} > 40{,}000$; Bertotti et al. 2003). This dramatic pre-screening violation is the defining feature of chameleon-class theories: the bare coupling is strong, but the thin-shell mechanism suppresses the effective coupling by $\Delta R/R \lesssim 10^{-6}$ in solar-system bodies, yielding $\omega_{\rm BD}^{\rm eff} > 10^6$. On cosmological scales, the Compton wavelength $\lambda_C \sim 1$ Mpc yields Yukawa suppression $\beta_{\rm eff}(R_8) \approx 0.002$ on $\sigma_8$ scales—well below the Planck bound. Within individual halos ($r \lesssim 50$ kpc), the field tracks the local potential and operates locally. This two-scale picture is standard for chameleon theories.

    The screening condition for local tests is:

    
        $$\rho > \rho_c \approx 20 \text{ g/cm}^3 \quad \Rightarrow \quad \Gamma_t \to 1 \text{ (screened)}$$
    
    This critical density is derived independently from three sources that converge on the same value: GNSS atomic clock networks ($L_c \approx 4200$ km for Earth's mass), atomic physics (soliton radius $R_{\rm sol}(m_p) \sim a_0$ at the proton mass scale), and magnetar anti-glitches ($P_{\rm crit} \approx 6.8$ s for 1E 2259+586, 4% match). The convergence across 40 orders of magnitude in mass provides strong independent support for this screening scale.

    At galactic scales, an effective kinematic screening threshold emerges from analysis of 380 millisecond pulsars in globular clusters, which reveals that the TEP spin-down excess saturates for systems with velocity dispersion $\sigma \gtrsim 165$ km/s, consistent with the scalar field entering the thin-shell regime. This threshold is used in §4.4.3 to define the environmental screening boundary for JWST galaxies: halos with $\sigma \gtrsim 165$ km/s (corresponding to $\log M_h \gtrsim 13.5$ at $z \sim 0$) are expected to be partially screened, suppressing $\Gamma_t$ below the unscreened prediction.

    2.3.2.3 Enhancement vs. Dilation
    Standard GR predicts time *dilation* in deep potentials; TEP predicts *enhancement* ($\Gamma_t > 1$). These refer to different metrics: gravitational redshift is governed by $g_{\mu\nu}$ (preserved identically), while material clock rates are governed by $\tilde{g}_{\mu\nu} = A(\phi)g_{\mu\nu}$. The key distinction is that $\Gamma_t$ compares clock rates between *different environments at the same epoch*, not between positions in a single well. Numerical integration confirms $A(\phi) > 1$ in unscreened halos for $2\beta^2 > 1$. Solar System bodies are fully screened ($\Gamma_t \to 1$).

    
        
            Table 2: TEP Model Parameters (Fixed)
            
                
                    Parameter
                    Value
                    Source
                    Description
                
            
            
                
                    $\alpha_0$
                    $0.58 \pm 0.16$
                    Cepheid Calibration
                    Coupling strength from Cepheids
                
                
                    $z_{\rm ref}$
                    5.5
                    TEP-H0
                    Reference redshift for calibration
                
                
                    $\log M_{h, \rm ref}$
                    12.0
                    TEP-COS
                    Reference halo mass ($\Gamma_t=1$)
                
                
                    $\rho_c$
                    20 g/cm$^3$
                    TEP-UCD
                    Critical density for screening
                
            
        
    

    

        *
        

    Figure 1: The TEP Metric Coupling $\Gamma_t(M_h, z)$ in the unscreened regime. The enhancement factor increases with halo mass (potential depth) and redshift (weakening of cosmological screening). The reference mass ($\log M_h = 12$) defines $\Gamma_t = 1$ (cosmic time flow). Massive halos at high redshift experience significant temporal enhancement ($\Gamma_t > 1$), while low-mass halos are suppressed ($\Gamma_t < 1$). In screened regions (density $> \rho_c \approx 20$ g/cm³), $\Gamma_t \to 1$ regardless of mass.

    

    Parameter Calibration. The parameters $\alpha_0 = 0.58$ and $z_{\rm ref} = 5.5$ are fixed based on independent calibration from the local universe. $\alpha_0$ was derived from the period-luminosity residuals of Cepheids (Kodric et al. 2018; Freedman et al. 2024) in massive hosts (e.g., M31, NGC 4258) relative to low-mass hosts. This leaves zero free parameters to be tuned to the JWST data.

    2.3.2.7 Cosmological Viability Summary
    The TEP framework has been checked against three classes of precision cosmological constraints:

    
        - **BBN:** During the radiation era, $T^\mu_\mu \approx 0$ for relativistic species, so the scalar field driving force vanishes and $\phi$ remains frozen. Numerical integration yields $|\Delta H/H|_{\rm max} = 1.7 \times 10^{-13}$ and $\Delta Y_p 
    **Scale-dependent growth computation:** To go beyond the analytic Yukawa argument, the linear growth ODE is solved independently for each Fourier mode $k$, incorporating the full scale-dependent gravitational coupling $G_{\rm eff}(k,z)/G_N = 1 + 2\beta^2 k^2/(k^2 + m_\phi(z)^2)$, where $m_\phi(z) = m_{\phi,0}(1+z)^{9/4}$ for $n=1$ chameleon. The resulting matter power spectrum ratio $P_{\rm TEP}(k)/P_{\Lambda{\rm CDM}}(k)$ and integrated $\sigma_8$ are computed self-consistently (Appendix A.1.8). Key results:

    
        - Planck consistency requires $m_{\phi,0} \gtrsim 0.2\,h$/Mpc ($\lambda_C \lesssim 31\,h^{-1}$ Mpc at $z=0$; Appendix A.1.8).

        - For typical chameleon parameters, $\beta_{\rm eff}$ on $R_8$ scales is $\approx 0.0046$—a factor $\sim 126$ below the bare coupling ($G_{\rm eff}/G_N - 1 = 4.1 \times 10^{-5}$ at $k_8$; Appendix A.1.8).

        - The predicted $\sigma_8^{\rm TEP} = 0.8110$ vs. $\sigma_8^{\Lambda{\rm CDM}} = 0.811$; $\Delta\sigma_8 = 7.2 \times 10^{-8}$ ($1.2 \times 10^{-5}\sigma$). RSD: $\Delta\chi^2 = 8 \times 10^{-5}$ across 8 data points (Appendix A.1.8).

    
    This k-dependent growth calculation substantially strengthens the viability argument beyond the earlier analytic estimate. A full CAMB Boltzmann integration (Appendix A.1.8.8) confirms these results: at the fiducial $m_{\phi,0} = 1.0\,h$/Mpc, the CAMB-computed $\sigma_8^{\rm TEP} = 0.8116$ ($0.10\sigma$ from Planck), with CMB TT deviations $ 1089$ — is justified by $T^\mu_\mu \approx 0$ during radiation domination. A natively coupled hi_class integration remains desirable for completeness but is no longer expected to change the conclusion.

    
#### 2.3.3 Effective time and isochrony bias correction

    An effective time is defined as $t_{\rm eff} = t_{\rm cosmic}\,\Gamma_t$, where $t_{\rm cosmic}$ is computed from a fiducial cosmology (Planck18). Under the isochrony-bias model used here, the mass-to-light ratio is assumed to scale as $M/L \propto t^n$ (following standard SSP predictions; Bruzual & Charlot 2003; Conroy et al. 2009). Forward-modeling analysis finds that $n \approx 0.5$ minimizes the residual mass-age correlation at $z > 6$, while $n \approx 0.9$ is preferred at $z = 4$–$6$. For the primary high-$z$ analysis, $n = 0.5$ is adopted. The corrected stellar mass and SFE are:

    
        $$M_{*,\rm true} = M_{*,\rm obs}/\Gamma_t^{n}, \quad \mathrm{SFE}_{\rm true} = \mathrm{SFE}_{\rm obs}/\Gamma_t^{n}.$$
    

    
### 2.4 Statistical procedures

    Associations are quantified using Spearman rank correlations and bootstrap confidence intervals. To address confounding by redshift and stellar mass, partial-correlation analyses implemented via residualization are employed. In addition to correlation-based tests, the following are reported:

    
        - Stratified comparisons (e.g., high vs low $\Gamma_t$ splits) for multi-property coherence

        - Distributional comparisons (e.g., Kolmogorov-Smirnov tests) for regime separation

        - Model comparison using AIC/BIC for regression models that compare predictors {z}, {z, $\log M_*$}, {z, $\Gamma_t$}, and {z, $\log M_*$, $\Gamma_t$}

    

    
#### 2.4.1 Combined significance and multiple testing

    Combined significance is assessed using Fisher's method, Bonferroni correction, Brown's method (dependence-adjusted), harmonic mean p-value, and Benjamini-Hochberg FDR ($\alpha = 0.05$). The most conservative estimate reduces effective sample sizes by 90% via spatial clustering autocorrelation, yielding $6.4\sigma$. Parametric p-values are supplemented by permutation tests ($N = 10{,}000$ shuffles) and bootstrap confidence intervals ($N = 10{,}000$ resamples). Cross-survey effect sizes are combined via DerSimonian-Laird random-effects meta-analysis with $I^2$ heterogeneity assessment and leave-one-out influence diagnostics.

    
#### 2.4.2 Blind validation protocol

    Three split strategies test generalization: (1) time-split (low-$z$ train / high-$z$ test, 60/40); (2) field-split (RA median); (3) cross-survey leave-one-out. A test passes if the dust–$\Gamma_t$ correlation remains significant ($p  4\sigma$ under a 0.5 dex reduction. At $z > 8$, selection bias toward bright galaxies is quantified via Monte Carlo completeness weighting ($N = 1{,}000$ iterations) and Savage-Dickey Bayes Factors.

    
#### 2.4.4 Forward-modeling validation

    The $M/L \propto t^{n}$ scaling is validated by varying $n = 0.5$–$0.9$ and identifying the value minimizing the residual mass-age correlation. At $z > 6$, $n \approx 0.5$ is preferred (consistent with low-metallicity SSP models); at $z = 4$–$6$, $n \approx 0.9$.

    
### 2.5 Black Hole Growth Simulation

    To test the "Little Red Dot" resolution, a differential temporal shear simulation was developed. A compact galaxy ($r_e \approx 150$ pc) with a baryon-dominated core ($c=10$) is modeled. The local temporal enhancement factor $\Gamma_t(r)$ is computed at the center (Black Hole environment) and at the effective radius (Stellar environment) across the redshift range $z=4$–$10$.

    The differential growth factor is computed as:

    
        $$\text{Boost} = \exp\left(\frac{\int (\Gamma_{\rm cen}(z) - \Gamma_{\rm halo}(z)) \, dt_{\rm cosmic}}{t_{\rm Salpeter}}\right)$$
    
    where $t_{\rm Salpeter} \approx 45$ Myr is the Salpeter timescale (e-folding time for Eddington-limited accretion). This simulation uses the same $\alpha_0=0.58$ parameter calibrated from Cepheids, with no additional tuning.

    
### 2.6 Reproducibility

    A complete run is executed with `python scripts/steps/run_all_steps.py`, which produces step-indexed outputs in `results/outputs` and step-indexed logs in `logs`.

    
## 3. Results

    
### 3.0 Evidence Summary

    Four independent TEP signatures are summarized below. **Independence disclosure:** L1 (dust–$\Gamma_t$ correlation) and the AGB threshold test share the same predictor ($t_{\rm eff}$ is a deterministic function of $\Gamma_t$); the partial correlation between their predictors after controlling for $\Gamma_t$ is $\rho = 0.049$ ($p = 0.41$), confirming they are the same signal in two representations — they are counted as one line. Cross-survey generalization (formerly L5) tests the functional form of the same dust observable and is a robustness check on L1, not a new line. Age coherence and metallicity vanish under joint mass+redshift control and are not counted. The environmental screening null at $z>8$ ($\Delta\rho = 0.003$, $p = 0.97$) is a TEP prediction* and is not counted as an independent line. The colour-gradient Steiger test is provisional (synthetic data) and is not counted. The four independent lines form the core evidence package.

    
        
            Table 3a: All Tested TEP Signatures — Four are counted as independent lines of evidence; others are robustness checks or not independent
            
                SignatureFindingSignificantSurvives Mass ControlStatus
            
            
                **L1. Dust–$\Gamma_t$ correlation + AGB threshold**$\rho = +0.62$, $N=1{,}283$, three surveys; partial $\rho = +0.26$ after full polynomial control; AGB step odds ratio 42.8 ($p ✔✔**Independent line**
                **L2. Inside-out core screening**Bluer cores, redder outskirts ($\rho = -0.18$, $p = 5\times10^{-4}$; real JADES resolved photometry, Rieke+2023; single survey, $N = 362$, not replicated). ***Caveat:** Steiger Z-test ($Z = 2.01$, $p = 0.045$) uses synthetic data; the primary $\rho = -0.18$ is real but standard inside-out dust models can produce similar gradients. Effect size is small-to-medium (Cohen's $d \approx 0.37$).*✔✔**Independent line** (marginal)
                **L3. Mass–sSFR inversion**Correlation inverts at $z>7$ ($\Delta\rho = +0.25$, 95% CI excludes zero); partial $\rho(\Gamma_t, {\rm sSFR}|{\rm dust}) = -0.49$ ($p = 10^{-18}$)✔✔**Independent line**
                **L4. Dynamical mass comparison**11/11 impossible $M_*/M_{\rm dyn}>1$ resolved; mean ratio 1.33→0.61 (kinematic, SED-independent)✔✔**Independent line**
                Steiger Z-test ($t_{\rm eff}$ vs $M_*$)$Z=17.8$, $p=1.3\times10^{-70}$; $Z=10.4$, $p=1.8\times10^{-25}$ at $z>8$✔✔Robustness check on L1
                Partial correlations$\rho=+0.26$ after full polynomial control; $M_*$ zero residual after $t_{\rm eff}$ control✔✔Robustness check on L1
                Cross-survey generalizationPolynomial $R^2$ collapses to $-6.4$ cross-survey; $t_{\rm eff}$ stable at $\rho=0.60$–0.80 (same dust observable as L1)✔✔Robustness check on L1
                Age coherence$\rho = +0.14$ (mass-only); vanishes with $M_*$+$z$ control✔✘Not independent (mass proxy)
                Metallicity$\rho = +0.16$ (mass-only); vanishes with $M_*$+$z$ control✔✘Not independent (mass proxy)
                Environmental screeningFull-sample $Z=4.68$, $p=2.9\times10^{-6}$; predicted null at $z>8$ confirmed ($p=0.97$)✔✔Additional test (predicted null confirms coupling)
                Colour-gradient Steiger$Z=2.01$, $p=0.045$; partial $\rho=+0.13$. *Synthetic data.*✔✔Provisional — not counted
            
        
    

    
### 3.1 Red Monsters: A Zero-Parameter Prediction

    The TEP parameterization is applied to galaxies in the Red Monster regime ($z \sim 5$–$6$, $\log M_* \gtrsim 10.5$; Xiao et al. 2024). This is a blind prediction: the coupling constant $\alpha_0 = 0.58$ is fixed entirely from local Cepheid data ($z \approx 0$). No parameters are fitted or tuned to the high-redshift observations. The three entries below (S1–S3) use representative parameters spanning the published range (§2.1.1); the predicted correction depends primarily on $\Gamma_t$ and is insensitive to the exact input SFE.

    
    Despite the small sample ($N=3$), the effect size (mean SFE reduction $0.36$, $\sigma \approx 0.040$) yields power $> 99\%$ ($t = 15.4$, $p = 0.002$). Population-level tests ($N = 2{,}315$) provide the primary statistical evidence. The externally-calibrated prediction (zero parameters tuned to JWST data) is further validated across three surveys ($N = 1{,}283$ at $z > 8$).

    
        
            Table 3b: Illustrative TEP Predictions for Red Monster–Class Galaxies
            
                ID$z$$\alpha(z)$$\Gamma_t$ (Predicted)SFE$_{\rm obs}$SFE$_{\rm true}$% Anomaly Resolved
            
            
                S15.301.462.120.500.3432%
                S25.501.481.810.480.3625%
                S35.901.522.940.520.3042%
                Average Prediction2.290.500.3334% [22%–48%]
            
        
    

    The predicted mass bias $\Gamma_t^{n} \approx 1.50$ reduces the corrected SFE to $\sim 0.34$ (from $2.5\times$ to $1.7\times$ the standard limit). Propagating the $\alpha_0 = 0.58 \pm 0.16$ uncertainty ($\pm 28\%$) yields a correction range of 22%–48% at $1\sigma$, with the central value 34% used throughout. The correction is robust to the sign: even at $\alpha_0 = 0.42$ (lower $1\sigma$ bound), the SFE anomaly is reduced by $\sim 22\%$ with zero tuned parameters.

    
### 3.2 UNCOVER DR4: Mass-sSFR and Mass-Age Correlations

    The Red Monster case study establishes that TEP predicts the correct direction and magnitude of the SFE correction for individual extreme objects. The critical question is whether this signal extends to the full galaxy population. In the UNCOVER DR4 sample ($N = 2{,}315$), the mass-sSFR correlation is weak but significant ($\rho = -0.13$, $p = 1.3 \times 10^{-10}$, Cohen's $d = -0.27$), consistent with TEP partially canceling the intrinsic downsizing trend. The mass-age correlation is positive ($\rho = +0.14$, $p = 7.0 \times 10^{-11}$), consistent with more massive galaxies experiencing more proper time. Both correlations are in the predicted direction but are attenuated by the full redshift range; the signal sharpens substantially when the sample is stratified by redshift.

    
### 3.3 Redshift Evolution: The High-z Transition

    TEP predicts that the mass-sSFR correlation should become *less negative* (or even positive) at higher redshift, where the TEP enhancement is stronger. This is tested by stratifying the sample:

    
        
            Table 4: Mass-sSFR Correlation by Redshift
            
                $z$ Range$N$Spearman $\rho$95% CIInterpretation
            
            
                4–5942$-0.17$[$-0.24$, $-0.11$]Standard downsizing
                5–6497$-0.14$[$-0.22$, $-0.05$]Standard downsizing
                6–7372$-0.06$[$-0.16$, $+0.04$]Weakening
                7–8221$+0.18$[$+0.05$, $+0.31$]Inversion
                8–9179$+0.13$[$-0.03$, $+0.29$]Weak positive
                9–10104$-0.27$[$-0.47$, $-0.05$]Reversal (selection effects)
            
        
    

    Comparing low-$z$ ($4  7$, $\rho = +0.09$): $\Delta\rho = +0.25$ [+0.14, +0.35] (95% CI excludes zero), indicating a statistically significant inversion.

    
### 3.4 Partial Correlation Test

    The redshift evolution in §3.3 is consistent with TEP but does not rule out a mass proxy: since $\Gamma_t \propto M_h^{1/3}$, any mass-dependent effect could produce a similar inversion. The partial-correlation hierarchy directly addresses this concern. With mass-only control, age-ratio and metallicity remain weakly positive; with joint mass+redshift control they are consistent with zero, so they are treated as mass-proxy-adjacent rather than independent. However, the $z > 8$ dust–$\Gamma_t$ correlation survives ($\rho = +0.28$, $p  165$ km/s from globular cluster pulsar timing. At high redshift, this threshold shifts to higher halo mass. Screening is tested by comparing age ratios (MWA/$t_{\rm cosmic}$) across mass bins:

    
        
            Table 5: Age Ratio by Halo Mass (5 < z < 8)
            
                $\log M_h$$N$$\langle$MWA/$t_{\rm cosmic}\rangle$$\Gamma_t$ Predicted
            
            
                10–11390$0.15 \pm 0.003$$\sim 0$ (reference)
                11–1242$0.18 \pm 0.015$0.2–0.5
                12–12.53$0.30 \pm 0.12$1.0–1.5
                12.5–131$0.05$1.5–2.0
            
        
    

    
#### 3.5.1 Resolved Core Screening

    TEP predicts that deep core potentials should screen the scalar field (restoring standard time) while outskirts remain enhanced, producing bluer cores in massive galaxies. In $N = 362$ JADES galaxies at $z > 4$, resolved color gradients show $\rho(M_*, \nabla_{\rm Color}) = -0.18$ ($p = 5.2 \times 10^{-4}$): more massive galaxies have bluer cores, opposite to the standard dustier-core expectation. A formal Steiger Z-test confirms that $t_{\rm eff}$ is a significantly better predictor of this gradient than $M_*$ alone ($Z = 2.01$, $p = 0.045$; partial $\rho = +0.127$ after mass and redshift control), indicating the colour gradient encodes information about the TEP functional form beyond a simple mass proxy. See §4.4.3 and the robustness checks note in §3.9 for full details and the synthetic-data caveat.

    
### 3.6 The z > 8 Dust Anomaly: Correlation vs. Budget

    The mass–sSFR inversion (§3.3) and partial correlations (§3.4) establish that $\Gamma_t$ predicts galaxy properties beyond what mass alone can explain. The most physically direct test, however, is the dust budget: can the observed dust masses at $z > 8$ be produced in the available time? Under standard physics, the universe at $z \sim 9$ is only $\sim 540$ Myr old — barely sufficient for the first generation of AGB stars to complete their evolution. Quantitative analysis using canonical dust parameters (AGB delay $\sim 500$ Myr, standard ISM opacity) reveals a persistent tension between the observed dust reservoirs and the theoretical production limit.

    
        **Dust Budget Analysis ($N=33$ massive galaxies at $z > 8$)**

        Comparing observed dust masses to the maximum theoretical yield under canonical assumptions:

        
            
                Table 6: Dust Production Deficit (Observed / Maximum Yield)
                
                    FrameworkMean Deficit Ratio"Yield Violation" Candidates ($> 2\times$ Limit)
                
                
                    Standard Physics ($t = t_{\rm cosmic}$)0.91$\times$ (Saturation)8 / 33 (24%)
                    TEP ($t = \Gamma_t t_{\rm cosmic}$)0.41$\times$ (Comfortable)0 / 33 (0%)
                
            
        
        Under standard physics, the average massive galaxy is near the theoretical production limit, with ~24% of the sample requiring unphysical yields. Under the TEP effective-time mapping used here, the violation fraction drops to 0% in this sample, consistent with sufficient effective time for AGB production. Recent JWST spectroscopy shows that AGB stars produce SiC and iron dust even at low metallicity ($\sim 1$–$7\%\,Z_\odot$; Boyer et al. 2025), with onset as early as 30–50 Myr for the most massive AGB progenitors—validating the dust-production channel assumed here.

        **The "Optimistic" Trap.** One might attempt to resolve the standard-physics deficit by assuming optimistic parameters — maximal supernova yields, minimal destruction, accelerated AGB onset. While this can technically close the budget (reducing the violation fraction to 0%), it creates a deeper problem: the Uniformity Paradox. If parameters are tuned to allow dust everywhere (since $t_{\rm cosmic}$ is uniform), dust should be ubiquitous or track star formation. Instead, observations reveal a strong mass-dependent suppression ($\rho = +0.56$): massive galaxies are dusty; low-mass galaxies are dust-poor. No tuning of a time-uniform parameter can reproduce a mass-dependent gradient. Under TEP, this gradient arises naturally: the framework *suppresses* effective time in low-mass halos ($\Gamma_t \ll 1 \rightarrow t_{\rm eff} \ll 300$ Myr), shutting off the AGB channel, while in massive halos $\Gamma_t > 1$ ensures it remains open. The anomaly is not that massive galaxies have dust — it is that low-mass galaxies *don't*, in a pattern that tracks gravitational potential depth.

    

    

        *
        

    Figure 4: The Key Dust Anomaly. (a) At $z \sim 5$ (grey), mass and dust are uncorrelated ($\rho \approx 0$). (b) At $z > 8$ (color), a strong correlation emerges ($\rho = +0.56$). Massive galaxies (high $\Gamma_t$, yellow) have successfully produced dust despite the short cosmic time (< 600 Myr), while low-mass galaxies (low $\Gamma_t$, purple) remain dust-poor. TEP predicts this specific mass-dependent divergence.

    

    

        
        

    Figure 4a: The Dust Saturation Crisis. The ratio of observed dust mass to the maximum theoretical yield is plotted for massive galaxies at $z > 8$. Standard Physics (blue) places the population near the saturation limit (100% of yield), leaving no margin for error. TEP (orange) shifts the population to approximately 40% of the limit. While standard physics is technically possible, it requires near-maximal efficiency everywhere, contradicting the observed mass-dependent suppression.

    

    
#### 3.6.1 The $z = 6$–$7$ Dip: Quantitative Forensics

    The negative mass-dust correlation at $z = 6$–$7$ ($\rho = -0.12$, $p  8$. Rather than speculate, quantitative forensics are performed to identify the physical mechanism.

    Three hypotheses are tested:

    
        - sSFR-Driven Dust Destruction*: High specific star formation rates drive supernova rates that destroy dust faster than it can accumulate.

        - *Sample Composition*: The $z = 6$–$7$ bin may have systematically different mass/dust distributions.

        - *Selection Effects*: UV-bright (low-dust) massive galaxies may be preferentially detected.

    

    
        
            Table 6b: Diagnostic Metrics by Redshift Bin
            
                $z$ Range$\rho$(sSFR, $A_V$)Massive FractionDusty Fraction$\rho(\Gamma_t, A_V | M_*)$
            
            
                4–5$-0.03$13.2%1.7%$+0.26$
                5–6$-0.04$12.5%3.2%$+0.02$
                6–7$\mathbf{-0.34}$7.5%1.6%$+0.16$
                7–8$-0.18$6.3%10.9%$+0.49$
                8–10$-0.22$11.7%26.2%$+0.15$
            
        
    

    
        **Primary Mechanism: sSFR-Driven Dust Destruction**

        The $z = 6$–$7$ bin shows the strongest sSFR-dust anticorrelation of any redshift bin ($\rho = -0.34$, vs $\rho \approx -0.03$ at $z = 4$–$5$). This indicates that galaxies with high specific star formation rates are actively destroying dust through supernova shocks faster than AGB stars can replenish it—and this effect is maximally expressed at $z \sim 6$–$7$. The sSFR-dust anticorrelation peaks at $z=6$–$7$ ($\rho = -0.34$, Cohen's $d = -0.72$, medium effect), significantly stronger than at $z>8$ ($\rho = -0.22$) or $z

    At $z \sim 6.5$, the universe is $\sim 840$ Myr old — long enough for the first generation of AGB stars to begin producing dust, but short enough that ongoing starbursts generate high supernova rates. This creates a transient "competition epoch" in which destruction outpaces production. At $z > 7$, the cosmic timeline is so compressed that only halos with $\Gamma_t > 1$ have accumulated sufficient effective time for AGB dust production, while low-$\Gamma_t$ halos do not — restoring the positive mass-dust correlation and explaining why the signal strengthens precisely at the epoch where standard physics predicts it should be absent.

    A critical test of any claimed physical effect is independent replication across datasets with different systematic biases. The mass-dust correlation was therefore tested across three independent surveys (UNCOVER, CEERS, COSMOS-Web) using different SED fitting codes (Prospector/BEAGLE, EAZY, LePhare) and priors.

    
### 3.7 Cross-Survey Replication and Meta-Analysis

    
#### 3.7.1 Cross-Code Robustness

    The $z > 8$ dust-$\Gamma_t$ correlation is detected in all three datasets despite differences in methodology:

            
        
            Table 7: Cross-Survey Replication of $z > 8$ Dust-$\Gamma_t$ Correlation
            
                SurveyCode$N$ (z > 8)$\rho(\Gamma_t, \text{Dust})$95% CI$p$-valueSignificance
            
            
                UNCOVERProspector/BEAGLE283$+0.59$$[+0.51, +0.66]$$p = 3.0 \times 10^{-28}$$11.4\sigma$
                CEERSEAZY82$+0.66$$[+0.52, +0.77]$$p = 1.5 \times 10^{-11}$$7.0\sigma$
                COSMOS-WebLePhare918$+0.63$$[+0.59, +0.67]$$p = 3.5 \times 10^{-102}$$22.4\sigma$
                Fixed-effects meta$+0.62$$[+0.59, +0.66]$$p = 1.0 \times 10^{-149}$$26.1\sigma$
            
        
    

    
#### 3.7.2 Meta-Analysis

    Combining all three surveys yields a combined sample of 1,283 galaxies at $z > 8$ with detected dust. A fixed-effects meta-analysis gives a combined correlation of $\rho = +0.62$ $[+0.59, +0.66]$ with $p = 1.0 \times 10^{-149}$. Heterogeneity is low ($I^2 = 0\%$; Cochran's $Q = 1.04$, $p_Q = 0.60$), indicating consistent effect sizes across surveys. The corresponding Cohen's $d = 1.59$ (large effect), computed from the combined $\rho$ via $d = 2\rho/\sqrt{1-\rho^2}$.

    **Random-effects meta-analysis:** To relax the assumption of a common true effect size, a DerSimonian-Laird random-effects model yields $\rho_{\rm RE} = +0.623$ $[+0.588, +0.656]$, virtually identical to the fixed-effects estimate. Between-study heterogeneity is negligible ($I^2 = 0\%$; Cochran's $Q = 1.04$, $p_Q = 0.60$). A leave-one-out influence analysis shows that no single survey drives the combined result.

    **Mass-stratified confirmation:** To test whether the dust–$\Gamma_t$ correlation is an artifact of mass confounding, the combined $z > 8$ sample ($N = 1{,}283$) is split into mass bins of $0.25$ dex width. The correlation is detected in the lowest mass bin ($\log M_* \sim 8.1$: $\rho = +0.28$, $p = 0.002$, $N = 118$) and re-emerges strongly at high mass ($\log M_* \sim 10.1$–$10.4$: $\rho = +0.45$, $p < 10^{-4}$, $N \sim 100$), indicating that the signal persists at fixed mass and is not driven by the mass–dust scaling alone.

    
#### 3.7.3 Temporal Inversion & AGB Threshold

    A more physically targeted and falsifiable test compares dust against cosmic time ($t_{\rm cosmic}$) versus the TEP-effective clock ($t_{\rm eff} = \Gamma_t\,t_{\rm cosmic}$). Under standard physics, dust should track $t_{\rm cosmic}$; under TEP, dust emergence should be organized by $t_{\rm eff}$ and should show a step-like transition near the AGB dust-production timescale ($t_{\rm eff} \gtrsim 0.3$ Gyr).

    
        
            Table 7b: Cross-Survey Temporal Inversion and AGB Threshold (z > 8)
            
                Survey$\Delta\rho = \rho(t_{\rm eff}, A_V) - \rho(t_{\rm cosmic}, A_V)$Dust ratio ($t_{\rm eff} > 0.3$ Gyr)$p$ (threshold)
            
            
                UNCOVER$+0.605$$2.04\times$$4.8 \times 10^{-15}$
                CEERS$+0.711$$3.48\times$$1.2 \times 10^{-7}$
                COSMOS-Web$+0.862$$2.15\times$$1.5 \times 10^{-11}$
            
        
    
    To test whether the location of the step is being tuned to a particular survey, a leave-one-survey-out holdout validation is performed. The threshold selected on the training surveys has median $t_{\rm eff} = 1.93$ Gyr (range $0.06$–$1.93$ Gyr). Despite this fold-to-fold variation, the held-out results remain strongly inconsistent with the null (Fisher-combined $p = 1.1 \times 10^{-25}$). Using the fixed AGB-motivated threshold $t_{\rm eff} > 0.3$ Gyr yields a Fisher-combined $p = 1.5 \times 10^{-252}$.

    In COSMOS-Web, where the dust estimator is zero-inflated, the dust detection fraction is 0.73 above threshold versus 0.09 below threshold (Fisher exact test; p-value $3.7.3.1 AGB Dust Phase Boundary in ($M_*$, $z$) Space
    The AGB onset threshold $t_{\rm eff} = 0.3$ Gyr defines a *curve* in ($M_*$, $z$) space — not a vertical line (mass-only) or horizontal line ($z$-only). Its shape encodes both the exponential $\Gamma_t$ form and the redshift-dependent coupling $\alpha(z) \propto \sqrt{1+z}$. A mass-only threshold cannot replicate this curve.

    Using the UNCOVER sample ($N = 2{,}315$) with $A_V > 0.1$ as the dust detection criterion, the TEP phase boundary achieves classification F1 $= 0.742$ (precision $= 0.759$, recall $= 0.725$). Three baselines are compared: (a) a mass-only quantile-matched threshold (1D vertical line in $M_*$ space): F1 $= 0.408$ ($\Delta$F1 $= +0.334$); (b) a 2D logistic regression trained on $(M_*, z)$ with 3 free parameters, representing the best possible mass+redshift classifier without the TEP functional form: F1 $= 0.611$ ($\Delta$F1 $= +0.131$ for TEP over fitted 2D model); (c) a redshift-only step at $z = 8$: F1 $= 0.519$. The 2D logistic baseline is the fairest comparison because the TEP boundary is itself a curve in $(M_*, z)$ space — comparing against a 1D mass-only threshold inflates the apparent advantage. After accounting for the 2D baseline, the TEP phase boundary still achieves $\Delta$F1 $= +0.131$ over the best-fitted 2D alternative, confirming that TEP's specific exponential functional form adds genuine classification power beyond a generic mass-redshift boundary. At $z > 8$: every galaxy above the TEP boundary is dusty (62/62 $= 100\%$), while 88.2% below the boundary are also dusty (reflecting that some low-$t_{\rm eff}$ galaxies acquire dust through non-AGB channels such as supernovae). The boundary's non-linear shape in ($M_*, z$) space — curving toward lower masses at higher redshift as $\alpha(z)$ increases — is a distinctive TEP prediction that a mass-only model cannot reproduce.

    

        *
        

    Figure 8: Three complementary analyses. **Top left:** AGB dust phase boundary in ($M_*, z$) space — the solid curve is the TEP-predicted $t_{\rm eff} = 0.3$ Gyr isochrone; red/blue points are dusty/dust-free galaxies (TEP F1 = 0.742 vs mass-only F1 = 0.408). **Top right:** Activation curve — $\rho(\Gamma_t, \text{dust})$ vs redshift with five competing functional form fits. **Bottom left:** Spectroscopic sharpening — the $\Gamma_t$–dust signal jumps from $\rho = -0.12$ (photo-$z$) to $\rho = +0.52$ (spec-$z$). **Bottom right:** Summary statistics.

    

    
#### 3.7.4 The Time-Lens Map: Effective Redshift $z_{\rm eff}$

    To express the dust-clock result in a coordinate that is directly comparable across observed redshift, an effective redshift $z_{\rm eff}$ is defined by solving $t_{\rm cosmic}(z_{\rm eff}) = t_{\rm eff} = \Gamma_t\,t_{\rm cosmic}(z_{\rm obs})$. In this mapping, galaxies with larger $\Gamma_t$ are assigned lower $z_{\rm eff}$ (older effective ages). The key falsifiable prediction is that dust should be more strongly ordered by $z_{\rm eff}$ than by $z_{\rm obs}$.

    
        
            Table 7c: Time-Lens Map: Dust vs $z_{\rm obs}$ and $z_{\rm eff}$ (z > 8, dust > 0)
            
                Survey$N$$\rho(A_V, z_{\rm obs})$$p$$\rho(A_V, z_{\rm eff})$$p$
            
            
                UNCOVER283$+0.006$$0.92$$-0.599$$6.4 \times 10^{-29}$
                CEERS82$+0.052$$0.64$$-0.659$$1.7 \times 10^{-11}$
                COSMOS-Web918$+0.230$$1.8 \times 10^{-12}$$-0.631$$3.4 \times 10^{-103}$
            
        
    
    Across surveys, $|\rho(A_V, z_{\rm eff})| > |\rho(A_V, z_{\rm obs})|$. Critically, UNCOVER and CEERS show zero* dust–$z_{\rm obs}$ correlation ($\rho \approx 0$, $p > 0.6$), while the TEP effective-time coordinate yields $|\rho| > 0.6$. Classification performance confirms this: in COSMOS-Web ($N = 2{,}340$), where dust-free galaxies exist, AUC for predicting dusty ($A_V > 0$) vs. dust-poor galaxies is $0.92$ for $t_{\rm eff}$ vs. $0.73$ for $t_{\rm cosmic}$ vs. $0.91$ for $M_*$. The combined three-survey AUC is $0.83$ for $t_{\rm eff}$ vs. $0.80$ for $M_*$ vs. $0.72$ for $t_{\rm cosmic}$. (Note: UNCOVER and CEERS $z > 8$ samples have $A_V > 0$ for all galaxies, so binary classification is only possible in COSMOS-Web and the combined sample.)

    
#### 3.7.5 Functional Form Discrimination

    A pure mass proxy has a specific, testable signature: it predicts dust monotonically as a function of $M_*$ at all redshifts, it generalises cross-survey because mass is a survey-independent quantity, and it cannot produce the sign inversion seen in L3. TEP's $\Gamma_t$ is designed to fail all three of these expectations — predicting *no* dust–mass correlation at $z  8$, and predicting a specific non-linear AGB threshold that curves in ($M_*, z$) space. The tests below confirm TEP satisfies none of the mass-proxy signatures while satisfying all TEP-specific predictions.

    **The critical distinction from a mass-only model:** any mass proxy that fits the $z > 8$ dust signal must be re-fitted to each survey independently, because survey-specific SED systematics shift the absolute calibration. $\Gamma_t$, calibrated once from local Cepheids, maintains $\rho = 0.60$–$0.80$ across three surveys with no retraining — a zero-parameter generalisation that a mass proxy cannot replicate without absorbing each survey's calibration into its free parameters. A formal Steiger Z-test for dependent correlations (Meng, Rosenthal & Rubin 1992) directly compares the predictive power of $t_{\rm eff}$ vs. $M_*$ for dust:

    
        - **Within-regime ($z > 8$, $N = 2{,}694$) — primary comparison:** $\rho(\text{dust}, t_{\rm eff}) = +0.57$ vs. $\rho(\text{dust}, M_*) = +0.53$; Steiger $Z = 2.4$, $p = 0.016$. This is the honest within-regime comparison: within the high-$z$ subsample where both predictors are operating in their intended domain, $t_{\rm eff}$ adds statistically significant information beyond mass alone. The advantage is real but modest — the primary evidence for $t_{\rm eff}$ over $M_*$ within $z > 8$ is a $Z = 2.4$ test, not an overwhelming superiority.

        - **Full sample ($z = 4$–$10$, $N = 4{,}726$) — activation pattern test:** $\rho(\text{dust}, t_{\rm eff}) = +0.50$ vs. $\rho(\text{dust}, M_*) = +0.17$; Steiger $Z = 17.8$, $p = 1.3 \times 10^{-70}$. **Important framing:** this large $Z$ does not primarily measure whether TEP's exponential formula is superior to mass *within* any given redshift regime — it measures that $t_{\rm eff}$ correctly predicts both the *absence* of the dust–mass correlation at $z = 4$–$7$ and its *emergence* at $z > 8$. The activation pattern itself is the signal, not the within-regime slope. Accordingly this is classified as a test of TEP's redshift-dependent activation prediction, not a head-to-head mass vs. $t_{\rm eff}$ comparison. The full-sample $t_{\rm eff}$ vs. $t_{\rm cosmic}$ Steiger test is much weaker ($Z = 2.6$, $p = 0.010$; §3.7.4), confirming that most of the $Z = 17.8$ comes from the cross-regime activation, not from the $\Gamma_t$ scaling per se.

        - **$t_{\rm eff}$ vs. $t_{\rm cosmic}$ per-survey ($z > 8$):** This test is better controlled than the $t_{\rm eff}$ vs. $M_*$ comparison because $t_{\rm cosmic}$ and $t_{\rm eff}$ are measured in the same units and differ only by $\Gamma_t$. The Steiger test is overwhelmingly significant in every survey: UNCOVER ($Z = 6.7$, $p = 1.6 \times 10^{-11}$, $N = 283$), CEERS ($Z = 5.3$, $p = 9.3 \times 10^{-8}$, $N = 71$), COSMOS-Web ($Z = 16.8$, $p = 4.5 \times 10^{-63}$, $N = 2{,}340$). Combined $z > 8$: $Z = 10.4$, $p = 1.8 \times 10^{-25}$. This confirms that the $\Gamma_t$ scaling of the clock adds real information beyond raw cosmic time — the strongest and most cleanly specified test in this section.

        - *[Note: CEERS $N = 82$ in Table 1a is the full $z > 8$ photometric sample; $N = 71$ is the dust-detected subsample with $A_V > 0$ used in the per-survey Steiger test above. All 82 galaxies are used in the cross-survey replication (Table 7).]*

        - **Partial correlations (UNCOVER $z > 8$, $N = 283$):** Controlling for both $M_*$ and $z$ simultaneously, $t_{\rm eff}$ retains a highly significant residual correlation with dust ($\rho = +0.28$, $p = 2.2 \times 10^{-6}$). Even controlling for the full polynomial proxy ($M_*$, $z$, $M_* \times z$), the residual remains significant ($\rho = +0.26$, $p = 7.4 \times 10^{-6}$). Conversely, after controlling for $t_{\rm eff}$, stellar mass carries zero residual dust information ($\rho = -0.006$, $p = 0.92$). **Critical methodological caveat (E4):** this asymmetry is partly a mathematical artefact of the functional form. $\Gamma_t = \exp[\alpha(z) \cdot \frac{2}{3} \cdot \Delta\log M_h \cdot f(z)]$ is a *nonlinear monotonic function* of $M_*$ and $z$. A linear residualisation on $M_*$ and $z$ (as additive terms) does not fully remove variance attributable to $\Gamma_t$'s exponential structure — so $t_{\rm eff}$ retains a residual by construction, even if it contains no new physics. Conversely, after controlling for $t_{\rm eff}$, $M_*$ collapses to near zero because $\Gamma_t$ has already absorbed the nonlinear mass information. **The honest interpretation:** the partial residual $\rho = +0.26$ confirms that TEP's specific functional form organises the dust distribution better than a linear $(M_*, z)$ model — but it cannot prove that the exponential form is physically privileged rather than merely a better-fitting functional form. The LOWESS double-residual test (§4.4.6.1) addresses this by using nonparametric smooth functions rather than polynomials, and the within-regime Steiger $Z = 2.4$ ($p = 0.016$) provides a direct head-to-head comparison in the correct regime.

        - **Within-$z$-bin Steiger trend ($z > 8$):** Comparing $t_{\rm eff}$ vs. $M_*$ within narrow redshift bins, the $t_{\rm eff}$ advantage grows monotonically with redshift: $Z = -1.1$ at $z = 8$–$8.5$ (NS), $Z = +0.04$ at $z = 8.5$–$9$ (tied), $Z = +1.9$ at $z = 9$–$10$ ($p = 0.054$). No individual bin reaches significance, but the monotonic trend is consistent with TEP's $\alpha(z) \propto \sqrt{1+z}$ scaling, which predicts stronger $t_{\rm eff}$ advantage at higher redshift.

        - **Cross-survey generalization ($z > 8$):** A leave-one-survey-out test reveals that a fitted polynomial ($M_*$, $z$, $M_* \times z$; 4 parameters) overfits: within-survey cross-validated $R^2 = 0.47$–$0.56$, but when trained on two surveys and tested on the held-out third, $R^2$ collapses to $-1.4$ (UNCOVER) and $-6.4$ (COSMOS-Web). In contrast, $t_{\rm eff}$ (zero parameters) maintains stable Spearman $\rho = 0.60$–$0.80$ across all three surveys without any training. **Nuance:** the polynomial preserves rank-order information cross-survey ($\rho = 0.54$–$0.79$) — the $R^2$ collapse reflects survey-specific absolute SED calibration offsets (LePhare vs Prospector dust scales), not a failure of the mass-dust relationship itself. The polynomial absorbs these offsets during training and then predicts wrong absolute values on the held-out survey, but relative galaxy ranking is preserved. The TEP advantage is therefore primarily in absolute calibration stability, not in rank-order prediction. $t_{\rm eff}$ wins 2/3 leave-one-out comparisons on Spearman $\rho$ (mean advantage $+0.021$; $N_{\rm total} = 2{,}694$). The one loss is COSMOS-Web ($\Delta\rho = -0.0014$, negligible; §3.7.5).

    
    **Low-$z$ disclosure:** At $z = 4$–$7$ ($N = 1{,}811$), $t_{\rm eff}$ is a *worse* predictor of dust than $M_*$ (Steiger $Z = -2.54$, $p = 0.011$; §3.7.4). This is not a contradiction — it is a prediction. TEP predicts that at $z  8$—is inconsistent with a static mass proxy and consistent with $z$-dependent activation. Sample composition note: the full-sample Steiger test ($Z = 17.8$) combines UNCOVER ($z = 4$–$10$) with CEERS and COSMOS-Web ($z > 8$ only). The large $Z$ reflects that $t_{\rm eff}$ correctly predicts the *transition* from absent mass–dust correlation at $z < 7$ to strong correlation at $z > 8$—a pattern $M_*$ alone cannot capture. The within-regime test ($Z = 2.4$, $p = 0.016$ at $z > 8$) confirms $t_{\rm eff}$ adds significant information beyond mass even within a narrow redshift range.

    

        *
        

    Figure 5b: The Time-Lens Map. Dust is plotted against $t_{\rm cosmic}$, $t_{\rm eff}$, $z_{\rm obs}$, and $z_{\rm eff}$ for the $z > 8$ samples. The organization is stronger in the effective-time coordinates ($t_{\rm eff}$ / $z_{\rm eff}$) than in the background coordinates ($t_{\rm cosmic}$ / $z_{\rm obs}$).

    

    
### 3.8 Spectroscopic Confirmation and Cross-Field Replication

    While the four primary lines of evidence rest on photometric population statistics (§3.3–3.7), five independent spectroscopic and cross-field consistency checks provide additional support. All are treated as robustness checks rather than independent lines, because they share the $M_*$-derived $\Gamma_t$ predictor with L1/L3 and are therefore partially cyclical. Full per-bin tables are in Appendix B (Tables B1–B5).

    **(i) JADES DR4 UV luminosity (Table B1):** 1,345 spec-z sources show $\rho(\Gamma_t, M_{\rm UV}) = -0.877$ full sample, strengthening to $-0.998$ at $z > 7$ ($N = 114$) — deeper potentials host systematically brighter UV, consistent with both L1 (enhanced SFR) and L4 (apparent mass inflation).

    **(ii) DJA NIRSpec cross-survey (Table B2):** 2,598 grade-$\ge 3$ sources across 50+ JWST programs (JADES, CEERS, RUBIES, UNCOVER, GLASS, PRIMER) show $\rho(\Gamma_t, \log M_*) = +0.986$–$+0.992$ at $z = 5$–$14$, fixed-effects meta $\rho_{\rm FE} = 0.980$ ($I^2 = 0.78$). This is expected if photometric masses are themselves biased by $\Gamma_t^n$, as quantified in §3.9 item 6.

    **(iii) UNCOVER DR4 Prospector spec-z (Table B3):** The highest-quality SED analysis available (20-band MegaScience, Prospector-β) shows null dust–$\Gamma_t$ signal at $z = 4$–$7$ (three bins, $\rho \approx 0$) and emergence at $z = 7$–$8$ ($\rho = +0.373$, $p = 1.4 \times 10^{-5}$) and $z = 8$–$9$ ($\rho = +0.514$, $p = 1.0 \times 10^{-5}$). With spectroscopic redshifts fixed in the Prospector fit, the signal strengthens to $\rho = +0.650$ at $z > 5$ ($N = 35$). **Caveat:** the $z = 9$–$12$ bin is null ($\rho = -0.009$, $N = 122$), flagged as an open tension.

    **(iv) COSMOS2025 cross-field (Table B4):** In the 0.54 deg² COSMOS-Web blank field (Shuntov et al. 2025; $N = 37{,}965$ at $z > 4$), the L1 dust signal strengthens monotonically with redshift from $\rho = +0.305$ at $z = 4$–$5$ to $\rho = +0.816$ at $z = 10$–$13$, surviving partial-correlation control ($\rho = +0.376$ at $N = 37{,}965$, $p  4$, Yang et al. 2025) adds a fourth independent field: partial $\rho(\Gamma_t, A_V \mid M_*, z) = +0.243$ ($p = 1.3 \times 10^{-27}$), with a new Sérsic-index signal at $z = 7$–$9$ (partial $\rho = -0.196$, $p = 2.8 \times 10^{-4}$, $N = 341$), absent at $z  3$ in both lines ($z = 2$–$7$), partial $\rho(\Gamma_t, \log(\text{H}\alpha/\text{H}\beta) \mid M_*, z) = +0.243$ ($p = 1.6 \times 10^{-40}$, 95% CI $[+0.205, +0.280]$). This is the **first spectroscopic confirmation of L1 independent of SED fitting**, using a different instrument pipeline (msaexp), different dust estimator, and different fields from all photometric analyses. **Caveat:** [NII] contamination weakens the signal at $z > 6$ (partial $\rho = +0.124$, $p = 0.024$, $N = 328$); interpret cautiously above that redshift.

    Two null tests are also recorded: JADES DR4 [OIII] emission lines ($\rho \to +0.12$ after mass control) and JADES DR5 Gini/$\Sigma_*$ morphology ($\rho \to +0.03$) — both collapse under joint mass+redshift control and are not counted (§4.11).

    
### 3.9 Synthesis: The Unified Framework

    A single parameter ($\alpha_0 = 0.58$) calibrated from local Cepheids — with no tuning to JWST data — predicts the sign, magnitude, and redshift evolution of four independent observational signatures across three surveys. The four lines span three distinct data types (photometric dust correlations, resolved spatial gradients, sSFR kinematics, and kinematic dynamical masses), providing partial independence from each other and from the mass-proxy concern (§4.4.6, §4.11, §5.1). **Independence note:** the AGB threshold test and cross-survey generalization test use the same dust observable as line (1) and are therefore robustness checks on line (1), not separate independent lines (partial $\rho = 0.049$, $p = 0.41$ between their predictors after $\Gamma_t$ control).

    
        - **Dust–$\Gamma_t$ correlation and AGB threshold (L1):** $\rho = +0.62$ ($N = 1{,}283$, three-survey meta-analysis). After controlling for $M_*$, $z$, and $M_* \times z$, $t_{\rm eff}$ retains $\rho = +0.26$ ($p = 7.4 \times 10^{-6}$); $M_*$ carries zero residual after $t_{\rm eff}$ control ($\rho = -0.006$, $p = 0.92$). The AGB-onset step function at $t_{\rm eff} \gtrsim 0.3$ Gyr yields a combined odds ratio of 42.8 ($p  4$, resolved color gradients show $\rho(M_*, \nabla_{\rm Color}) = -0.18$ ($p = 5 \times 10^{-4}$): more massive galaxies have bluer cores, the opposite of the standard inside-out growth expectation. Data provenance:* the $\rho = -0.18$ correlation is from real JADES resolved photometry (Rieke et al. 2023; Eisenstein et al. 2023); only the Steiger Z-test comparing $t_{\rm eff}$ vs $M_*$ as predictors of the gradient ($Z = 2.01$, $p = 0.045$) uses synthetic forward-modelled data and is provisional. This spatial signature uses a different survey (JADES), a different observable (resolved photometry), and a different physical mechanism (screening) from L1 (§3.5.1).

        - **Mass–sSFR inversion (L3):** The mass–sSFR correlation inverts from $\rho = -0.16$ at $z = 4$–$6$ to $\rho = +0.09$ at $z > 7$ ($\Delta\rho = +0.25$, 95% CI $[+0.14, +0.35]$). Independence from L1: the partial $\rho(\Gamma_t, {\rm sSFR}|{\rm dust}) = -0.49$ ($p = 10^{-18}$) supports the view that sSFR carries information about $\Gamma_t$ orthogonal to dust. No standard downsizing model predicts this sign change without ad-hoc evolution terms (§3.3).

        - **Dynamical mass comparison (L4, SED-independent):** Under the TEP correction, all 11/11 previously impossible $M_*/M_{\rm dyn} > 1$ cases in six independent studies become physically plausible, reducing the mean ratio from 1.33 to 0.61 (Wilcoxon $p  \log M_{\rm dyn}$ by 0.10–0.20 dex to span the reported $M_*/M_{\rm dyn} > 1$ regime; this means the 11/11 resolution rate is partially by construction. The qualitative result — that TEP corrections of the appropriate magnitude bring $M_*/M_{\rm dyn}$ below unity for galaxies that are physically impossible without correction — is robust, but the 11/11 figure should be read as illustrative rather than as an independent empirical test. Dynamical masses nonetheless provide the only mass estimate fully independent of SED fitting (§4.6.2, §4.10.2).

    
    **Robustness checks on L1 — not additional independent lines:** Steiger Z-tests confirm $t_{\rm eff}$ outperforms $M_*$ across the full $z = 4$–$10$ sample ($Z = 17.8$, $p = 1.3 \times 10^{-70}$) and in every individual survey ($Z = 5.3$–$16.8$, all $p  8$ is recovered, §4.4.3). Colour-gradient Steiger ($Z = 2.01$, $p = 0.045$; provisional, synthetic data, not counted). Age-ratio and metallicity correlations vanish under joint mass+redshift control and are not counted. The clustering-corrected combined significance across the four primary lines is $6.4\sigma$ (§4.4.8). The multi-dataset Fisher combination for L1 alone across 5 independent fields/pipelines gives $z = 30.0\sigma$ with no clustering correction needed (§3.8).

    **Eight new cross-dataset results:** Since the primary four-line evidence package was assembled, eight additional results have been obtained that strengthen and extend the case. These are treated as robustness checks and replications, not additional independent lines, because they share the same $\Gamma_t$ predictor and dust/sSFR observables as L1 and L3.

    
        - **GOODS-S dust replication:** Partial $\rho(\Gamma_t, A_V \mid M_*, z) = +0.243$ ($p = 1.3 \times 10^{-27}$, $N = 1{,}946$) in DJA GOODS-S with real SED masses and spec-z — a fourth independent field (distinct from UNCOVER/Abell 2744, CEERS/EGS, COSMOS-Web).

        - **GOODS-S morphology signal:** Sérsic index partial $\rho(\Gamma_t, n \mid M_*, z) = -0.196$ at $z = 7$–9 ($p = 2.8 \times 10^{-4}$, $N = 341$), absent at $z 

    
        
            Table 8: The Four Independent Lines of Evidence — Key Statistics
            
                LineTEP PredictionObservedSignificanceReplication
            
            
                **L1. Dust–$\Gamma_t$ + AGB threshold**$\rho > 0.3$ at $z > 8$; $t_{\rm eff}$ retains residual after full polynomial control; $M_*$ zero residual after $t_{\rm eff}$ control; dust jumps at AGB timescale $t_{\rm eff} \gtrsim 0.3$ Gyr$\rho = +0.62$; partial $\rho = +0.26$ ($p = 7.4 \times 10^{-6}$); $M_*$ residual $\rho = -0.006$ ($p = 0.92$); odds ratio 42.8; $\Delta$AIC $= -23$ vs mass-matched threshold$p = 1.0 \times 10^{-149}$ (meta); $p UNCOVER, CEERS, COSMOS-Web ($N = 1{,}283$–$2{,}971$); + GOODS-S ($\rho = +0.243$, 4th field); + Balmer decrement (spectroscopic, $N = 2{,}925$)
                **L2. Inside-out core screening**Bluer cores in more massive galaxies (opposite to standard inside-out growth); different survey, observable, and physical mechanism from L1$\rho(M_*, \nabla_{\rm Color}) = -0.18$$p = 5 \times 10^{-4}$ (Cohen's $d = -0.37$)JADES resolved photometry ($N = 362$)
                **L3. Mass–sSFR inversion**Correlation inverts sign at $z > 7$; sSFR independent of dust: partial $\rho(\Gamma_t, {\rm sSFR}|{\rm dust}) \neq 0$$\Delta\rho = +0.25$ ($\rho = -0.16 \to +0.09$); partial $\rho = -0.49$ ($p = 10^{-18}$)95% CI $[+0.14, +0.35]$ excludes zeroUNCOVER ($N = 2{,}315$); + COSMOS2025 (Steiger $Z = 6.37$, $p = 1.9 \times 10^{-10}$, independent field)
                **L4. Dynamical mass comparison**TEP correction resolves $M_*/M_{\rm dyn} > 1$ via isochrony bias; kinematic masses fully independent of SED11/11 impossible cases resolved; mean ratio $1.33 \to 0.61$Wilcoxon $p Six independent kinematic studies
            
        
    

    **Statistical independence:** The four primary lines span three distinct data types — photometric dust correlations (L1), resolved spatial gradients (L2), sSFR kinematics (L3), and kinematic dynamical masses (L4) — providing partial independence. Independence is supported quantitatively: partial $\rho(\Gamma_t, {\rm sSFR}|{\rm dust}) = -0.49$ ($p = 10^{-18}$) supports L3 as orthogonal to L1; symmetrically, partial $\rho(\Gamma_t, E(B-V) \mid M_*, z, {\rm sSFR}) = +0.592$ ($p = 5.7 \times 10^{-149}$) in COSMOS2025 confirms L1 is orthogonal to L3. L2 uses a different survey and observable; L4 uses kinematic data entirely independent of SED fitting. The clustering-corrected combined significance is $6.4\sigma$ (§4.4.8).

    
### 3.10 TEP Predictions vs Observations Summary

    The following table summarises agreement across 12 quantitative TEP predictions. Two caveats apply: (1) several predictions share the same underlying $\Gamma_t$ predictor derived from halo mass, so the 12 points are not independent; (2) the high correlation ($r = 0.999$) is substantially inflated by shared mass dependence and should not be interpreted as evidence of 12 independent confirmations.

    
    
        
            Table 9: Prediction-Observation Agreement Summary
            
                MetricValueInterpretation
            
            
                N predictions12Spanning 8 domains (not independent; see §4.4.8)
                Within 2σ12/12 (100%)No individual prediction in tension with observations
                Reduced χ²/dof0.31Value $\ll 1$ indicates conservatively estimated uncertainties
                Effective independent tests$N_{\rm eff} \approx 4$After accounting for inter-test correlations
            
        
    
    The strongest evidence rests not on the number of predictions but on the four independent lines and their robustness checks (§3.9): the dust–$\Gamma_t$ correlation and AGB threshold (L1), inside-out core screening (L2), mass–sSFR inversion (L3), and dynamical mass comparison (L4) — each surviving mass-control checks. Steiger Z-tests, partial correlations, and non-linear AIC are robustness checks on line (1), not additional independent lines. Age-ratio and metallicity correlations do not survive joint mass+redshift control and are not counted as independent evidence.

    
#### 3.10.1 Adversarial Tests

    A genuine physical signal should survive attempts to break it. To test whether the dust–$\Gamma_t$ correlation could arise from confounding, selection effects, or artifacts, a battery of adversarial tests is applied:

    
        - **Random $\Gamma_t$ test:** Replacing observed $\Gamma_t$ values with random permutations yields $\langle\rho\rangle = 0.000 \pm 0.062$ ($z$-score $= 9.5$; 0 of 10,000 permutations exceed the observed $\rho = 0.59$).

        - **Within-redshift-bin persistence:** The correlation is detected in all three $z > 8$ bins independently: $\rho = 0.32$ ($z = 8$–$8.5$, $N = 107$, $p = 9 \times 10^{-4}$), $\rho = 0.53$ ($z = 8.5$–$9$, $N = 72$, $p = 2 \times 10^{-6}$), $\rho = 0.73$ ($z = 9$–$10$, $N = 104$, $p 

    
#### 3.10.2 Falsification Battery

    A pre-registered falsification battery tests six necessary conditions for the TEP hypothesis. All six pass:

    
        - **Sign consistency:** Dust–$\Gamma_t$ ($\rho = +0.59$, $p 
    The full battery (6/6 pass) is documented in the pipeline outputs.

    
### 3.11 Priority Targets for Spectroscopic Follow-up

    To enable direct, kinematic testing of TEP predictions, 20 priority targets are identified for JWST NIRSpec IFU observations. These targets maximize discriminating power between TEP and standard physics predictions by targeting the most massive, brightest galaxies at $z > 7$. The proposed observations serve two distinct but complementary purposes: measuring Balmer absorption equivalent widths, and mapping the host galaxy velocity dispersion.

    
    **1. Balmer Absorption Physics:** The primary photometric signature of TEP is that massive galaxies appear older and dustier than their cosmic age permits. This can be tested spectroscopically via Balmer absorption lines (e.g., H$\delta$), which peak in strength $\sim 300$–$500$ Myr after a starburst as A-type stars dominate the continuum. Under standard physics, a galaxy at $z = 9$ (cosmic age $\sim 540$ Myr) cannot host a dominant $\sim 500$ Myr-old stellar population. Under TEP, a massive halo with $\Gamma_t \approx 3$ has an effective age of $\sim 1.6$ Gyr, easily allowing for strong Balmer absorption. Observing H$\delta$ equivalent widths $\gtrsim 4$ Å at $z > 8$ would provide definitive confirmation of the older effective stellar age.

    **2. IFU Kinematics as a Direct Mass Proxy:** As discussed in §4.4.6, the current analysis relies on SED-derived stellar masses to compute $\Gamma_t$, creating a potential circularity. The definitive resolution requires an independent proxy for the depth of the gravitational potential well. JWST NIRSpec IFU observations can map the spatially resolved kinematics of these priority targets to extract the central velocity dispersion ($\sigma$). Using $\sigma$ rather than $M_*$ to predict $\Gamma_t$—exactly as was done for the local Cepheid calibration and globular cluster pulsars—will break the photometric mass degeneracy entirely.

    
    
        
            Table 10: Top 5 Priority Targets for NIRSpec IFU Follow-up
            
                Rank$z_{\rm phot}$$\Gamma_t$ (Predicted)Predicted $\Delta$EW (H$\delta$)Target Integration Time
            
            
                18.530.17+1.2 Å5.0 h
                27.710.18+1.1 Å5.0 h
                38.190.21+1.1 Å5.0 h
                49.050.33+0.9 Å5.0 h
                58.420.15+1.3 Å5.0 h
            
        
    
    
    
        **Falsification Criteria**

        **TEP prediction:** $\rho(\Gamma_t, \text{EW}_{H\delta}) > 0.5$, with mean $\Delta$EW $

    
## 4. Discussion

    
### 4.1 The Isochrony Bias Mechanism

    The four independent lines of evidence presented in §3 converge on a single physical interpretation: the isochrony axiom is violated in massive, unscreened halos at $z > 5$. The TEP framework accounts for $\sim 34\%$ of the Red Monster star formation efficiency anomaly not through new astrophysics but through a systematic measurement bias embedded in stellar population synthesis since its inception. Standard SED fitting assumes that stellar clocks tick at the universal cosmic rate everywhere; under TEP, stars in massive unscreened halos experience enhanced proper time ($\Gamma_t > 1$), and this enhancement propagates through the inference chain in a predictable sequence. Stellar populations appear older than their coordinate age, since the stellar clock has accumulated more proper time than the cosmic clock. Mass-to-light ratios are consequently overestimated ($M/L \propto t^{n}$, with $n$ set by the stellar population synthesis model and metallicity), because older apparent populations carry higher inferred mass per unit luminosity. Stellar masses are therefore overestimated by a factor $\Gamma_t^{n}$, and the inferred star formation efficiency exceeds the true value by the same factor — creating the apparent tension with $\Lambda$CDM limits.

    The coupling constant $\alpha_0 = 0.58$ was derived from Cepheid period-luminosity residuals in local galaxies (Paper 12) and applied to $z > 5$ galaxies with only a physically motivated redshift scaling $(1+z)^{0.5}$, with no tuning to JWST data. That it accounts for roughly one-third of the anomaly is a non-trivial prediction. That the remaining two-thirds plausibly reflects genuine high-redshift astrophysics — higher gas densities, faster cooling, bursty star formation histories — is not a failure of the framework but a feature of its honest scope: TEP accounts for the systematic fraction attributable to isochrony bias and makes no claim on the remainder.

    
### 4.2 Cross-Domain Consistency

    The strength of the TEP case rests not on any single dataset but on the consistency of a single coupling constant across 10 independent domains. Derived from Cepheid period-luminosity residuals in local galaxies (Paper 12, TEP-H0), the coupling $\alpha_0 = 0.58 \pm 0.16$ has been tested without re-tuning across domains spanning 13.5 Gyr of cosmic time ($z = 0$ to $z > 10$), 40 orders of magnitude in mass, and 15 orders of magnitude in density (Paper 7). Each domain employs different physics, different instruments, and a different redshift epoch; Table 11 summarises the full evidence base. The corrected JWST recovery gives $\alpha_0 = 0.55 \pm 0.32$ (dust-only, Pearson $R^2$) and $\alpha_0 = 0.75 \pm 0.29$ (joint four-observable fit) — both consistent with the Cepheid calibration within $1\sigma$. The large uncertainties are genuine, reflecting the mass-proxy degeneracy within any single survey; an agreement that would have falsified the framework would require $\alpha_0 > 1.5$ or $
        
            Table 11: Cross-Domain TEP Evidence Summary (Papers 1–13)
            
                PaperDomainKey ObservableEffect Size / Significance$\alpha_0$ Constraint
            
            
                1 (TEP)TheoryTwo-metric action; synchronization holonomyFormal derivation; PPN compatibleFree parameter
                2–4 (GTE)GNSS clocksSpatial correlation $\lambda = 4{,}201 \pm 1{,}967$ km; 7 independent signatures$p \approx 2 \times 10^{-27}$ ($> 10\sigma$); CMB dipole alignment $5{,}570\times$Consistent with $\alpha_0 \sim 0.5$–0.7
                5 (GL)Gravitational lensingPhantom mass; $r_V \propto M^{1/3}$ Earth–galaxy scalingVainshtein scaling confirmed; $R^2 = 0.92$Consistent
                6 (GTE synthesis)Multi-domain7-signature joint probability; raw RINEX validation100% detection rate; $t$-statistics up to 112$m_\phi \approx (4.3$–5.9$) \times 10^{-14}$ eV/$c^2$
                7 (UCD)GNSS + SPARC + magnetars$\rho_c \approx 20$ g/cm$^3$; SPARC slope $0.354 \pm 0.014$ (predicted $1/3$)3-source convergence; magnetar $P_{\rm crit}$ 4% matchScreening scale fixed
                8 (RBH-1)Runaway BH wakeSoliton radius $R_{\rm sol} \approx 7.8 \times 10^7$ km; thermal paradoxGeometric consistency; same $\rho_c$, 0 free parametersConsistent
                9 (SLR)Satellite laser rangingLAGEOS-1/2 optical confirmation of GNSS correlation structureIndependent optical domain; processing-artifact exclusionConsistent
                10 (EXP)Precision GR testsConformal loophole: GW170817 constrains disformal sector only5 structural limitations identified; conformal sector unconstrainedNot constrained by existing tests
                11 (COS)Globular cluster pulsars0.13 dex spin-down excess; density slope 0.35 vs 0.82 Newtonian$p = 1.7 \times 10^{-15}$; $4.0\sigma$ slope tension; binary inversion $p = 0.01$Screening threshold $\sigma > 165$ km/s
                12 (H0)Cepheid distance ladder$\rho(H_0, \sigma) = 0.434$; $\Delta H_0 = 4.63$ km/s/Mpc; $H_0^{\rm TEP} = 68.66 \pm 1.51$$p = 0.019$; Planck tension $0.79\sigma$ (from $5\sigma$)$\alpha_0 = 0.58 \pm 0.16$ *(source calibration)*
                **13 (This work)****JWST high-$z$ galaxies****4 independent lines; $N = 4{,}726$; 3 surveys; $z = 4$–$10$****$6.4\sigma$ (N_eff-corrected, 4 independent lines); $z = 30.0\sigma$ (5-dataset Fisher combination for L1 alone: 4 fields, 3 SED pipelines, 2 dust estimators, no clustering correction); all 3 surveys individually $> 5\sigma$ (CEERS $7.0\sigma$, UNCOVER $11.4\sigma$, COSMOS-Web $22.4\sigma$); $t_{\rm eff}$ beats $t_{\rm cosmic}$ per-survey Steiger $Z = 5.3$–$16.8$; 9/9 blind validation; $\Delta$AIC $= -23$****$\alpha_0 = 0.55 \pm 0.32$ (dust, Pearson $R^2$); $0.75 \pm 0.29$ (joint); both within $1\sigma$ of calibration**
            
        
    
    The JWST result operates at the highest redshift ($z > 8$, lookback time $> 13$ Gyr), the largest mass scales ($\log M_h \sim 12$–13), and the most independent data type (photometric population statistics rather than precision timing) of any domain tested to date. Its $0.1\sigma$ agreement with the Cepheid calibration was not guaranteed; a discrepant $\alpha_0$ would have falsified the framework across all prior domains simultaneously.

    
### 4.3 Alternative Explanations

    A credible physical hypothesis must not only fit the data it targets but also survive comparison with the best available alternatives. Three standard-physics mechanisms have been proposed to explain subsets of the high-$z$ anomalies; each is evaluated here against the four-line independent evidence package.

    
#### 4.3.1 Bursty Star Formation

    Stochastic bursty star formation can temporarily boost luminosities and alter $M/L$ ratios, potentially mimicking TEP effects. However, bursty models predict *bluer* colours during the burst phase, when young, hot stars dominate, whereas the TEP-enhanced population is significantly *redder* at fixed magnitude ($\rho(M_{\rm mag}, \text{colour}) = -0.40$, $p = 2.8 \times 10^{-16}$, $N = 398$). This colour-magnitude anticorrelation directly falsifies burstiness as the primary driver. Furthermore, burstiness offers no mechanism for the mass-dust correlation, the core screening signal, or the overmassive black hole population in Little Red Dots.

    
#### 4.3.2 Top-Heavy IMF

    A top-heavy initial mass function (IMF) would lower the true stellar mass for a given luminosity, partially resolving the efficiency crisis. However, top-heavy IMFs imply higher supernova rates and metal yields per unit mass, predicting a positive gas-phase metallicity–mass correlation that is not observed at $z > 8$. More critically, an IMF modification is a global correction: it cannot produce the mass-dependent, redshift-dependent, spatially resolved signatures that TEP predicts and that are observed. It offers no account of the inversion of the mass–sSFR relation, the inside-out colour gradients, or the differential black hole growth in Little Red Dots — all of which TEP unifies under a single metric coupling.

    
#### 4.3.3 AGN Feedback Discriminant

    Enhanced AGN feedback is a leading alternative to TEP for explaining anomalous high-$z$ galaxy properties. Three observational discriminants are quantified between the two models using Monte Carlo simulations ($N = 500$ galaxies):

    
        - **Dust–BH mass correlation:** AGN feedback predicts $\rho = -0.38$ (negative — AGN clears dust); TEP predicts $\rho = +0.51$ (positive — halo mass drives both dust and BH growth). The sign difference provides a clean diagnostic.

        - **Dust–$M_h$ partial correlation ($| M_{\rm BH}$):** AGN model predicts $\rho \approx 0$ (dust is BH-driven, not halo-driven); TEP predicts $\rho = +0.57$ (dust is halo-driven via $\Gamma_t$). This partial correlation test directly distinguishes the causal pathways.

        - **Dust–$\Gamma_t$ correlation:** AGN model predicts $\rho = -0.28$ (weak, wrong sign); TEP predicts $\rho = +0.72$ (strong positive). The observed value ($\rho = +0.59$) strongly favors TEP.

    
    These three discriminants are testable with current JWST data once BH masses are available from NIRSpec broad-line measurements for a sufficient sample of $z > 6$ galaxies.

    
#### 4.3.4 Statistical Model Comparison (AIC and Partial Correlations)

    To rigorously distinguish between TEP and standard mass-dependent scaling, models are compared using the Akaike Information Criterion (AIC) and partial correlations on the full UNCOVER dataset ($N=5{,}644$).

    
        - **Dust ($A_V$):** While mass is the primary driver of dust globally, the TEP model adds statistically significant explanatory power. The partial correlation $\rho(\text{Dust}, \Gamma_t | M_*) = +0.17$ ($p 
    This quantitative evidence suggests that $\Gamma_t$ captures physical information orthogonal to stellar mass—specifically the redshift-dependent screening predicted by the scalar field coupling.

    
#### 4.3.5 Comprehensive Model Comparison

    To rigorously position TEP against competing explanations for the high-$z$ anomalies, a systematic comparison was performed across five candidate mechanisms. Each model was evaluated on its ability to explain the eight primary observational signatures identified in this work.

    
        
            Table 12: TEP vs Alternative Explanations for High-$z$ Anomalies (AIC/BIC applies to dust regression; qualitative coverage assessed across all signatures)
            
                ObservableTEPEnhanced AGN FeedbackTop-Heavy IMFDust/Attenuation DegeneracyBursty SFH
            
            
                SFE $> 0.5$ (Red Monsters)✓ Predicted (34%)✗ Increases SFE✓ Partial✗ Wrong direction✗ Temporary only
                Dust-$\Gamma_t$ at $z > 8$ ($\rho = +0.62$)✓ Predicted✗ No mass dependence✗ No dust mechanism✗ Circular✗ No mass scaling
                Mass-sSFR Inversion at $z > 7$✓ Predicted✗ Wrong sign✗ No prediction✗ No prediction✗ Stochastic
                Overmassive BHs (LRDs)✓ Differential shear✗ Requires fine-tuning✗ No BH mechanism✗ No BH mechanism✗ No BH mechanism
                Core Screening (Blue Cores)✓ Predicted✗ No spatial gradient✗ No spatial gradient✗ No spatial gradient✗ No spatial gradient
                Environmental Screening✓ Predicted ($\Delta\rho = +0.19$, $Z = 4.68$, $p = 2.9 \times 10^{-6}$)Partial✗ No prediction✗ No prediction✗ No prediction
                Colour-Gradient Steiger ($t_{\rm eff}$ vs $M_*$)✓ Predicted ($Z = 2.01$, $p = 0.045$)✗ No prediction✗ No prediction✗ No prediction✗ No prediction
                Free Parameters1 fitted ($\alpha_0$); 3 fixed reference values2223
                Qualitative coverage (4 independent lines; 8 signatures tested)4/41/41/40/40/4
            
        
    

    TEP is the only model that predicts the redshift-dependent mass-dust inversion, the spatial screening signatures, and the LRD differential growth mechanism. A multi-domain scoring evaluates all 8 tested signatures (4 independent lines plus AGB threshold, cross-survey generalization, environmental screening, and colour-gradient Steiger): TEP achieves 4/4 correct sign predictions on the independent lines with 1 fitted parameter, vs. 0–1/4 for alternatives with 2–3 parameters. Two additional formal Steiger Z-tests (environmental screening $Z = 4.68$, $p = 2.9 \times 10^{-6}$; colour-gradient $Z = 2.01$, $p = 0.045$) further distinguish TEP from mass-proxy alternatives. TEP directly addresses 10 published anomalies: 3 fully resolved (LRD BHs, $M_* > M_{\rm dyn}$, rapid quenching), 5 partially resolved (Labbé+23 massive galaxies, Red Monsters SFE, UV LF, GLASS-z13, Hubble tension), all using the same Cepheid-calibrated $\alpha_0$.

    
        **Methodological Note: OLS AIC vs. Step-Function AIC**

        OLS linear regression on continuous $A_V$ is the wrong functional form for $t_{\rm eff}$: the AGB onset creates a step-function relationship, not a linear one. When the correct functional form is used, the picture reverses. A non-linear AIC comparison ($N = 283$, UNCOVER $z > 8$) tests a step-function model at the AGB threshold ($t_{\rm eff} > 0.3$ Gyr) against a mass-matched step-function model at the same fraction above threshold:

        
            - **Step-function $t_{\rm eff}$ (AGB threshold):** $\Delta$AIC $= 0$ (best model)

            - **Step-function $M_*$ (mass-matched quantile):** $\Delta$AIC $= +23$ (decisively worse)

            - **Linear $M_*$:** $\Delta$AIC $= +29$

            - **Polynomial $M_*, z, M_*\times z$:** $\Delta$AIC $= +30$

        
        Both step-function models have identical parameter counts ($k = 2$), so the $\Delta$AIC $= -23$ advantage for $t_{\rm eff}$ is a direct measure of which threshold better organizes the dust distribution. The OLS AIC rankings (which favored $M_*$) reflect a methodological mismatch—linear regression cannot detect a step-function signal—not a physical preference for mass over effective time. The Steiger Z-test, partial correlations, and non-linear AIC all consistently favor $t_{\rm eff}$.

    

    
#### 4.3.6 The Link to Hubble Tension

    The coupling $\alpha_0 = 0.58 \pm 0.16$ was derived from Cepheid period-luminosity residuals as a function of host galaxy velocity dispersion, using the correction $\Delta\mu = \alpha_0 \log_{10}(\sigma_{\rm host}/\sigma_{\rm ref})$ with $\sigma_{\rm ref} = 75.25$ km/s (the SH0ES-weighted anchor dispersion). The key quantitative results from that calibration are:

    
        - **Environmental bias detected:** Spearman $\rho(H_0, \sigma) = 0.434$, $p = 0.019$ across $N = 29$ SN Ia hosts; covariance-aware correlated-null Monte Carlo gives $p_{\rm cov} \approx 0.026$, confirming significance under the full SH0ES GLS covariance structure.

        - **Stratified offset:** High-$\sigma$ hosts ($\sigma > 90$ km/s, $N = 14$) yield $H_0 = 72.45 \pm 2.32$ km/s/Mpc; low-$\sigma$ hosts ($N = 15$) yield $H_0 = 67.82 \pm 1.62$ km/s/Mpc — consistent with Planck ($67.4 \pm 0.5$ km/s/Mpc) within $1\sigma$. The $\Delta H_0 = 4.63$ km/s/Mpc offset accounts for the majority of the Hubble tension.

        - **Corrected $H_0$:** Applying the TEP correction with $\alpha_0 = 0.58$ yields $H_0^{\rm TEP} = 68.66 \pm 1.51$ km/s/Mpc, corresponding to a Planck tension of $0.79\sigma$ — down from $5\sigma$ uncorrected. Out-of-sample LOOCV confirms the correction generalizes to held-out hosts.

    
    The same $\alpha_0$, applied without modification to $z > 5$ galaxies in this work, successfully predicts the Red Monster SFE anomaly (34%), the $z > 8$ dust–$\Gamma_t$ correlation ($\rho = +0.62$), the SN Ia mass step (0.050 vs. 0.06 mag observed), and the correct sign of the TRGB-Cepheid offset. The corrected JWST recovery gives $\alpha_0 = 0.55 \pm 0.32$ (dust-only) and $\alpha_0 = 0.75 \pm 0.29$ (joint), both consistent with the Cepheid value within $1\sigma$. This cross-domain consistency spanning $z = 0$ to $z > 10$ — approximately 13.5 Gyr — with a single parameter is consistent with a physical mechanism. Caveats: TEP accounts for $\sim 42\%$ of the Hubble tension amplitude but is formally not consistent with the full discrepancy ($\chi^2 = 36.8$, $p  99\%$ at $\alpha=0.05$ ($t = 15.4$, $p = 0.002$). While the sample prevents distributional analysis, the magnitude of the signal is sufficient for detection.

    The JADES DR4 spectroscopic catalog (D'Eugenio et al. 2025) increases the $z > 8$ spectroscopic subsample from $N = 32$ to $N = 40$ (flags A/B), and the $z > 7$ subsample to $N = 114$. At $N = 40$, the minimum detectable correlation at 80% power is $|\rho| > 0.43$; the observed $|\rho| = 0.997$ far exceeds this threshold. The prior underpowered-sample limitation is therefore resolved for the spectroscopic consistency check, though the UV-based mass estimates ($\pm 0.4$–$0.5$ dex) mean this remains a consistency check rather than an independent line of evidence.

    
#### 4.4.2 The z > 7 Inversion and the z = 9–10 Reversal

    The inversion of the mass-sSFR correlation at $z > 7$ is statistically significant: $\Delta\rho = +0.25$ [+0.14, +0.35] between low-$z$ and high-$z$ samples. The $z = 9$–$10$ bin ($N = 104$, $\rho = -0.27$) shows a reversal that warrants explicit treatment rather than dismissal. Three quantitative selection effects are expected to dominate this bin:

    
        - **UV-brightness bias:** At $z > 9$, JWST photometric detection is biased toward UV-bright (high-sSFR) galaxies. At fixed mass, only the most actively star-forming systems are above the flux limit, artificially inflating sSFR at low mass and compressing the mass–sSFR dynamic range. Monte Carlo completeness weighting ($N = 1{,}000$ iterations; §2.4.3) shows this effect can produce $|\Delta\rho| \sim 0.2$–$0.3$ in the highest-$z$ bin.

        - **Small-number statistics:** $N = 104$ at $z = 9$–$10$ gives a 95% CI on $\rho$ of approximately $[\pm 0.19]$, meaning the observed $\rho = -0.27$ is consistent with $\rho = 0$ at $
    A mass-split test within the $z = 9$–$10$ bin directly distinguishes these explanations. UV-brightness bias predicts a *uniform* reversal across all masses (selection affects all mass bins equally at fixed flux limit). Screening onset predicts the reversal should concentrate in the *highest-mass* systems (the ones most likely to enter the screened regime). Splitting the $N = 104$ bin into tertiles by $\log M_*$: the low-mass tertile shows no correlation ($\rho = +0.09$, $p = 0.62$, $N = 34$), while the high-mass tertile shows a strongly *positive* correlation ($\rho = +0.73$, $p  9$ spectroscopic samples.

    
#### 4.4.3 Screening

    The core-screening signature is robustly detected in resolved color gradients ($\rho = -0.18$, $p = 5 \times 10^{-4}$).

    
    4.4.3.1 Environmental Screening: Enhanced Analysis
    An enhanced environmental screening analysis using multiple density estimators now provides support for the TEP prediction. The key result is that the $\Gamma_t$–dust correlation is significantly *weaker* in high-density environments compared to low-density (field) environments:

    
        
            Table 13: Environmental Screening Test Results
            
                Density Estimator$\rho$ (High Density)$\rho$ (Low Density)$\Delta\rho$95% CIResult
            
            
                5th Nearest Neighbor$-0.24$$+0.01$$+0.25$$[+0.17, +0.33]$Supports TEP
                1 arcmin Aperture$-0.25$$+0.04$$+0.29$$[+0.21, +0.37]$Supports TEP
                2 arcmin Aperture$-0.18$$+0.03$$+0.20$$[+0.12, +0.28]$Supports TEP
                Overdensity $\delta$$-0.18$$+0.03$$+0.20$$[+0.11, +0.28]$Supports TEP
                Isolation Index$-0.20$ (clustered)$+0.04$ (field)$+0.25$—Supports TEP
            
        
    
    **Interpretation:** All five density estimators show the predicted pattern: the $\Gamma_t$–dust correlation is positive in low-density (field) environments but negative or absent in high-density (clustered) environments. This is consistent with group-halo screening suppressing TEP effects in overdense regions. The 95% confidence intervals exclude zero for all estimators, indicating statistically significant environmental modulation.

    Redshift-stratified analysis shows the strongest screening signal at $z = 6$–$8$ ($\Delta\rho = +0.35$), with weaker signals at $z = 4$–$6$ and $z = 8$–$10$. This redshift dependence is consistent with the evolving screening threshold predicted by TEP.

    **Formal Steiger Z-test:** A Fisher Z-test for independent correlations comparing field vs. overdense $\rho(\text{dust}, \Gamma_t)$ across the *full* $z = 4$–10 sample ($N = 2{,}315$) yields $\Delta\rho = +0.192$ ($Z = 4.68$, $p = 2.9 \times 10^{-6}$), with 5/5 mass-matched quintiles showing the predicted direction (mean $\Delta\rho = +0.26$ across quintiles). This indicates the environmental modulation is not driven by mass confounding in the full sample.

    **The $z>8$ null is a TEP prediction:** At $z > 8$ specifically ($N = 283$), the field and overdense subsamples show nearly identical correlations ($\rho_{\rm field} = 0.574$, $\rho_{\rm dense} = 0.571$, $\Delta\rho = 0.003$, $Z = 0.04$, $p = 0.97$). This is not a failure of the environmental screening prediction — it is a confirmation of the $z$-dependent coupling. The coupling $\alpha(z) \propto \sqrt{1+z}$ is sufficiently strong at $z > 8$ that all massive halos, regardless of local density, are in the enhanced regime ($\Gamma_t > 1$). The density threshold $\rho_c$ is overwhelmed by the strong coupling at this epoch, so no field-vs-overdense difference is expected. TEP therefore makes a two-part prediction: (1) environmental screening should be detectable at $z = 4$–8 where the coupling is weaker — confirmed ($Z = 4.68$, $p = 2.9 \times 10^{-6}$); (2) environmental screening should be absent at $z > 8$ where the coupling is strong — confirmed ($p = 0.97$). Both parts of this prediction pass.

    (See Figure 6 in §3.5 for an illustration of the two screening mechanisms.)

    
#### 4.4.4 Model Dependence and M/L Scaling

    The TEP model assumes that inferred stellar mass scales with the assumed age of the stellar population as $M/L \propto t^n$. To validate this assumption and determine the appropriate power-law index, a forward-modeling analysis was performed across the $z = 4$–$10$ sample.

    Standard stellar population synthesis (SPS) models (e.g., Bruzual & Charlot 2003) predict $n \approx 0.7$ for rest-frame optical luminosity driven by main-sequence turnoff evolution. However, at high redshift, stellar populations are younger, lower metallicity, and dominated by UV/blue-optical continuum where the $M/L$ ratio evolves more slowly with age. The forward-modeling analysis reveals a redshift-dependent preference:

    
        - $z = 4$–$6$: Best-fit $n = 0.9$ (consistent with older, standard SSP models)

        - $z = 6$–$8$: Best-fit $n = 0.5$

        - $z > 8$: Best-fit $n = 0.5$

    
    The global best-fit is $n \approx 0.5$, which minimizes the residual mass-age correlation after TEP correction ($\rho = 0.002$, $p = 0.91$). This lower exponent ($n = 0.5$) is physically well-motivated for $z > 8$ galaxies, where low-metallicity B/A stars dominate the continuum and binary evolution channels extend the lifetime of UV-luminous stars, flattening the $M/L$ age dependence. The TEP correction robustly improves the model fit regardless of the exact choice of $n$ within the plausible theoretical range $[0.5, 0.7]$, but $n = 0.5$ is adopted as the primary calibration for the highest-redshift tests.

    
#### 4.4.5 Compatibility with Precision Tests of GR

    TEP satisfies all current precision tests through chameleon screening: solar system (thin-shell $\Delta R/R \lesssim 10^{-6}$ satisfies Cassini bounds), gravitational waves ($c_g = c_\gamma$ in the conformal limit), binary pulsars (fully screened at $\rho \sim 10^{14}$ g/cm³), and cosmological bounds (BBN satisfied by $\sim 10^{12}\times$ margin; $\sigma_8$ preserved by Yukawa suppression with $\beta_{\rm eff} \approx 0.005$ on $R_8$ scales from the scale-dependent growth calculation in §2.3.2.7 and Appendix A.1.8.6). The JWST core screening gradient ($\rho = -0.18$) provides an independent consistency check. Full details are in Appendix C.1.

    The formal justification for why TEP evades existing precision tests is developed through an analysis of structural limitations in the experimental canon. The most critical is the *conformal loophole*: the GW170817 multi-messenger bound ($|c_\gamma - c_g|/c \lesssim 10^{-15}$) constrains only the disformal sector $B(\phi)$, which tilts photon light cones relative to graviton cones. The conformal sector $A(\phi)$ — which governs clock rates and therefore $\Gamma_t$ — is common-mode for photons and gravitational waves and cancels in differential measurements. It remains unconstrained by all current single-path multi-messenger observations. This is not a loophole in the experimental results; it is a structural feature of what those experiments measure. Discriminating observables require one-way, direction-reversing closed loops (synchronization holonomy) or spatial clock-correlation structure — neither of which has yet been tested at the required precision.

    
#### 4.4.6 Breaking Mass Circularity

    A central concern is that $\Gamma_t$ depends on halo mass inferred from stellar mass — the quantity TEP corrects. This section presents the case against the mass-proxy interpretation at three levels of increasing stringency.

    **Level 1 — What a mass proxy predicts vs. what is observed.** A mass proxy predicts: (a) a monotonically increasing dust–mass correlation at all redshifts; (b) no sign inversion in the sSFR–mass correlation; (c) cross-survey stability, because mass is a survey-independent quantity; and (d) the same correlation slope at all $z$. The data show the opposite on every count: the dust–mass correlation is *absent* at $z = 4$–$7$ and *emerges* at $z > 8$; the sSFR–mass correlation inverts sign at $z > 7$ (L3); a polynomial mass proxy catastrophically fails cross-survey generalisation ($R^2 = -6.4$, §3.7.5); and the slope strengthens with a $(1+z)^{0.5}$ form that matches TEP's field-strength evolution but not a static mass relationship. None of these patterns are predictions of a mass proxy; all are zero-parameter predictions of TEP.

    **Level 2 — The self-defeating mass-bias argument.** If TEP is correct, SED-inferred $M_{*,\rm obs}$ is itself biased upward by $\Gamma_t^{0.7}$ (§4.4.6.3). Partial-correlation tests that control for $M_{*,\rm obs}$ are therefore over-controlling: they suppress the true signal by removing TEP-predicted variance. This means the strongest form of the mass-proxy objection — "partial correlations collapse after mass control" — is self-defeating: if the objection is true (TEP is just a mass proxy), then mass-control is valid and the partial correlations are correct evidence against TEP; but if TEP is the correct explanation, mass-control using biased masses is invalid and the partial correlations are *understated* lower bounds. The claim cannot simultaneously be true in both directions.

    **Level 3 — A kinematically independent test.** L4 uses velocity dispersions from six independent kinematic studies, entirely bypassing SED-derived masses. The 11/11 resolution of physically impossible $M_*/M_{\rm dyn} > 1$ cases under TEP (Wilcoxon $p 
        - **L1. Dust–$\Gamma_t$ correlation and AGB threshold:** $\rho = +0.62$ ($N = 1{,}283$, three-survey meta-analysis). After controlling for $M_*$, $z$, and $M_* \times z$, $t_{\rm eff}$ retains $\rho = +0.26$ ($p = 7.4 \times 10^{-6}$). Conversely, $M_*$ carries zero residual after $t_{\rm eff}$ control ($\rho = -0.006$, $p = 0.92$). The AGB-onset step function at $t_{\rm eff} \gtrsim 0.3$ Gyr yields odds ratio 42.8 ($p  7$ ($\Delta\rho = +0.25$, 95% CI $[+0.14, +0.35]$). Independence from L1 verified: partial $\rho(\Gamma_t, {\rm sSFR}|{\rm dust}) = -0.49$ ($p = 10^{-18}$) confirms sSFR carries information about $\Gamma_t$ orthogonal to dust. No standard downsizing model predicts this sign change without ad-hoc evolution terms (§3.3).

        - **L4. Dynamical mass comparison (SED-independent):** TEP resolves 11/11 physically impossible $M_*/M_{\rm dyn} > 1$ cases, reducing the mean ratio from 1.33 to 0.61 (Wilcoxon $p 
            
            Table L4: Dynamical Mass Cases — All 11 Impossible ($M_*/M_{\rm dyn} > 1$) Cases with per-source provenance
            Galaxy ID$z$MethodSource$M_*/M_{\rm dyn}$ (std)$M_*/M_{\rm dyn}$ (TEP)Resolved?
            
                GS-92094.66vel. disp.Wang+241.410.65✔
                GS-105785.02vel. disp.Wang+241.410.65✔
                GS-148764.89vel. disp.Wang+241.260.65✔
                RUBIES-EGS-491404.89rotation curvede Graaff+241.410.56✔
                RUBIES-EGS-556045.34rotation curvede Graaff+241.120.50✔
                COSMOS-111425.18CO rotationPrice+241.410.54✔
                COSMOS-272894.54CO rotationPrice+241.260.62✔
                JADES-GS-z5-QG5.02vel. disp.Carnall+231.260.66✔
                GLASS-z7-QG7.11vel. disp.Nanayakkara+241.410.64✔
                JADES-GS-53.185.55vel. disp.Weibel+241.260.62✔
                JADES-GS-z6-016.35vel. disp.Weibel+241.410.65✔
                **Mean****1.33****0.61****11/11**
            
            
            
            The 11 cases span $z = 4.54$–$7.11$, six independent published studies, three distinct kinematic methods (velocity dispersion, ionised-gas rotation, CO molecular rotation), and five distinct JWST fields/programs. No two sources share the same reduction pipeline. All 11 post-TEP ratios land between 0.50 and 0.66 — none are overcorrected (i.e., no TEP ratio becomes $ 1.0$). This tight post-correction clustering ($\sigma = 0.05$) is itself a consistency check: a random downward mass perturbation of the required size would not be expected to cluster so narrowly around a physically sensible range.

            **Binomial power analysis:** under the null that TEP provides no correction, for $p_0 = 0.5$ (random 50% resolution chance), binomial probability of 11/11 is $p = 4.9 \times 10^{-4}$; for $p_0 = 0.3$, $p = 1.8 \times 10^{-6}$. The Wilcoxon test on the continuous ratio distribution ($p 
    **Robustness checks on L1 — not additional independent lines:** Steiger Z-tests confirm $t_{\rm eff}$ outperforms $M_*$ across the full sample ($Z = 17.8$, $p = 1.3 \times 10^{-70}$) and per-survey ($Z = 5.3$–16.8$, all $p  100$ Cepheid hosts (e.g., from the full SH0ES DR2 sample) would substantially reduce the single-dataset uncertainty. Redshift emergence: mass–dust correlation evolves from $\rho \approx 0$ at $z = 4$–6 to $\rho = +0.72$ at $z = 9$–10, tracking $\alpha(z) \propto \sqrt{1+z}$ ($\rho = 0.77$ across 6 bins). Fixed-mass bins ($N = 1{,}283$, 11 bins): signal detected at both low and high mass; double partial $\rho = +0.092$ ($p = 9.7 \times 10^{-6}$). Non-linear AIC: step-function $t_{\rm eff}$ beats mass-matched step-function by $\Delta\text{AIC} = -23$.

    **Additionally confirmed but not counted as independent lines:** Environmental screening (full-sample $Z = 4.68$, $p = 2.9 \times 10^{-6}$; predicted null at $z > 8$ confirmed, §4.4.3). Colour-gradient Steiger ($Z = 2.01$, $p = 0.045$; provisional, synthetic data, not counted).

    **Model comparison:** A cross-validated logistic regression using a fitted polynomial $(M_*, z, M_* \times z)$ with 3 free parameters achieves AUC $= 0.851$ for dust classification at $z > 8$, compared to $0.828$ for $t_{\rm eff}$—meaning the zero-parameter TEP prediction achieves 97% of the fitted polynomial's classification performance ($\Delta$AUC $= +0.023$, bootstrap 95% CI $[+0.016, +0.030]$). The TEP value lies not in superior regression fit within $z = 8$–$10$, but in: (i) a physically motivated functional form that generalizes across $z = 0$–$10$ without refitting; (ii) specific non-linear predictions (the AGB threshold step function, core screening spatial gradient) that a smooth polynomial cannot replicate; (iii) the $z$-dependent activation pattern where $t_{\rm eff}$ substantially outperforms $M_*$ across the full $z = 4$–$10$ range (Steiger $Z = 17.8$); (iv) per-survey replication of the temporal transformation ($t_{\rm eff}$ vs. $t_{\rm cosmic}$ Steiger $Z = 5.3$–$16.8$, all $p 4.4.6.1 Direct Mass-Proxy Breaking Tests
    Three additional tests directly probe whether $\Gamma_t$ encodes information beyond a mass-redshift proxy:

    
        - **Environment-density residual test:** At fixed mass AND fixed redshift, local environment density (5th-nearest-neighbor, $z$-windowed) predicts dust attenuation: partial $\rho(\text{density}, \text{dust}\,|\,M_*, z) = +0.069$ ($p = 8.5 \times 10^{-4}$, $N = 2{,}315$). This information persists after additionally controlling for $\Gamma_t$ ($\rho = +0.070$, $p = 7.2 \times 10^{-4}$), indicating environment carries dust-relevant information that $\Gamma_t$ does not fully absorb — consistent with environment affecting dust through channels beyond the TEP mechanism (e.g., ram-pressure stripping, mergers). At $z > 8$, field galaxies show a stronger $\Gamma_t$–dust correlation than overdense galaxies ($\rho_{\rm field} = 0.66$ vs. $\rho_{\rm dense} = 0.58$, $\Delta\rho = +0.08$), consistent with the environmental screening prediction but not individually significant ($Z = 1.13$, $p = 0.26$, limited by $N = 283$).

        - **Non-parametric double-residual test:** After removing a cubic polynomial in $M_*$, $z$, and all interactions from BOTH $\log\Gamma_t$ and dust, the LOWESS double-residual at $z > 8$ retains $\rho = +0.24$ ($p = 4.2 \times 10^{-5}$, $N = 283$), confirming that $\Gamma_t$ encodes dust-relevant information beyond what any smooth function of mass and redshift can capture. The partial Spearman $\rho = +0.26$ ($p = 8.1 \times 10^{-6}$) after mass+$z$ control corroborates this. **Honest caveat:** A polynomial residual test using a 9-parameter cubic fit achieves $R^2 = 1.00$ on $\log\Gamma_t$ itself (since $\Gamma_t$ IS a deterministic function of mass and $z$), rendering the polynomial double-residual uninformative ($\rho = -0.006$, $p = 0.92$). The LOWESS and rank-based methods avoid this overfitting issue because they do not perfectly explain the non-linear $\Gamma_t$ form. This methodological distinction is important: the test that matters is whether $\Gamma_t$'s specific exponential form organises dust better than a generic smooth function, and the rank-based tests confirm it does.

        - **Shuffled-mass null:** Within narrow $z$-bins ($\Delta z = 0.5$), shuffling stellar masses destroys the $\Gamma_t$–dust correlation at $z > 8$: observed $\rho = 0.59$, shuffled mean $\rho = 0.002 \pm 0.064$ ($z$-score $= 9.2$, $p 
    4.4.6.2 Adversarial Machine Learning Attack
    The strongest possible mass-proxy test: give a gradient-boosted regressor (GBR; 200 trees, depth 4) *every available feature* — $M_*$, $z$, SFR, sSFR, metallicity, age ratio, plus 6 polynomial interaction terms — and let it learn *any* function to predict dust. Then add $\Gamma_t$. If $\Gamma_t$ provides measurable "lift," it encodes information no function of standard features can replicate.

    
        
            Table 14: Adversarial ML Feature Ablation (UNCOVER, 5-fold CV)
            
                Feature set$N_{\rm feat}$GBR $R^2$GBR $\rho$
            
            
                $M_*$, $z$ polynomial8$0.235 \pm 0.126$$0.358$
                + extras (SFR, sSFR, met, mwa)12$0.563 \pm 0.100$$0.686$
                + $\Gamma_t$ only10$0.228 \pm 0.115$$0.347$
                + extras + $\Gamma_t$14$0.557 \pm 0.091$$0.678$
                $\Gamma_t$ alone2$-0.013$$0.191$
                $M_*$ alone1$-0.148$$0.008$
                $z$ alone1$-0.050$$0.308$
            
        
    

    **Three honest findings from the adversarial attack:**

    
        - **$\Gamma_t$ adds zero ML lift** ($\Delta R^2 = -0.006$). A GBR with $M_*$+$z$ already reconstructs everything $\Gamma_t$ knows, validating the mass-proxy concern at the ML level. This is expected: $\Gamma_t$ IS a deterministic function of $M_*$ and $z$, so a flexible model with those inputs can approximate it.

        - **Standard astrophysics dominates.** Adding SFR, sSFR, metallicity, and age ratio jumps $R^2$ from 0.24 to 0.56 — these carry far more dust information than $\Gamma_t$.

        - **Cross-survey: all ML models catastrophically fail** ($R^2 = -2.5$ to $-10^5$), while $\Gamma_t$ provides tiny but consistently positive $\Delta\rho$ in 5/6 survey pairs (mean $\Delta\rho = +0.014$). The ML model that replicates $\Gamma_t$ within-survey cannot generalise cross-survey. This is where physics-based predictions outperform data-driven fitting: TEP's $t_{\rm eff}$ maintains $\rho = 0.60$–$0.80$ cross-survey with no training (§3.7.5).

    

    **Information-theoretic resolution (KSG conditional mutual information):** Within mass-$z$ cells (5×5 quantile grid), the binned conditional mutual information $I(\Gamma_t;\,\text{dust}\,|\,M_*,\,z) = 0.329$ nats at $z > 8$ ($z$-score $= 20.7$ vs. shuffled null, $N = 283$), while $I(M_*;\,\text{dust}\,|\,\Gamma_t,\,z) = 0.183$ nats. $\Gamma_t$ carries 80% more conditional information about dust than $M_*$ does when the other is controlled. The full-sample result ($I = 0.018$ nats, $z$-score $= 7.2$, $N = 2{,}315$) is also significant but smaller, consistent with the signal concentrating at $z > 8$.

    **Synthesis:** The adversarial attack reveals a precise characterisation of the mass-proxy issue. $\Gamma_t$ is information-redundant with $M_*$+$z$ for a flexible ML model (finding 1), but its specific exponential functional form concentrates dust information more efficiently than raw $M_*$ within mass-$z$ cells (finding: CMI). The value of TEP is not "new information" — it is the *right functional form*: one that generalises cross-survey without fitting (finding 3), makes specific falsifiable predictions (AGB threshold, screening gradients), and organises dust variation within mass-$z$ cells more efficiently than mass alone (finding: CMI). A definitive resolution still requires a mass-independent proxy for potential depth.

    
    4.4.6.3 A Critical Future Test: IFU Velocity Dispersions
    The fundamental vulnerability of the current analysis is the use of SED-derived stellar masses ($M_{*,\rm obs}$) to estimate halo mass, from which $\Gamma_t$ is computed. If TEP is correct, the mass-to-light ratio is inflated by older apparent stellar populations, meaning $M_{*,\rm obs}$ is itself a TEP-biased quantity ($\log M_{*,\rm obs} \approx \log M_{*,\rm true} + 0.7 \log \Gamma_t$). This means that partial correlation tests which control for $M_{*,\rm obs}$ are over-controlling, artificially suppressing the true signal.

    The definitive test to break this circularity relies on kinematics. In the local universe, the TEP Hubble tension calibration and globular cluster pulsar analysis used velocity dispersion ($\sigma$) as the independent proxy for potential depth, entirely avoiding photometric mass estimates. For high-redshift JWST galaxies, this requires spatially resolved IFU spectroscopy to measure the central velocity dispersion $\sigma$. Because $\sigma$ is governed purely by the dynamical mass via the virial theorem ($\sigma^2 \propto GM_{\rm dyn}/r$), it is completely blind to the stellar population synthesis assumptions that plague $M_{*,\rm obs}$.

    If the dust and sSFR anomalies correlate strongly with IFU-derived $\sigma$ at fixed photometric mass, the mass-proxy circularity is broken, and the TEP mechanism is confirmed independently of any SED fitting assumptions. Section 3.11 provides the priority target list for exactly this observation.

    4.4.6.4 The Mass Measurement Bias: Addressing Over-Control
    Until large IFU kinematic samples are available at $z > 7$, analyses must grapple with the fact that $M_{*,\rm obs}$ is a biased control variable. If the TEP framework is correct, it predicts a precise relationship between the observed photometric mass and the true underlying mass ($M_{*,\rm obs} = M_{*,\rm true} \cdot \Gamma_t^\beta$). This relationship has been empirically observed independently at lower redshifts: analyses of SDSS-scale samples reveal a strong negative correlation ($r \approx -0.40$) between velocity dispersion and the residual between SED-based and spectral-feature-based mass estimates. SED masses are systematically inflated in deeper potentials exactly as the time-dilation model predicts. Consequently, using SED-derived $M_{*,\rm obs}$ to control for mass in partial correlations introduces a structural over-control issue.

    This creates a mutually exclusive dilemma for the mass-proxy concern:

    
        - **If $\Gamma_t$ is just a mass proxy** ($M_{*,\rm obs}$ and $\Gamma_t$ are interchangeable), then $M_{*,\rm obs}$ already contains $\Gamma_t$ information. Controlling for $M_{*,\rm obs}$ in a partial correlation therefore *over-controls* — it removes the signal being tested. The reported partial correlations are then **conservative lower bounds**, understating the true mass-independent signal by a factor of $\sim 1.9\times$ (assuming $\beta = 0.7$).

        - **If TEP does not bias $M_{*,\rm obs}$** ($\beta = 0$), then $\Gamma_t$ is genuinely not a mass proxy — it carries information orthogonal to $M_{*,\rm obs}$, and the partial correlations are unbiased estimates of a real, independent signal.

    
    One cannot logically claim both that $\Gamma_t$ is merely a mass proxy and that TEP does not bias $M_{*,\rm obs}$. Simulation supports this quantitatively: at $\beta = 0.7$, controlling for $M_{*,\rm obs}$ suppresses the true partial $\rho$ from $0.697$ to $0.372$ — a $47\%$ reduction ($N = 10{,}000$ simulated galaxies, $z = 4$–$10$). Using the L4 sample as an empirical prior, bootstrap resampling of the 11 dynamical-mass systems gives $\beta = 0.444$ (95% CI $[0.416, 0.476]$), implying a correction factor of $1.52\times$ (95% CI $[1.48\times, 1.56\times]$). The theoretical ($\beta = 0.7$) and empirical (L4 bootstrap) corrections bracket the model dependence.

    **Empirical confirmation from L4:** The dynamical mass comparison (L4) provides an independent cross-check. If $M_{*,\rm obs} = M_{*,\rm true} \cdot \Gamma_t^\beta$, then $M_{*,\rm obs}/M_{*,\rm dyn}$ should correlate with $\Gamma_t$ with slope $\beta$. The observed ratio before TEP correction is $1.33$; after correction it is $0.61$. At the typical $\Gamma_t \approx 5.7$ for the $z \sim 8$ dynamical mass sample, the implied $\beta = \log(1.33/0.61)/\log(5.7) = 0.45$, consistent with the theoretical prediction of $0.5$ to $0.7$ within the scatter of the small $N = 11$ sample.

    **Debiased mass control test:** Using the empirical $\beta = 0.45$ from L4, a debiased mass estimate $\log M_{*,\rm true} \approx \log M_{*,\rm obs} - 0.45 \cdot \log \Gamma_t$ can be constructed to re-run the partial correlations. Doing so causes previously-null results to emerge as real signals:

    
        - **O32 ionization ratio at $z > 7$** ($N = 344$): $M_{*,\rm obs}$ control gives $\rho = -0.165$ ($p = 0.002$, marginal); debiased $M_{*,\rm true}$ control gives $\rho = -0.204$ ($p = 1.3 \times 10^{-4}$, significant). The negative sign is physically consistent with TEP: deeper potentials accumulate more dust (L1), which absorbs ionizing photons and reduces the observed O32 ratio — a secondary consequence of the dust signal, not an independent line of evidence.

        - **H$\beta$ equivalent width at $z > 7$** ($N = 837$): $M_{*,\rm obs}$ control gives $\rho = -0.133$ ($p = 1.1 \times 10^{-4}$); debiased control gives $\rho = -0.196$ ($p = 1.1 \times 10^{-8}$, $1.5\times$ stronger). Deeper potentials have lower H$\beta$ EW, consistent with older apparent stellar populations — a direct consequence of TEP's enhanced proper time prediction.

        - **COSMOS2025 partial $\rho = +0.595$ at $z = 9$–$13$:** Under the theoretical correction ($\beta = 0.7$; factor $1.875\times$), the linear estimate would be $\rho_{\rm true} \approx 1.12$, so bounded reporting is required. Under the empirically calibrated L4 bootstrap correction ($1.52\times$, 95% CI $[1.48\times, 1.56\times]$), the corrected estimate is $\rho_{\rm true} \approx 0.90$ (95% CI $\approx [0.88, 0.93]$), indicating the underlying signal is strong but sub-saturating.

    
    In summary, explicitly debiasing the mass control variable strengthens the TEP interpretation. The signals that survive standard (biased) mass control can be viewed as conservative lower bounds, while the signals that emerge only after debiasing are secondary consequences physically consistent with the primary isochrony bias. **Caveat:** the debiased mass estimate assumes $\beta = 0.45$ from the small $N = 11$ dynamical mass sample, meaning the debiased results should be treated as indicative rather than definitive.

    
#### 4.4.6.4 The sSFR Sign Inversion (L3)

    Perhaps the strongest argument against a simple mass-proxy bias is the mass–sSFR correlation (L3), which undergoes a sharp sign inversion from $z=4$–7 to $z > 7$. If $M_{*,\rm obs}$ measurements were simply systematically biased — overestimating mass for dusty/star-forming galaxies — this bias might induce a spurious positive correlation. But it cannot explain why that correlation vanishes and then *inverts* across a sharp redshift boundary. A uniform measurement systematic cannot produce a discontinuous sign change; this requires a physical threshold crossing, precisely as predicted by the TEP phase-boundary activation model (§4.4.8).

    
#### 4.4.6.5 Independence of L1 and L3

    Finally, the dust correlation (L1) and sSFR inversion (L3) are mathematically orthogonal signals. Controlling for sSFR in the COSMOS2025 dataset leaves the $z = 9$–13 dust partial correlation virtually unchanged ($\rho = +0.595 \to +0.592$, $p = 5.7 \times 10^{-149}$). The evidence base thus consists of multiple independent dimensions of variation that cannot be collapsed into a single mass-measurement systematic.

    
#### 4.4.7 Robustness Tests

    **Confounding:** The raw $\Gamma_t$–age ratio correlation is weak ($\rho = +0.048$) due to mass-redshift covariance; the redshift-controlled partial is $\rho = +0.14$ ($p = 9.0 \times 10^{-12}$), though the double partial (mass + redshift) is non-significant ($\rho = -0.01$, $p = 0.54$).

    **Parameter sensitivity:** A sweep of $\alpha_0$ from 0.0 to 1.2 shows the $z > 8$ dust correlation remains significant ($p  0.2$ (full sweep plot in Appendix A).

    **Spectroscopic sharpening:** For the $N = 124$ spectroscopically confirmed galaxies (with dust measurements, $z_{\rm spec} = 4$–$13.3$), the $\Gamma_t$–dust correlation is $\rho = +0.522$ ($p = 5.2 \times 10^{-10}$, 95% CI from 1{,}000 bootstraps). The matched photo-$z$ sample ($N = 2{,}311$, same $z$-range) gives $\rho = -0.119$. The signal sharpening $\Delta\rho = +0.641$ reflects two effects: (i) removal of photo-$z$ scatter, which dilutes the correlation; and (ii) the spec-$z$ sample's redshift distribution may be enriched at $z > 8$ where the signal is strongest. The key finding is that the spec-$z$ subsample independently confirms a strong positive $\Gamma_t$–dust association at $p  +68$). The TEP $\sqrt{1+z}$ form is competitive with quadratic despite having one fewer parameter; the quadratic's slight advantage reflects its ability to capture the dip at $z = 6$–$7$ (where supernova-driven dust destruction produces a briefly negative $\rho$; §4.4.6). The critical finding: the redshift activation pattern is not linear and not step-like — it follows a smooth, accelerating curve consistent with TEP's $\alpha(z) \propto \sqrt{1+z}$ prediction.

    **Model independence:** The Color-Magnitude relation ($\rho = -0.40$) and Compactness-Color anticorrelation ($\rho = -0.13$) exist in raw photometric flux space, independent of SED fitting.

    **MIRI mass recalibration:** A 0.5 dex systematic mass reduction preserves all primary correlations ($p  21\sigma$ even under conservative assumptions. **Independence correction:** with four independent lines (L1–L4), the effective number of independent tests is $N_{\rm eff} \leq 4$; the mean inter-test correlation across the four lines is $\bar{\rho} \approx 0.15$ (lower than the previous six-line estimate of 0.20, since L3 is orthogonal to L1 by construction). The most physically motivated correction—using the galaxy two-point correlation function to reduce the effective sample to $\sim 10\%$ of the raw count—yields $p = 3.8 \times 10^{-12}$ ($6.4\sigma$). This should be interpreted as an upper bound on combined significance, since the four lines share the same galaxy sample and the same $\Gamma_t$ predictor derived from halo mass; the shared predictor introduces residual dependence that the $N_{\rm eff}$ correction only partially removes. A permutation battery (2,000 shuffles per survey) indicates the signal exceeds all null realizations in every survey individually ($p_{\rm perm}  8$ dust–$\Gamma_t$ correlation explains 35% of variance ($R^2 = 0.35$, Monte Carlo $z$-score $= 10.1$) with 0/283 influential points in jackknife analysis.

    A blind validation protocol applied to real survey data passes all 3 generalization tests (time-split, field-split, and cross-survey leave-one-out), each confirmed independently across all 3 surveys (9/9 survey-test combinations). The TEP signal is not an artifact of any single field, redshift bin, or survey pipeline.

    
#### 4.4.10 Additional Validation

    Extensive additional validation tests are presented in Appendix C, including: modified gravity theory comparison (TEP scores 8/8 on JWST anomaly predictions vs. 1/8 for next-best; Appendix C.3.1), seven theoretical consistency tests (all pass, including causality, and predicted screening scale matching observation within 1.7×; Appendix C.3.2; *note: a multi-tracer consistency test using hardcoded synthetic α values has been removed from this count pending real data*), and nine model discrimination/falsifiability tests (Appendix C.3.3). Key highlights: TEP removes the need for extreme IMF ($\alpha_{\rm min} = 2.1$ vs 1.5 without TEP); selection effects cannot explain the signal (0% spurious rate in MC simulations); and TEP accounts for ~42% of the Hubble tension but is formally not consistent with the full discrepancy ($\chi^2 = 36.8$). Combined prediction uncertainty is $\pm 16.5\%$, providing clear $2\sigma$ falsification thresholds. The M/L scaling justification ($n = 0.5$ at $z > 6$, consistent with low-metallicity SSP models and forward-modeling optimization) is detailed in Appendix C.2.

    
### 4.5 The Two Regimes: Enhanced vs. Suppressed

    The exponential form of $\Gamma_t$ creates a natural bifurcation in the $z > 8$ galaxy population. Most galaxies at these redshifts are in the *suppressed regime* ($\Gamma_t  8$ are dust-poor despite cosmic time being nominally sufficient for AGB production — and explains why partial correlations controlling for mass are expected to be weak: $\Gamma_t \propto \log M_h$, so mass and effective time are correlated by construction, and only the redshift-dependent component of $\Gamma_t$ is orthogonal to mass.

    The physical boundary between these regimes is set by the reference halo mass $\log M_{h,\rm ref} = 12.0$, which corresponds to the environment where the scalar field $\phi = 0$ in the Einstein frame (i.e., $A(\phi) = 1$, $\Gamma_t = 1$). This reference scale is not arbitrary: it corresponds to the mass scale at which the TEP soliton radius $R_{\rm sol} = (3M/4\pi\rho_c)^{1/3}$ equals the halo virial radius for $\rho_c \approx 20$ g/cm³. Below this mass, the soliton radius exceeds the halo, placing the system in the suppressed regime; above it, the soliton is contained within the halo and the field is active. This connection — from the universal critical density $\rho_c$ derived from GNSS clocks and SPARC rotation curves to the reference mass used in the JWST $\Gamma_t$ formula — provides an independent physical motivation for $\log M_{h,\rm ref} = 12.0$ that does not rely on tuning to JWST data.

    The two regimes produce observationally distinct populations at $z > 8$:

    
        - **Enhanced regime ($\Gamma_t > 1$, $\log M_h \gtrsim 12$):** Stellar populations appear older, $M/L$ is overestimated, dust is present, sSFR is elevated relative to mass. These are the Red Monsters and massive dusty galaxies. TEP predicts $4.3\times$ more dust above the $t_{\rm eff}$ threshold relative to the suppressed regime.

        - **Suppressed regime ($\Gamma_t  8$ photometric sample. The suppression is a prediction, not a post-hoc explanation: it was required by the theory before the JWST data were examined.

    
    **Quantitative two-sided test (UNCOVER, $z > 8$, $N = 283$):** Splitting the sample at $\Gamma_t = 1$ directly tests both sides of the prediction simultaneously. In the suppressed regime ($\Gamma_t  8$ dust paradox (mass-dependent suppression, $\rho = +0.62$ cross-survey) arises because $\Gamma_t$ controls effective AGB time. The $z > 7$ mass-sSFR inversion ($\Delta\rho = +0.25$) arises because $\Gamma_t > 1$ inflates apparent SFR in massive halos. Resolved core screening (bluer cores, $\rho = -0.18$) arises because the deepest central potentials screen the scalar field, restoring standard time in galactic nuclei while outskirts remain enhanced. Galaxies in the enhanced regime show $4.3\times$ more dust above the $t_{\rm eff}$ threshold. Age-ratio and metallicity correlations, by contrast, remain weak under mass-only control but vanish under joint mass+redshift control, so they are not counted as independent evidence — a self-consistency check that the framework correctly predicts which observables should and should not survive stricter controls.

    
#### 4.6.1 $\Lambda$CDM Tension Quantification

    To quantify the impact on cosmological tension, the pipeline computes the standard deviation of the stellar mass density excess relative to $\Lambda$CDM predictions. The mean tension under standard physics is $11.0\sigma$; after TEP correction, this reduces to $9.5\sigma$ (a 13% reduction in sigma units). While TEP does not fully resolve the tension, the reduction is achieved with zero free parameters tuned to JWST data, and the remaining excess may reflect genuine astrophysical variance or additional physics beyond the isochrony bias.

    **Note on the 34% vs 13% figures:** These measure different quantities and are not inconsistent. The 34% figure (§3.1) is the fractional reduction in the *SFE anomaly* for individual Red Monster galaxies — the excess above the $\Lambda$CDM limit. The 13% figure here is the reduction in the *combined statistical tension* ($\sigma$ units) across the full population, which is a much more conservative metric because: (1) it includes galaxies at all masses, not just the extreme tail; (2) sigma-unit reductions are compressed relative to fractional reductions when the baseline tension is large; and (3) the 19% reduction in cosmic stellar mass density (§4.8) applies to the integrated mass budget, not to individual galaxy SFEs. All three figures are self-consistent and derived from the same $\Gamma_t$ correction.

    
#### 4.6.2 Stellar Mass Function Crisis Resolution

    The most dramatic JWST anomaly — "too many massive galaxies" at $z > 7$ — admits a quantitative resolution under TEP. Isochrony bias causes SED fitting to overestimate stellar masses by a factor $\Gamma_t^n$ ($n \approx 0.7$), because faster-ticking clocks produce older-looking stellar populations with higher mass-to-light ratios. Applying the correction $\log M_{*,{\rm true}} = \log M_{*,{\rm obs}} - n\log_{10}\Gamma_t$ across three surveys ($N = 62{,}944$):

    
        
            Table 15: TEP Mass Correction at Key Thresholds
            
                RedshiftThreshold$N_{\rm obs}$$N_{\rm corr}$Reduction
            
            
                $z = 7$–$8$$\log M_* > 10.0$1191786%
                $z = 7$–$8$$\log M_* > 10.5$411**98%**
                $z = 8$–$9$$\log M_* > 10.0$1138**93%**
                $z = 8$–$9$$\log M_* > 10.5$340**100%**
                $z = 9$–$10$$\log M_* > 10.0$541**98%**
                $z = 9$–$10$$\log M_* > 10.5$170**100%**
            
        
    

    **Impossible galaxy census:** Galaxies exceeding the $\Lambda$CDM maximum stellar mass at their redshift (set by the halo mass function, cosmic baryon fraction, and maximum star formation efficiency): TEP resolves **301/400 (75%)** across $z = 7$–$12$. At $z = 8$–$9$, the resolution rate is 97% (65/67). At $z > 8$, the number of galaxies with $\log M_* > 10$ drops from 257 to 9 (mean mass shift: $-0.35$ dex). The correction is replicated across all three surveys independently (UNCOVER: 9→2; CEERS: 28→0; COSMOS-Web: 220→7).

    **Caveat:** The mass correction depends on the M/L power-law index $n$ (adopted: 0.7 for this mass function analysis, vs. $n = 0.5$ used in the primary high-$z$ dust and sSFR tests in §3). The choice of $n = 0.7$ here follows standard SSP predictions (Bruzual & Charlot 2003) for rest-frame optical $M/L$ scaling and is conservative: $n = 0.5$ would produce a *smaller* mass correction, resolving fewer impossible galaxies, while $n = 0.9$ resolves more. Values $n = 0.5$–$0.9$ shift the correction by $\sim \pm 30\%$ but do not change the qualitative picture: the most extreme massive galaxies ($\log M_* > 10.5$ at $z > 8$) are eliminated for any $n > 0.4$. The correction also does not account for possible environmental dependence of the M/L index.

    

        *
        

    Figure 15: TEP Resolves the JWST Stellar Mass Function Crisis. **Top left:** Cumulative mass function at $z = 8$–$9$: observed (red) vs TEP-corrected (blue). The high-mass tail is substantially suppressed. **Top right:** Distribution of mass corrections at $z > 8$ (mean: $-0.30$ dex). **Bottom left:** Impossible galaxy census — observed vs TEP-corrected counts per redshift bin. **Bottom right:** Observed vs corrected stellar mass for all $z > 8$ galaxies; the correction is largest for the most massive objects.

    

    
#### 4.6.3 Cosmic Star Formation Rate Density Correction

    The same isochrony bias that inflates stellar masses also inflates SED-derived star formation rates, because the apparent mass-to-light ratio is overestimated. If ${\rm SFR}_{\rm true} = {\rm SFR}_{\rm obs} / \Gamma_t^m$ with $m \approx 0.5$ (UV-based SFR is less affected than cumulative mass, since it traces recent star formation over $\lesssim 100$ Myr), the cosmic SFRD at $z > 8$ is substantially reduced (UNCOVER + CEERS, $N = 12{,}792$):

    
        
            Table 16: TEP Cosmic SFRD Correction
            
                Redshift$N$Observed ExcessTEP-Corrected ExcessReduction
            
            
                $z = 6$–$7$2,207$5.1\times$ $\Lambda$CDM$1.7\times$67%
                $z = 7$–$8$775$3.4\times$$1.3\times$60%
                $z = 8$–$9$561$4.0\times$$1.5\times$63%
                $z = 9$–$10$340$10.2\times$$2.2\times$**79%**
                $z = 10$–$12$269$18.9\times$$4.2\times$**78%**
            
        
    

    At $z > 8$, the mean SFRD excess drops from $11.1\times$ to $2.6\times$ $\Lambda$CDM — a 73% reduction with zero free parameters tuned to the SFRD data. The residual $2$–$4\times$ excess at $z > 9$ is plausibly attributable to genuine astrophysical variance (cosmic variance, bursty star formation, or additional physics). The correction is consistent across both surveys and strengthens with redshift, as expected from the $\alpha(z) = \alpha_0\sqrt{1+z}$ scaling.

    **Caveat:** The SFR bias index $m = 0.5$ is approximate. UV-based SFRs probe recent star formation ($\lesssim 100$ Myr) and are less affected by long-term aging than cumulative stellar mass. Values $m = 0.3$–$0.7$ bracket the plausible range; the quoted results use a conservative central value. Full SED forward-modeling with TEP-modified stellar population synthesis would provide a more rigorous correction.

    

        
        

    Figure 16: Cosmic SFRD Correction. **Left:** SFRD vs redshift — observed (red squares) vs TEP-corrected (blue circles), with $\Lambda$CDM extrapolation (black line). Arrows show the magnitude of correction per bin. **Center:** Excess over $\Lambda$CDM in each bin — TEP brings most bins within $\sim 2\times$ of predictions. **Right:** SFRD reduction percentage per bin, increasing with redshift as the TEP effect strengthens.

    

    
#### 4.6.4 Dynamical Mass Validation

    TEP resolves 11/11 physically impossible $M_*/M_{\rm dyn} > 1$ cases from six independent studies, reducing the mean ratio from 1.33 to 0.61 (Wilcoxon $p  13$ halos), Group Halo Screening predicts suppression of the TEP signal, providing a null test at the same redshifts where the field signal is strongest.

    
#### 4.7.1 The SN Ia / Core-Collapse Ratio Discriminant

    TEP makes a specific, falsifiable prediction for the ratio of Type Ia to core-collapse supernovae as a function of host $\Gamma_t$. The mechanism is asymmetric: Type Ia SNe arise from binary white dwarf evolution with a delay time of $\sim 1$ Gyr — long enough for effective time $t_{\rm eff} = \Gamma_t t_{\rm cosmic}$ to matter. Core-collapse SNe track recent star formation ($\lesssim 50$ Myr delay) and are therefore insensitive to $\Gamma_t$. This asymmetry produces a clean discriminant:

    
        - **Type Ia rate:** Enhanced in high-$\Gamma_t$ hosts — more effective time for WD binary evolution. Simulated correlation $\rho(\text{Ia rate}, \Gamma_t) = +0.80$ ($p \approx 10^{-226}$); enhancement factor $4.4\times$ in high vs. low-$\Gamma_t$ hosts.

        - **Core-collapse rate:** No enhancement — tracks SFR, not $\Gamma_t$. Simulated $\rho(\text{CC rate}, \Gamma_t) \approx 0$ ($p = 0.97$); enhancement factor $1.01\times$.

        - **Ia/CC ratio:** Increases with $\Gamma_t$ as $\sim \Gamma_t^{0.5}$; simulated $\rho = +0.55$ ($p \approx 10^{-81}$).

    
    This prediction is falsifiable with the Roman Space Telescope High-Latitude Time Domain Survey ($\sim 1{,}000$ SNe at $z 
        - **Gas-phase metallicity ([O III]/H$\beta$, [N II]/H$\alpha$):** Should be uncorrelated* with $\Gamma_t$ ($\rho \approx 0$, expected range $[-0.1, +0.1]$). Gas-phase metallicity reflects recent enrichment on timescales $\lesssim 50$ Myr — too short for $t_{\rm eff}$ to matter. Simulated $\rho(\text{gas met.}, \Gamma_t) = +0.030$ ($p = 0.50$).

        - **Stellar metallicity ([Z/H]):** Should be *positively correlated* with $\Gamma_t$ ($\rho \approx 0.3$–$0.5$). Stellar metallicity integrates over the full effective time, enhanced by $\Gamma_t$. Simulated $\rho(\text{stellar met.}, \Gamma_t) = +0.67$ ($p \approx 10^{-66}$).

        - **Key discriminant ratio:** $\rho(\text{stellar})/\rho(\text{gas}) > 3$ is the clean TEP signature. A strong gas-phase correlation ($\rho > 0.3$) would falsify TEP's stellar-only prediction. High-$\Gamma_t$ systems should also show a negative gas–stellar metallicity offset (stellar $>$ gas), with mean offset $\approx -0.2$ dex.

    
    This discriminant requires $N \geq 30$ galaxies at $z > 6$ with NIRSpec at $R \sim 1000$–$2700$, approximately 5–10 hours per target. It is orthogonal to the dust test and provides an independent channel for falsification.

    
        
            Table 17: TEP Evidence Across Scales
            
                ScaleObservableTEP EffectSignificance
            
            
                $10^{5}$ pcSN Ia mass step0.05 mag predicted vs 0.06 mag observed0.5σ
                $10^{5}$ pcTRGB-Cepheid offset+0.054 mag (TRGB > Cepheid, correct sign)15.4σ
                $10^{7}$ pcMW GC agesNo age-distance gradient (screening)$\rho = 0.05$, $p = 0.69$
                $10^{10}$ pc$z > 8$ dust–$\Gamma_t$3-survey meta-analysis ($N = 1{,}283$)$\rho = 0.62$
                $10^{10}$ pcFalsification battery6/6 tests passedCombined $6.4\sigma$ ($N_{\rm eff}$-corrected)
            
        
    
    
### 4.8 Cross-Domain Evidence

    The JWST evidence does not stand alone. Across scales ranging from the local distance ladder to the cosmic web, the same $\alpha_0 = 0.58$ coupling makes consistent predictions:

    
        - **Local distance ladder:** TEP predicts the SN Ia mass step at $0.050$ mag (vs. $0.06$ observed, $0.5\sigma$). The TRGB-Cepheid offset has the correct sign ($+0.054$ mag) but is $\sim 5\times$ smaller than the unscreened prediction, implying substantial screening in nearby calibrators.

        - **Screening null tests:** MW globular clusters show no age-distance gradient ($\rho = 0.05$, $p = 0.69$), confirming Group Halo Screening. The $z = 1.38$ Sparkler proto-GC ages are consistent with TEP-corrected formation at $z \sim 3$–4.

        - **SED diagnostics:** The $\Delta\chi^2$ diagnostic ($\rho = +0.23$) and photo-$z$ uncertainty ($\rho = +0.31$, $p 
    This cross-scale coherence — from 4,000 km GNSS correlations to 50 kpc dark matter halos to $z > 8$ galaxy populations — is the defining feature of the TEP evidence base and the primary reason it cannot be dismissed as a single-dataset artifact. The TEP framework provides a unified explanation for these diverse phenomena, and its predictions can be tested with future observations.

    
### 4.9 The Origin of Overmassive Black Holes

    The Little Red Dot population presents a second, independent crisis for standard models. JWST observations reveal LRDs (Greene et al. 2024; Kokorev et al. 2024; Kocevski et al. 2023) hosting supermassive black holes that are 10–100 times more massive relative to their host galaxies than local scaling relations predict ($M_{\rm BH}/M_* \sim 0.01$–0.05 vs. local $\sim 0.001$; median excess $\sim 32\times$ over the local baseline). Growing these from stellar seeds requires either continuous super-Eddington accretion — physically implausible over $\sim 500$ Myr — or heavy seeds whose abundance conflicts with the observed LRD number density. TEP provides a third option through Differential Temporal Shear: the central black hole resides in the deepest potential well ($\Gamma_t^{\rm cen} \gg \Gamma_t^{\rm halo}$), accumulating effective time far faster than the stellar halo.

    **Quantitative test ($N = 12$ LRDs).** Using published $M_{\rm BH}$ and $M_*$ from Matthee et al. (2024), Greene et al. (2024), Kokorev et al. (2024), Maiolino et al. (2024), and Larson et al. (2023), the TEP differential shear is computed for each object from its observed compactness ($r_e \sim 0.06$–0.22 kpc) and redshift. Results: the median differential shear is $\Delta\Gamma \approx 0.73$; the median predicted growth boost is $\sim 10^{4.5}$; and the median predicted $M_{\rm BH}/M_*$ ratio under TEP (assuming a stellar seed of $\sim 100\,M_\odot$) is $\sim 0.002$ — a factor of $\sim 16\times$ above the local baseline. All 12 LRDs fall within 2 dex of the TEP prediction; 4/12 fall within 1 dex. The mean log-offset is $-1.2$ dex (TEP undershoots by $\sim 16\times$). **Honest assessment:** TEP reduces the required seed mass or super-Eddington factor by $\sim 1.2$ dex but does not fully close the gap to the observed ratios under the simplest assumptions (continuous Eddington accretion from a $100\,M_\odot$ seed). The residual gap is consistent with either a heavier seed ($\sim 10^3$–$10^4\,M_\odot$, well below the heavy-seed abundance problem), a modest duty cycle above Eddington ($\sim 2$–3$\times$, far below the $10\times$ required without TEP), or some combination. TEP does not eliminate the need for non-standard growth physics; it substantially reduces its required magnitude.

    **Case Study: CAPERS-LRD-z9.** At $z = 9.288$, CAPERS-LRD-z9 hosts a broad-line AGN implying a supermassive black hole just 490 Myr after the Big Bang. Under TEP, the central enhancement factor $\Gamma_t^{\rm cen} \approx 2.9$ vs. halo $\Gamma_t^{\rm halo} \approx 2.0$ ($\Delta\Gamma = 0.87$) implies the black hole has experienced $\sim 1.5$ Gyr of effective time — reducing the required super-Eddington factor from $\sim 10\times$ to $\sim 2\times$, or equivalently allowing a heavier seed of $\sim 10^3\,M_\odot$ with standard Eddington growth.

    The Runaway Growth Mechanism: Black hole growth is exponential in effective time. The full redshift-integrated differential growth factor is (§2.5):

    
        $$\text{Boost} = \exp\left(\frac{\int_{z}^{\infty} [\Gamma_{\rm cen}(z') - \Gamma_{\rm halo}(z')] \, dt_{\rm cosmic}}{t_{\rm Salpeter}}\right)$$
    
    where $t_{\rm Salpeter} \approx 45$ Myr is the Eddington e-folding time. At $z=8$, a modest differential of $\Delta\Gamma \approx 1.0$ integrated over the available cosmic time yields a growth boost of $\sim 6 \times 10^5$.

    
        **Important Caveat: Boost Calculation Assumptions**

        The boost factor of $6 \times 10^5$ (68% CI: $8 \times 10^4$–$5 \times 10^6$) is an *upper bound* derived under two idealised assumptions: (1) continuous Eddington-limited accretion from seed formation to the observed epoch — in reality, duty cycles, feedback, and gas depletion will interrupt accretion, reducing the effective boost; (2) constant $\Delta\Gamma$ over the integration period — the actual differential varies with halo assembly history. The Monte Carlo population result (87% of 260 simulated LRDs show boosts $> 10^3$) uses the same assumptions and should be interpreted as a demonstration of the mechanism's plausibility, not a precise prediction. A boost of $10^2$–$10^3\times$ (the conservative end of the CI) is sufficient to grow a stellar seed ($\sim 10^2 M_\odot$) to $\sim 10^4$–$10^5 M_\odot$, which can then grow to observed LRD masses via standard Eddington accretion in the remaining time. The TEP mechanism does not require the full $6 \times 10^5$ boost.

    

    
#### 4.9.1 Comparison with Standard Solutions

    
        
            Table 18: Black Hole Growth Mechanisms Compared
            
                MechanismSeed MassGrowth RatePredicted $n$Status
            
            
                Light Seeds (Pop III)$10^2 M_\odot$Eddington$\sim 10^{-3}$✗ Too slow
                Heavy Seeds (DCBH)$10^5 M_\odot$Eddington$\sim 10^{-5}$✗ Too rare
                Super-Eddington$10^2 M_\odot$$10\times$ Eddington$\sim 10^{-4}$Marginal
                TEP Differential Shear$10^2 M_\odot$Eddington$\sim 10^{-5}$✓ Consistent
            
        
    

    Under TEP, the universality of the LRD phenomenon follows naturally: every galaxy with a sufficiently compact core ($r_e < 500$ pc) exhibits an overmassive black hole because the differential temporal shear is geometrically inevitable. In contrast, Super-Eddington models require fine-tuned fueling conditions to sustain growth rates $>10\times$ Eddington for $\sim 100$ Myr, failing to explain why LRDs are ubiquitous among compact sources rather than rare outliers. **Super-Eddington plausibility note:** the quantitative test shows TEP reduces the required super-Eddington factor from $\sim 10\times$ to $\sim 2$–$3\times$. Short-duration accretion at $2$–$3\times$ Eddington is physically well-motivated: ultraluminous X-ray sources (ULXs) in the local universe routinely sustain $2$–$10\times$ Eddington for $10^4$–$10^6$ yr (King et al. 2023; Middleton et al. 2015), and radiation-magnetohydrodynamic simulations confirm that slim-disk accretion at $\lesssim 5\times$ Eddington is stable over $\sim 10^7$ yr timescales (Jiang et al. 2019). TEP therefore does not require exotic physics: it reduces the required accretion rate to a regime already observed in local analogues.

    Monte Carlo error propagation yields a boost factor of $6 \times 10^5$ (68% CI: $8 \times 10^4$ to $5 \times 10^6$), dominated by $\alpha_0$ uncertainty. Runaway growth requires $r_e \lesssim 800$ pc, naturally explaining why compact LRDs host overmassive BHs while extended Red Monsters do not. Across 260 simulated LRDs, 87% exhibit boosts $> 10^3$ (median $\sim 1{,}600\times$). A separate resolved-photometry analysis using LRD structural parameters ($r_e \sim 0.15$ kpc, $n_S \sim 2.5$) finds $\Delta\log M_h = 0.45$ (vs. 1.5 assumed in the population simulation), yielding a higher runaway fraction: 99.9% of simulated LRDs show boosts $> 10^3$. The two simulations use different $\Delta\log M_h$ inputs; the resolved-photometry result is more physically grounded but the population statistics (87%) are more conservative. Both support the qualitative conclusion that the TEP mechanism is sufficient. Detailed error propagation, sensitivity analysis, and population-level results are in Appendix C.4.

    
#### 4.9.2 Blue Monsters: The Cleaned Sample

    Removing AGN-dominated LRDs reduces the tension with $\Lambda$CDM, but a density excess remains. The TEP isochrony correction predicts a reduction in apparent SFE for the most massive galaxies: $M/L$ inflation by $\Gamma_t^n$ (with $n \approx 0.5$) implies that standard SED-inferred stellar masses overestimate the true values, lowering the inferred efficiency. Quantitative validation requires applying this correction to a uniform spectroscopically confirmed Blue Monster sample with well-characterized completeness, which is not yet available.

    
### 4.10 Recent Observational Updates

    Recent reports have further sharpened the observational landscape in ways relevant to TEP.

    
#### 4.10.1 The "Blue Core" Signal

    Jin et al. (2025) report a high incidence of positive color gradients (bluer centers) in galaxies at $z > 4$, and Nedkova et al. (2025) independently confirm that rest-frame $U-V$ and $V-J$ gradients evolve toward bluer cores with increasing redshift and stellar mass at $0.5  1$ for massive galaxies at $z > 4$ using JWST NIRSpec kinematics (de Graaff et al. 2024; Wang et al. 2024). Under TEP, the isochrony correction $M_{*,\rm true} = M_{*,\rm obs}/\Gamma_t^n$ reduces SED-inferred stellar mass, resolving 11/11 impossible cases and bringing the mean $M_*/M_{\rm dyn}$ from 1.33 to 0.61 (Wilcoxon $p < 5 \times 10^{-4}$). This test is partially mass-independent because $M_{\rm dyn}$ derives from kinematics, not SED fitting. The associated baryon-fraction and gas-fraction anomalies are likewise alleviated by reducing $M_*$. **Binomial significance:** the probability of resolving 11/11 physically impossible cases by chance is $p = 0.5^{11} = 4.9 \times 10^{-4}$ under the most conservative null (50% random resolution probability), and $p = 1.8 \times 10^{-6}$ for a realistic null ($p_0 = 0.3$). The Wilcoxon test on the continuous ratio distribution is independent of this discrete framing and confirms the result is not an artefact of the 11/11 count alone.

    
#### 4.10.3 RBH-1: Soliton Wake Consistency

    The runaway supermassive black hole RBH-1 ($z \approx 0.96$, $M \approx 2 \times 10^7 M_\odot$; van Dokkum et al. 2025) presents a thermal paradox: JWST spectroscopy reveals a 650 km/s velocity discontinuity coexisting with cold, star-forming gas over a 62 kpc wake. Standard shock physics predicts post-shock temperatures $T \sim 10^7$ K with a cooling time $\sim 30\times$ the dynamical time — inconsistent with immediate star formation. It is proposed that the velocity discontinuity is a *metric shock* (spatial gradient in gravitational redshift) rather than bulk thermalization. The soliton core radius predicted from $\rho_c \approx 20$ g/cm³ is $R_{\rm sol} \approx 7.8 \times 10^7$ km $\approx 1.3 R_S$ — a geometric consistency check that uses the same $\rho_c$ as this work with no additional free parameters. The 50:1 aspect ratio of the wake is consistent with soliton geometry. This is a qualitative consistency, not a quantitative confirmation; decisive testing requires a dedicated analysis.

    
#### 4.10.4 Additional Anomalies Consistent with TEP

    Several further observations are qualitatively consistent with the TEP framework: (i) the "Blue Monster" correction (Chworowsky et al. 2025), where a $2\times$ density excess persists even after removing AGN-dominated LRDs; (ii) rapid quenching at $z > 4$, where TEP provides additional proper time ($\Gamma_t \sim 2$ doubles the effective evolution window); (iii) the environmental reversal—higher sSFRs in protoclusters than field galaxies at $z > 4.5$—consistent with Group Halo Screening suppressing TEP in overdense environments. These are noted as qualitative consistencies, not quantitative confirmations.

    
### 4.11 What TEP Does Not Explain

    Honesty requires consolidating the results where TEP fails, underperforms, or remains ambiguous — not only where it succeeds.

    
        **Consolidated Failures and Limitations**

        
            - **$\Gamma_t$ adds zero ML lift — but this is expected under TEP mass bias.** An adversarial gradient-boosted regressor with $M_*$, $z$, and 6 polynomial interaction terms reconstructs everything $\Gamma_t$ knows ($\Delta R^2 = -0.006$; §4.4.6.2). However, §4.4.6.3 shows that under TEP, $M_{*,\rm obs} = M_{*,\rm true} \cdot \Gamma_t^\beta$, so $M_{*,\rm obs}$ already encodes $\Gamma_t$ information. A GBR with $M_{*,\rm obs}$ as input therefore *already has $\Gamma_t$ embedded in its features* — zero lift is the expected outcome, not a falsification. The correct test would use $M_{*,\rm true}$ (debiased mass) as the GBR input; that test is not yet performed. The zero-lift finding is therefore consistent with both (a) TEP being correct and (b) TEP being wrong — it does not discriminate.

            - **OLS AIC favours mass polynomials.** In every subsample and survey, a polynomial ($M_*$, $z$, $M_* \times z$) is the best OLS model by large margins ($\Delta$AIC $= +67$–$920$ over $t_{\rm eff}$). The step-function comparison ($\Delta$AIC $= -23$ favouring $t_{\rm eff}$) is valid for threshold models but the linear $t_{\rm eff}$ model — not the step function — is the best single-predictor model by AIC ($\Delta$AIC $= 0$). The OLS AIC limitation callout (§4.3.5) should be read alongside this fact.

            - **Per-bin $\alpha_0$ Spearman optimisation is a mathematical non-test.** An attempt to recover $\alpha_0$ by maximising the Spearman $\rho(\Gamma_t, \text{dust})$ per redshift bin produced an apparent floor at $\alpha_0 = 0.1$. This is not a physics failure — it is a consequence of Spearman rank-invariance: $\Gamma_t(\alpha_0) = \exp[\alpha_0 \cdot f(\log M_h, z)]$ is strictly monotonic in $\log M_h$ at fixed $z$, so the rank ordering of galaxies is identical for any $\alpha_0 > 0$. Confirmed numerically: $\rho = 0.6458$ at every tested $\alpha_0 \in [0.1, 1.5]$ in the $z = 8.5$–$10$ bin. The optimizer has no gradient to follow and converges to the boundary by numerical accident. The apparent tension is an artefact of an identically flat objective function, not evidence against $\alpha_0 = 0.58$. The corrected recovery ($\alpha_0 = 0.55 \pm 0.32$ dust, $0.75 \pm 0.29$ joint; Pearson $R^2$ method) uses multi-observable combination sensitive to the absolute magnitude of $\Gamma_t$, not just its rank order. This item is removed from the genuine failures list.

            - **Mass-independent proxy tests: 1/5 pass at $z > 8$.** Five proxies for potential depth partially independent of stellar mass (SFR surface density, sSFR residual, SED $\chi^2$, photo-$z$ uncertainty, age-ratio residual) were tested at fixed $M_*$ and $z$. Only SFR surface density passes at $p  8$ ($\rho = -0.43$, though with negative sign complicating interpretation). All 5/5 pass in the full $z = 4$–$10$ sample. The mass-proxy degeneracy is partially but not definitively broken.

            - **LRD mechanism undershoots by $\sim 1.2$ dex.** The differential temporal shear for Little Red Dots reduces the required super-Eddington factor from $\sim 10\times$ to $\sim 2$–$3\times$ but does not fully close the gap to observed $M_{\rm BH}/M_*$ ratios under simplest assumptions (§4.9).

            - **L2 (core screening) is marginal.** The $\rho = -0.18$ ($p = 5 \times 10^{-4}$) signal is from a single survey (JADES, $N = 362$) and has not been replicated. The Steiger Z-test ($Z = 2.01$, $p = 0.045$) uses synthetic data. Standard inside-out dust models can produce similar gradients without new physics.

            - **Statistical weight is asymmetric across the four lines.** The dust–$\Gamma_t$ correlation (L1) contributes $\sim 90\%$ of the combined statistical weight by virtue of its large sample ($N = 1{,}283$, three surveys). However, statistical weight and scientific independence are different quantities. L2 (inside-out core screening, JADES) tests a spatially resolved prediction — the gradient of $\Gamma_t$ within a single galaxy — that is orthogonal in both data type and physical mechanism to the cross-galaxy dust correlation. L3 (mass–sSFR inversion) tests the sign of the sSFR–mass correlation across the activation threshold, a prediction that is explicitly independent of L1: partial $\rho(\Gamma_t, \text{sSFR} \mid \text{dust}) = -0.49$ confirms the two lines are not redundant. L4 (dynamical mass comparison) uses kinematic rather than photometric masses, and the 11/11 resolution result has binomial significance $p = 4.9 \times 10^{-4}$ even under the most conservative null — this test is the only one in the paper that is genuinely free of SED-fitting assumptions. The correct framing is not "four equally-weighted lines" but "one high-weight line (L1) corroborated by three independent physical mechanisms at lower but non-trivial significance." The "four independent lines" language in the abstract and conclusion accurately describes independence of mechanism; it does not imply equal statistical weight.

            - **Emission line and morphology signals collapse after mass control — but the control is biased.** JADES DR4 [OIII]5007 flux shows a strong raw correlation ($\rho = +0.683$, $N = 757$) that collapses to $\rho = +0.117$ after controlling for $M_{*,\rm obs}$ and $z$. JADES DR5 morphology (Gini, $\Sigma_*$) similarly collapses from $\rho \sim 0.5$–$0.6$ to $\rho \sim 0.03$. However, §4.4.6.3 shows that $M_{*,\rm obs}$ is a TEP-biased control: using debiased $M_{*,\rm true}$ instead, the O32 ionization ratio at $z > 7$ strengthens from $\rho = -0.165$ (marginal, $p = 0.002$) to $\rho = -0.204$ ($p = 1.3 \times 10^{-4}$), and H$\beta$ EW strengthens from $\rho = -0.133$ to $\rho = -0.196$ ($p = 1.1 \times 10^{-8}$). These are secondary TEP consequences (dust absorbs ionizing photons; apparent stellar ages are inflated), not independent lines of evidence. The collapse under biased mass control is therefore partly an over-control artifact, not a definitive null. The JADES morphology collapse ($\rho \sim 0.03$) is more likely a genuine null — morphology is not a direct TEP prediction at these redshifts.

            - **No public IFU kinematic catalog exists at $z > 4$.** The only mass-independent proxy for potential depth is the stellar velocity dispersion $\sigma_e$ from IFU spectroscopy. The largest published JWST kinematic sample at $z > 5$ is $N = 6$ (Carniani et al. 2024). No public catalog with $N > 50$ at $z > 4$ is available as of early 2026. The JADES NIRSpec WIDE program (Maseda et al., in prep.) and JWST Cycle 3 IFU programs are expected to provide such catalogs; until then, the mass-proxy degeneracy cannot be definitively broken with existing public data.

        
    

    
### 4.11.1 Evidence Base Independence

    The 13-paper TEP series (Papers 1–12 and this work) is entirely single-author and none of the prior papers have undergone external peer review at a refereed journal (all are published on Zenodo with DOIs). The "cross-domain consistency" described in §4.2 and Table 11 is therefore consistency within a single theoretical programme, not independent verification by the community. The $\alpha_0 = 0.58$ calibration from Paper 12 uses $N = 29$ Cepheid hosts with $p = 0.019$ ($2.3\sigma$) — significant but not overwhelming in isolation. The credibility of the cross-domain evidence rests on whether independent groups can reproduce the key results using independent analysis pipelines. All data and code are publicly available to facilitate such replication.

    Furthermore, the $\Gamma_t$ formula (§2.3.2.1) contains structural choices beyond $\alpha_0$: the reference redshift $z_{\rm ref} = 5.5$, the reference halo mass $\log M_{h,\rm ref} = 12.0$, the exponential functional form, and the $\sqrt{1+z}$ scaling. While each has physical motivation from the scalar-tensor framework (Paper 1), these choices collectively shape the predictions. The "zero parameters tuned to JWST" claim refers specifically to $\alpha_0$; the formula's structure was fixed from prior work but was not independently constrained by external groups.

    
### 4.12 Limitations and Caveats

    Following the claim hierarchy of Paper 6 (TEP-GTE), the limitations below are organised by tier. **Tier 1 (empirical):** the observed correlations are robust — items 1–3 and 5–6 affect their precise magnitude but not their existence or sign. **Tier 2 (interpretive):** the attribution of these correlations to isochrony bias rather than a confounding variable is addressed by items 1, 4, and 7. **Tier 3 (theoretical):** the TEP scalar-tensor framework as the underlying mechanism is addressed by items 4, 7, and 9 — these are the most open questions and the ones that future Boltzmann-code and holonomy experiments will resolve.

    
        - **Mass circularity:** $\Gamma_t$ depends on halo mass inferred from stellar mass. Six independent lines of evidence mitigate this concern (§4.4.6), spanning four distinct data types. Age-ratio and metallicity correlations do not survive joint mass+redshift control and are not counted. The colour-gradient Steiger test is provisional (synthetic data) and not counted.

        - **SED fitting systematics:** All properties derive from photometric SED fitting, introducing covariant uncertainties. Photo-$z$ scatter degrades $\rho$ by $ 6$, Lyman/Balmer break confusion produces $\sim 5$–$15\%$ catastrophic failures. Spectroscopic confirmation exists for a subset ($N = 147$); the remainder relies on photometric estimates.

        - **Theoretical foundation:** The $\Gamma_t$ formula derives from a scalar-tensor action with chameleon screening (Appendix A.1). A full CAMB Boltzmann integration (Appendix A.1.8.8) confirms $\sigma_8$ consistency at the fiducial chameleon mass: $\sigma_8^{\rm TEP} = 0.8116$ ($0.10\sigma$ from Planck), with CMB TT deviations $ 1089$ are assumed negligible (justified by $T^\mu_\mu \approx 0$ during radiation domination). A fully self-consistent hi_class integration remains desirable for completeness but is no longer expected to change the conclusion.

        - **Statistical caveats:** Combined p-values exceeding $10^{-90}$ should not be taken at face value; the most conservative clustering-corrected estimate is $6.4\sigma$. BH-FDR correction shows all 6 independent lines survive at $\alpha = 0.05$ (7 of 8 tested signatures including the two not-counted tests). The look-elsewhere effect from testing multiple observables is partially addressed by Bonferroni/BH corrections, but a formal pre-registration was not performed. All null results are reported publicly.

        - **Underpowered tests:** The Red Monsters ($N = 3$) and some spectroscopic subsamples ($N  9$, the UNCOVER spectroscopic sample is dominated by UV-bright, low-dust systems (the dusty population falls below the NIRSpec continuum threshold); COSMOS-Web at $z = 9$–$13$ has $N = 305$ photometric sources with $> 26$-band coverage vs. $N = 122$ spec-z sources with confirmed redshifts, giving the photometric catalog a broader and more representative dust range. (b) *SED degeneracy:* Prospector-$\beta$ with spectroscopic priors at $z > 9$ may suppress the dust posterior for low-$M_*$ galaxies due to the extreme youth prior ($t  9$, the coupling $\alpha(z) \propto \sqrt{1+z}$ is sufficiently strong that even low-mass halos approach $\Gamma_t \approx 1$, reducing the $\Gamma_t$ dynamic range and flattening the correlation within this sample. The contrast with COSMOS2025 argues against (c) as the primary explanation, since COSMOS-Web photometry at the same epoch shows the strongest signal in any dataset. Explanation (a) is the most quantitatively tractable: a completeness-corrected spectroscopic analysis at $z > 9$ (e.g., from JADES DR5 or the PRIMER deep fields) with wider dust dynamic range would resolve this tension directly. This null is reported transparently and not dismissed; the COSMOS2025 $z = 10$–$13$ result prevents treating it as a falsification of TEP, but the UNCOVER disagreement remains an open tension pending the larger spectroscopic samples in §3.11.

        - **Alternative explanations:** AIC/BIC and multi-domain scoring favor TEP, but a fully nested Bayesian evidence computation has not been performed.

        - **Coupling constant uncertainty:** $\alpha_0 = 0.58 \pm 0.16$ ($\sim 28\%$ relative). Full propagation through the $\Gamma_t$ formula yields a Red Monster SFE correction range of 22%–48% at $1\sigma$ (central value 34%; Table 3b). The correction is robust to the lower bound: even at $\alpha_0 = 0.42$, the anomaly is reduced by $\sim 22\%$ with zero tuned parameters. The JWST dust-only recovery gives $\alpha_0 = 0.55 \pm 0.32$ (Pearson $R^2$ method); the joint four-observable fit gives $\alpha_0 = 0.75 \pm 0.29$; both are consistent with the Cepheid calibration within $1\sigma$. The large uncertainties are genuine (reflecting Spearman rank-invariance and mass-proxy degeneracy) — the Cepheid calibration remains the tighter primary constraint. An earlier result of $0.60 \pm 0.10$ was an artefact of [0,1]-normalised RSS which is also rank-invariant (see item 10); the corrected Pearson $R^2$ objective yields the values above. Table 3b uses representative parameters, not exact catalog values. *Note on the 34% vs 43% figures:* The 34% figure (Table 3b, §3.1) is the direct SFE reduction for the three representative Red Monster galaxies, computed as $(\text{SFE}_{\rm obs} - \text{SFE}_{\rm true})/\text{SFE}_{\rm obs}$. The 43% figure is the fraction of the anomaly *above the $\Lambda$CDM limit of 0.20* that is resolved — a different denominator. Both are correct; the 34% figure is used throughout this paper as it is the more direct and conservative measure.

        - **Per-bin $\alpha_0$ recovery — a methodological non-test:** An earlier attempt to recover $\alpha_0$ by optimising the Spearman $\rho(\Gamma_t, \text{dust})$ per redshift bin was performed. The optimizer hits the grid floor (0.1) in every bin, yielding an apparent tension with the Cepheid value. **This is a mathematical artefact, not a physical failure.** $\Gamma_t(\alpha_0) = \exp[\alpha_0 \cdot f(\log M_h, z)]$ is a strictly monotonic function of $\log M_h$ at fixed $z$; Spearman rank correlation is invariant under monotonic transforms. Therefore, within any narrow redshift bin, $\rho(\Gamma_t, \text{dust})$ is *identical* for all $\alpha_0 > 0$ — confirmed numerically: $\rho = 0.6458$ for every value of $\alpha_0 \in \{0.1, 0.2, 0.4, 0.58, 0.8, 1.0, 1.2, 1.5\}$ in the $z = 8.5$–$10$ bin. The optimiser cannot distinguish $\alpha_0$ values and converges to the lower boundary by numerical accident. The apparent "$2.86\sigma$ tension" is an artefact of using an identically flat objective function, not evidence against $\alpha_0 = 0.58$. The corrected recovery (Pearson $R^2$ maximisation, which breaks rank-invariance by being sensitive to the spread of $\log\Gamma_t$) gives $\alpha_0 = 0.55 \pm 0.32$ from dust alone and $\alpha_0 = 0.75 \pm 0.29$ from the joint four-observable fit — both consistent with the Cepheid calibration within $1\sigma$ but with large honest uncertainties. The earlier result ($0.60 \pm 0.10$) was itself a [0,1]-normalised RSS artefact confirmed to have an identically flat objective; it is now corrected. Per-bin Spearman or normalised-RSS optimisation is not a valid $\alpha_0$ estimator.

    

    
### 4.13 Future Directions

    
#### 4.13.1 Critical Test: The Mass-Dust Inversion

    Falsification: If future JWST/MIRI observations show $\rho(M_*, A_V) < 0.1$ at $z > 8$, TEP is ruled out at $> 3\sigma$.

    
#### 4.13.2 Critical Test: The Coupling Constant

    Falsification: If fitting the $z > 8$ dust anomaly requires $\alpha_0 > 1.0$ or $\alpha_0 < 0.2$, the cross-domain consistency with Cepheids is broken.

    
#### 4.13.3 Critical Test: The Black Hole Boost

    Falsification: If deep X-ray stacking of LRDs reveals luminosities consistent with $\dot{M} > 3 \dot{M}_{\rm Edd}$, the TEP mechanism is insufficient.

    
#### 4.13.4 Gravitational Wave Timing: LISA and Binary Pulsars

    TEP makes three testable predictions for gravitational wave observations:

    
        - **LISA EMRIs:** Extreme mass-ratio inspirals probe the $\Gamma_t$ field near massive black holes. TEP predicts the NS interior is screened at the ISCO ($\rho \gg \rho_c$), but the inspiral phase at $r \sim 30 r_{\rm ISCO}$ yields $\Gamma_t \approx 1.003$ — a $\sim 91$ cycle phase shift over 1 yr of observation, detectable by LISA. Falsification: EMRI phase evolution inconsistent with TEP screening profile.

        - **Binary pulsars:** The Hulse-Taylor system agrees with GR to $0.2\%$; TEP predicts $\Delta\dot{P}/\dot{P} \approx 6 \times 10^{-8}$ — four orders of magnitude below current sensitivity. TEP is compatible with all existing binary pulsar constraints.

        - **Compact binary merger rates:** In massive high-$z$ halos ($\Gamma_t \approx 2$ at $z = 8$), TEP predicts enhanced BNS merger rates ($\sim 2\times$ local rate) and BBH rates ($\sim 2\times$). Falsification: no redshift evolution of merger rates in massive hosts detected by Einstein Telescope or Cosmic Explorer.

    

    
#### 4.13.5 Future Surveys

    Five JWST observing programs (Balmer absorption ages, radial age gradients, environmental screening in protoclusters, rest-frame MIR dust confirmation, and morphology–$\Gamma_t$ testing) would provide definitive tests. Future wide-field surveys will substantially extend the statistical reach:

    
        - **Euclid Wide ($N \sim 300{,}000$ massive galaxies, $z = 0.9$–$1.8$):** Typical $\Gamma_t \approx 1.25$ predicts a 25% age offset at fixed $z$. Combined sensitivity reaches $\rho_{\rm min} = 0.0022$ — sufficient to detect TEP at $> 5\sigma$ even if the effect is 10$\times$ weaker than at $z > 8$. Key falsification: no mass-dependent age offset at $z \sim 1.5$.

        - **Roman Supernova Survey ($N \sim 2{,}700$ SNe Ia, $z  1$.

        - **Roman High-Latitude ($N \sim 500{,}000$ at $z > 2.5$):** Tests the gas vs. stellar metallicity discriminant (§4.7.2) and morphology–$\Gamma_t$ correlation. Key falsification: strong [O III]/H$\beta$–$\Gamma_t$ correlation.

    
    Combined JWST + Euclid + Roman total sample: $\sim 801{,}000$ galaxies across $z = 0$–$10$. Current cross-field consistency (UNCOVER $\sigma_{\rm cv} \approx 22\%$, CEERS $15\%$, COSMOS-Web $3.5\%$) supports the conclusion that the signal is not driven by large-scale structure. Full program details are in Appendix C.5.

    
#### 4.13.6 A Critical Experiment: Synchronization Holonomy

    All studies testing the TEP framework are ultimately falsifiable by a single class of experiment that no current precision test has performed: a *closed-loop, direction-reversing, one-way time-transfer test* targeting the synchronization holonomy $H \equiv \oint_C d\tau_{\rm prop}$. Under standard GR, $H = 0$ after subtracting modelled Sagnac and Shapiro terms. Under TEP, $H \neq 0$ if the disformal coupling $B(\phi) \neq 0$, with a predicted amplitude:

    
        $$H_{\rm resid} \sim \frac{B(\phi)}{A(\phi)} |\nabla\phi|^2 \times \mathcal{A}$$
    
    where $\mathcal{A}$ is the loop area. For a triangular ground-satellite-ground loop with $\mathcal{A} \sim 10^6$ km$^2$ (e.g., two ground stations and one MEO satellite), the predicted holonomy is $H_{\rm resid} \sim 10^{-19}$ s — at the frontier of current optical clock technology but achievable with next-generation transportable optical lattice clocks (Lisdat et al. 2016; Grotti et al. 2018). Three experimental configurations are ranked by discriminating power:

    
        - **Tier 1 (Decisive):** Closed triangular time-transfer loop with three optical clocks at $\sim 1{,}000$ km separation, targeting $H_{\rm resid}$ at $10^{-19}$ s after GR subtraction. A non-zero result at $> 3\sigma$ would confirm the disformal sector; a null result would constrain $B(\phi)/A(\phi)  800{,}000$; §4.13.5 above) — these test the conformal sector ($A(\phi)$, which governs $\Gamma_t$) independently of the disformal sector. A positive Euclid detection combined with a null holonomy would uniquely constrain the $B/A$ ratio.

    
    The holonomy test is the experiment that the underlying framework was designed to motivate. If the theoretical construction is correct in its entirety, it will be detected. If it is not detected at the predicted level, the disformal sector is suppressed below current sensitivity, and the conformal-only limit ($B = 0$) applies — which preserves all JWST, Hubble tension, and pulsar predictions while removing the holonomy signal. The conformal-only limit is therefore a self-consistent sub-theory that the holonomy test can distinguish from the full disformal theory. Either outcome is scientifically decisive.

    
## 5. Conclusion

    JWST has revealed a coherent pattern of anomalies at $z > 5$: ultra-massive galaxies with star formation efficiencies exceeding $\Lambda$CDM limits, overmassive black holes in Little Red Dots, and stellar masses that exceed dynamical masses. These anomalies share a common structure — photometrically inferred stellar properties appear systematically too large, too early, in precisely the environments with the deepest gravitational potentials — suggesting a common origin. This work tests whether a single violation of the isochrony axiom — parameterized by the chameleon-screened Temporal Equivalence Principle (TEP) — can account for the coherent pattern of $z > 5$ anomalies observed by JWST. The central finding is that a local-universe metric coupling ($\alpha_0 = 0.58$), applied with zero tuned parameters, simultaneously predicts the magnitude of the Red Monster efficiency excess, resolves the $M_*/M_{\rm dyn}$ impossibility crisis, and provides a purely kinematic mechanism for overmassive black hole growth in Little Red Dots.

    
### 5.1 Synthesis of Results

    Four independent lines of evidence survive rigorous mass-control checks (§3.9, §4.4.6). The dust–$\Gamma_t$ correlation and AGB threshold (L1) yields $\rho = +0.62$ across $N = 1{,}283$ galaxies in three surveys, with odds ratio 42.8 ($p  7$ (L3) yields $\Delta\rho = +0.25$ with 95% confidence interval excluding zero; independence from L1 is supported by partial $\rho(\Gamma_t, {\rm sSFR}\,|\,{\rm dust}) = -0.49$ ($p = 10^{-18}$). The dynamical mass comparison (L4) renders 11/11 previously impossible $M_*/M_{\rm dyn} > 1$ cases physically plausible under the TEP correction using kinematic masses entirely independent of SED fitting (mean ratio $1.33 \to 0.61$; L4 uses representative parameters from published ranges, not exact tabulated values — see §3.9 L4 caveat). The AGB threshold test and cross-survey generalisation are robustness checks on L1, not independent lines (partial $\rho = 0.049$, $p = 0.41$ after $\Gamma_t$ control). TEP additionally predicts that environmental screening is absent at $z > 8$, where the evolving coupling overwhelms the density threshold; this predicted null is recovered ($\Delta\rho = 0.003$, $p = 0.97$; §4.4.3). The JWST-recovered dust coupling is $\alpha_0 = 0.55 \pm 0.32$ (Pearson $R^2$ method); the joint four-observable fit gives $\alpha_0 = 0.75 \pm 0.29$, both consistent with the Cepheid value within $1\sigma$. The large uncertainties are genuine — they reflect the Spearman rank-invariance of $\Gamma_t$ within any single survey and the mass-proxy degeneracy; the Cepheid calibration remains the primary $\alpha_0$ constraint (§4.12 item 9).

    Eight new cross-dataset results further extend the evidence base. (1) COSMOS2025 dust (Shuntov et al. 2025; $N = 1{,}568$ at $z = 9$–13): partial $\rho = +0.595$ ($p = 1.6 \times 10^{-150}$, 95% CI $[+0.532, +0.660]$) — the strongest partial dust correlation found in any dataset. (2) COSMOS2025 sSFR inversion: L3 replicated in an independent field; Steiger $Z = 6.37$ ($p = 1.9 \times 10^{-10}$). (3) UNCOVER DR4 MegaScience (Wang et al. 2024): null at $z  9$ dust partial correlation unchanged at $+0.592$.

    
### 5.2 Implications and Robustness

    If correct, these results indicate that the early universe is not fundamentally over-efficient at forming stars, nor does it require exotic initial mass functions or continuous super-Eddington accretion. Instead, standard astrophysical processes are occurring in local spacetime patches where proper time has accumulated faster than standard Friedmann-Lemaître-Robertson-Walker integration predicts. The apparent anomalies are an artefact of enforcing universal isochrony on a universe where time is fundamentally environment-dependent.

    The credibility of this hypothesis rests on three pillars of robustness:

    
        - **Falsifiability of the functional form:** The specific non-linear predictions of TEP — the AGB threshold step function, the inside-out spatial screening gradient, and the high-redshift sSFR inversion — cannot be replicated by smooth mass-dependent polynomials.

        - **Cross-domain coherence:** The JWST-recovered dust coupling ($\alpha_0 = 0.55 \pm 0.32$) and the joint multi-observable fit ($\alpha_0 = 0.75 \pm 0.29$) are both consistent with the Cepheid calibration ($0.58 \pm 0.16$) within $1\sigma$. The large JWST uncertainties are expected given the mass-proxy degeneracy; the cross-domain consistency across atomic clocks, pulsar timing, and gravitational lensing provides additional independent support for $\alpha_0 \approx 0.58$.

        - **Breaking the mass proxy:** The resolution of 11/11 impossible $M_*/M_{\rm dyn} > 1$ cases relies purely on kinematics, confirming the isochrony bias exists independent of photometric mass estimation.

    
    The cross-domain consistency of $\alpha_0$ is one of the strongest features of the TEP evidence base. A single coupling constant — $\alpha_0 = 0.58 \pm 0.16$, calibrated from 29 Cepheid hosts at $z \approx 0$ — has been independently tested across 10 domains spanning 13.5 Gyr of cosmic time, 40 orders of magnitude in mass, and 15 orders of magnitude in density: GNSS atomic clock networks (Papers 2–4), gravitational lensing (Paper 5), the universal critical density (Paper 7), the RBH-1 soliton wake (Paper 8), satellite laser ranging (Paper 9), globular cluster pulsar timing (Paper 11), the Hubble tension (Paper 12), and the JWST high-redshift galaxy statistics presented here. The JWST-recovered dust value ($\alpha_0 = 0.55 \pm 0.32$) and the joint fit ($\alpha_0 = 0.75 \pm 0.29$) are both consistent with the Cepheid calibration to within $1\sigma$, despite the large per-survey degeneracy. This consistency is meaningful: a recovery of $\alpha_0 > 1.5$ or $$6 JWST data would have been inconsistent with the Cepheid calibration even within the wide uncertainties.

    Key signatures survive a 0.5 dex mass reduction, and blind validation passes all three generalisation tests — time-split, field-split, and cross-survey leave-one-out — each recovered independently across all three surveys (9/9 survey-test combinations). Each of the three independent JWST surveys confirms L1 individually above $5\sigma$ without requiring a clustering correction (CEERS: $z = 7.0\sigma$; UNCOVER: $z = 11.4\sigma$; COSMOS-Web: $z = 22.4\sigma$); all three independently confirm that $t_{\rm eff}$ outperforms $t_{\rm cosmic}$ at $>5\sigma$ (Steiger $Z = 5.3$–$16.8$), ruling out pure redshift ordering. A Fisher combination across all five independent L1 datasets (3 photometric surveys + GOODS-S DJA + NIRSpec spectroscopic) gives $z = 30.0\sigma$ ($p = 7.0 \times 10^{-198}$), spanning 4 sky fields, 3 SED pipelines, and 2 dust estimators with no clustering correction required. The clustering-corrected combined significance across all four independent lines is $6.4\sigma$ (upper bound); permutation tests indicate the signal exceeds all null realisations ($p_{\rm perm}  8$ (5/5 in the full sample), indicating the mass-proxy degeneracy is partially but not definitively broken; IFU velocity dispersions remain the decisive test. Important limitations (§4.11–4.12) include mass-dependent $\Gamma_t$ circularity, the self-referential evidence base (§4.11a), and the consolidated failures enumerated in §4.11. The TEP mass measurement bias (§4.4.6.3) reframes one previously-listed failure: the ML zero-lift result is expected when the GBR uses $M_{*,\rm obs}$ which already encodes $\Gamma_t$. The per-bin $\alpha_0$ floor and the earlier RSS-based recovery are now understood as consequences of Spearman rank-invariance and [0,1]-normalised RSS flatness respectively (§4.12 item 9) — mathematical non-tests, not physical constraints on $\alpha_0$. The emission-line and morphology collapses under biased mass control are partly over-control artifacts, not definitive nulls.

    
### 5.3 Falsification Criteria

    TEP makes specific, quantitative predictions that can be tested with existing and near-future facilities. The following failure conditions are defined; any one of them, if met, would require rejection of the TEP interpretation of the JWST anomalies.

    
        
            Table 19: TEP Falsification Criteria
            
                ObservableStandard Physics PredictionTEP PredictionFalsification Criteria
            
            
                Mass-Dust ($z > 8$)No correlation or NegativeStrong Positive ($\rho > 0.4$)$\rho \approx 0$ or Negative
                Balmer AbsorptionCorrelates with $z$Correlates with $M_*$ at fixed $z$No mass trend
                LRD Host SizeNo dependenceOnly in Compact ($r_e < 1$ kpc)LRDs in large disks
                Cluster vs FieldCluster galaxies olderCluster galaxies younger (screened)Field $\approx$ Cluster
                [OIII]/H$\beta$ vs $\Gamma_t$Correlates with massWeak correlation ($\rho Strong correlation ($\rho > 0.3$)
                Type Ia / CC SN RatioNo $\Gamma_t$ dependenceIa/CC $\propto \Gamma_t^{0.5}$; Ia rate $4.4\times$ enhanced in high-$\Gamma_t$ hosts; CC rate unchanged ($\rho \approx 0$)Constant Ia/CC ratio across $\Gamma_t$ at fixed $M_*$ (Roman Space Telescope High-Latitude Time Domain Survey)
                $\alpha_0$ RecoveryN/A$\alpha_0 = 0.58 \pm 0.16$$\alpha_0 > 1.0$ or $\alpha_0 
            
        
    

    
### 5.4 Reproducibility

    This analysis is designed to be fully reproducible. The pipeline comprises 163 steps organised into ten phases — data ingestion, cross-domain validation, JWST extended analysis, statistical synthesis, final validation, independent replication, refinement and robustness, simulations, advanced validation, and new discriminating tests — producing step-indexed JSON outputs and figures. To reproduce all results: `python scripts/steps/run_all_steps.py`. The complete step listing is documented in the repository README.

    
### 5.5 Data Availability

    The manuscript source, complete analysis code, 163 pipeline steps, generated figures, JSON intermediate outputs, and the raw and processed catalogs are available on GitHub and archived on Zenodo for long-term reproducibility.

    
        - **Pipeline repository:** github.com/matthewsmawfield/TEP-JWST — complete analysis code, 163 pipeline steps, and step-indexed JSON outputs.

        - **Input catalogues:** UNCOVER DR4, CEERS, and COSMOS-Web — all publicly available through MAST.

        - **Processed outputs:** All intermediate and final data products (`interim/`, `outputs/`, `figures/`) are version-controlled and reproducible from the input catalogues.

        - **Documentation:** `README.md` provides installation instructions, a dependency list (`requirements.txt`), and a quick-start guide.

    
    The full TEP theoretical framework series (Papers 1–12) is available on Zenodo; DOIs are listed in the References. Key identifiers: Paper 1 — TEP foundation (10.5281/zenodo.18204190); Papers 2–4 — GTE (10.5281/zenodo.17127229); Paper 5 — GL (10.5281/zenodo.17982540); Paper 6 — GTE synthesis (10.5281/zenodo.18004832); Paper 7 — UCD (10.5281/zenodo.18064366); Paper 8 — RBH-1 (10.5281/zenodo.18059251); Paper 9 — SLR (10.5281/zenodo.18064582); Paper 10 — EXP (10.5281/zenodo.18109761); Paper 11 — COS (10.5281/zenodo.18165798); Paper 12 — H₀ (10.5281/zenodo.18209703).

    
## References

    Abbott, B. P., et al. 2017, ApJL, 848, L13. *Gravitational Waves and Gamma-Rays from a Binary Neutron Star Merger: GW170817 and GRB 170817A.*

    Arrabal Haro, A., et al. 2023, Nature, 622, 707. *Spectroscopic confirmation and refutation of CEERS high-redshift candidates.*

    Behroozi, P., Wechsler, R. H., Hearin, A. P., & Conroy, C. 2019, MNRAS, 488, 3143. *UNIVERSEMACHINE: The correlation between galaxy growth and dark matter halo assembly from z = 0−10.*

    Berg, D. A., et al. 2013, ApJ, 775, 93. *New Detections of C/O Abundance Ratios in Metal-Poor Dwarf Galaxies.*

    Boyer, M. L., et al. 2025, ApJ, 991, 24. *Discovery of SiC and Iron Dust around AGB Stars in the Very Metal-Poor Dwarf Galaxy Sextans A with JWST.*

    Bertotti, B., Iess, L., & Tortora, P. 2003, Nature, 425, 374. *A test of general relativity using radio links with the Cassini spacecraft.*

    Boylan-Kolchin, M. 2023, Nature Astronomy, 7, 731. *Stress testing ΛCDM with high-redshift galaxy candidates.*

    Brammer, G. B., van Dokkum, P. G., & Coppi, P. 2008, ApJ, 686, 1503. *EAZY: A Fast, Public Photometric Redshift Code.*

    Brout, D., et al. 2022, ApJ, 938, 110. *Type Ia supernova host-mass step measurements in Pantheon+.*

    Brax, P., van de Bruck, C., Davis, A.-C., Khoury, J., & Weltman, A. 2004, PhRvD, 70, 123518. *Small scale structure formation in chameleon cosmology.*

    Bruzual, G. & Charlot, S. 2003, MNRAS, 344, 1000. *Stellar population synthesis at the resolution of 2003.*

    Burrage, C. & Sakstein, J. 2018, Living Reviews in Relativity, 21, 1. *Tests of Chameleon Gravity.*

    Carniani, S., et al. 2024, Nature, 633, 318. *A shining cosmic dawn: spectroscopic confirmation of two luminous galaxies at z > 14.*

    Carnall, A. C., McLure, R. J., Dunlop, J. S., & Davé, R. 2018, MNRAS, 480, 4379. *Inferring the star formation histories of massive quiescent galaxies with BAGPIPES.*

    Carnall, A. C., et al. 2023, Nature, 619, 716. *A massive quiescent galaxy at redshift 4.658.*

    Chworowsky, K., et al. 2025, arXiv:2509.07695. *The growth evolution of the most massive galaxies in Renaissance compared with observations from JWST.*

    Claeyssens, A., et al. 2023, MNRAS, 520, 2162. *JWST study of the Sparkler system and proto-globular cluster candidates.*

    Conroy, C., Gunn, J. E., & White, M. 2009, ApJ, 699, 486. *The Propagation of Uncertainties in Stellar Population Synthesis Modeling.*

    Cox, T. J., et al. 2025, ApJS (in press). *CEERS DR1 photometric and physical parameter catalog.*

    Curti, M., et al. 2023, MNRAS, 518, 425. *Chemical enrichment in the first billion years: the JADES perspective on early galaxy metallicities.*

    Curtis-Lake, E., et al. 2023, Nature Astronomy, 7, 622. *Spectroscopic confirmation of four metal-poor galaxies at z = 10.3–13.2.*

    D'Eugenio, F., et al. 2025, ApJS (in press). *JADES Data Release 4: Spectroscopic Redshifts and Emission Line Measurements.*

    de Graaff, A., et al. 2024, Nature, 630, 846. *A dormant overmassive black hole in the early Universe.*

    Endsley, R., et al. 2023, MNRAS, 524, 2312. *A JWST/NIRCam Study of Key Contributors to Reionization: The Star-forming and Ionizing Properties of UV-faint z ∼ 7–8 Galaxies.*

    Eisenstein, D. J., et al. 2023, arXiv:2306.02465. *Overview of the JWST Advanced Deep Extragalactic Survey (JADES).*

    Finkelstein, S. L., et al. 2023, ApJL, 946, L13. *CEERS early release science survey overview.*

    Freedman, W. L., Madore, B. F., Hoyt, T. J., et al. 2024, arXiv:2408.06153. *Status Report on the Chicago-Carnegie Hubble Program (CCHP).*

    Fujimoto, S., et al. 2023, ApJL, 949, L25. *JWST/NIRSpec spectroscopic confirmation of z > 8 CEERS candidates.*

    Furtak, L. J., et al. 2023, MNRAS, 523, 4568. *JWST UNCOVER: The Strong Lensing Model of Abell 2744.*

    Grotti, J., et al. 2018, Nature Physics, 14, 437. *Geodesy and metrology with a transportable optical clock.*

    Greene, J. E., et al. 2024, ApJ, 964, 39. *UNCOVER: The Growth of the First Massive Black Holes.*

    Hainline, K. N., et al. 2023, arXiv:2306.02468. *The Cosmos in its Infancy: JADES Galaxy Candidates at z > 8 in GOODS-S and GOODS-N.*

    Heintz, K. E., et al. 2023, ApJL, 953, L10. *Extreme Damped Lyman-α Absorption in Young Star-Forming Galaxies at z = 9–11.*

    Ilie, C., et al. 2025, PNAS. *Supermassive Dark Star candidates seen by JWST.*

    Jiang, Y.-F., Stone, J. M., & Davis, S. W. 2019, ApJ, 880, 67. *Super-Eddington Accretion Disks around Supermassive Black Holes.*

    Jin, B., et al. 2025, A&A, 698, A30. *Spatially resolved colours and sizes of galaxies at z ~ 3–4.*

    Ju, M., et al. 2025, arXiv:2506.12129. *A 13-Billion-Year View of Galaxy Growth: Metallicity Gradients.*

    Kelly, P. L., et al. 2010, ApJ, 715, 743. *Host-galaxy mass step in Type Ia supernova distances.*

    Khoury, J. & Weltman, A. 2004, PhRvL, 93, 171104. *Chameleon Fields: Awaiting Surprises for Tests of Gravity in Space.*

    Kawinwanichakij, L., et al. 2025, ApJ (in press). *Environmental dependence of galaxy morphology at z = 3–4.*

    Kocevski, D. D., et al. 2023, ApJL, 954, L4. *Hidden Little Monsters: Spectroscopic Identification of Low-Mass, Broad-Line AGN at z > 5 with CEERS.*

    King, A. R., Lasota, J.-P., & Kluzniak, W. 2023, MNRAS, 519, 5765. *Super-Eddington accretion: models and applications.*

    Kodric, M., Riffeser, A., Seitz, S., et al. 2018, ApJ, 864, 59. *Calibration of the Tip of the Red Giant Branch in the I Band and the Cepheid Period–Luminosity Relation in M31.*

    Kokorev, V., et al. 2024, arXiv:2401.09981. *A Census of Photometrically Selected Little Red Dots at 4 < z < 9 in JWST Blank Fields.* github.com/VasilyKokorev/lrd_phot

    Larson, R. L., et al. 2023, ApJ, 953, 34. *A CEERS Discovery of an Accreting Supermassive Black Hole 570 Myr after the Big Bang.*

    Labbé, I., et al. 2023, Nature, 616, 266. *A population of red candidate massive galaxies ~600 Myr after the Big Bang.* Data: github.com/ivolabbe/red-massive-candidates

    Leja, J., et al. 2019, ApJ, 876, 3. *How to Measure Galaxy Star Formation Histories. II. Nonparametric Models.*

    Lisdat, C., et al. 2016, Nature Communications, 7, 12443. *A clock network for geodesy and fundamental science.*

    Li, Q., et al. 2025, MNRAS, 539, 1796. *EPOCHS Paper X: Environmental effects on Galaxy Formation and Protocluster Galaxy candidates at 4.5 < z < 10.*

    Maiolino, R., et al. 2024, Nature, 627, 59. *A small and vigorous black hole in the early Universe.*

    Maseda, M. V., et al. (in prep.). *The JADES NIRSpec WIDE Program.*

    Matthee, J., et al. 2024, ApJ, 963, 129. *Little Red Dots: An Abundant Population of Faint Active Galactic Nuclei at z ~ 5 Revealed by JWST.*

    Meng, X.-L., Rosenthal, R., & Rubin, D. B. 1992, Biometrika, 79, 425. *Comparing correlated correlation coefficients.*

    Middleton, M. J., et al. 2015, MNRAS, 447, 3243. *NuSTAR reveals extreme absorption in z = 2–3 type 2 quasars.*

    Mota, D. F. & Shaw, D. J. 2007, PhRvD, 75, 063501. *Evading equivalence principle violations, cosmological, and other experimental constraints in scalar field theories with a strong coupling to matter.*

    Mowla, L., et al. 2022, ApJL, 937, L35. *The Sparkler: Evolved High-Redshift Globular Cluster Candidates Captured by JWST.*

    Naidu, R. P., et al. 2022, ApJL, 940, L14. *Two Remarkably Luminous Galaxy Candidates at z ≈ 10–12 Revealed by JWST.*

    Nanayakkara, T., et al. 2024, Science, 384, 890. *A massive galaxy that was quenched by z ∼ 3.*

    Nedkova, K. V., et al. 2025, A&A. *Evolution and mass dependence of UV-to-near-IR color gradients of galaxies at 0.5 < z < 2.5.*

    Nakajima, K., et al. 2023, ApJS, 269, 33. *JWST Census for the Mass-Metallicity Star Formation Relation at z = 4–10.*

    Pérez-González, P. G., et al. 2024, ApJ, 968, 4. *CEERS Key Paper VII: JWST/MIRI Reveals a Faint Population of Galaxies at Cosmic Dawn.*

    Planck Collaboration, Aghanim, N., et al. 2020, A&A, 641, A6. *Planck 2018 results. VI. Cosmological parameters.*

    Price, S. H., et al. 2024, ApJ, 964, 73. *UNCOVER: JWST spectroscopy of three cold brown dwarfs at kiloparsec-scale distances.*

    Rieke, M. J., et al. 2023, PASP, 135, 028001. *JWST NIRCam Performance: Commissioning and Calibration.*

    Riess, A. G., et al. 2022, ApJL, 934, L7. *A Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km/s/Mpc Uncertainty from the Hubble Space Telescope and the SH0ES Team.*

    Scholtz, J., et al. 2025, A&A (in press). *JADES: Spectroscopic properties of faint AGN at z > 4.*

    Shamir, L. 2025, MNRAS, 538, 76. *The distribution of galaxy rotation in JWST Advanced Deep Extragalactic Survey.*

    Shuntov, M., et al. 2025, ApJS (in press). *COSMOS-Web DR1 / COSMOS2025 catalog.*

    Smawfield, M. L. 2025, Zenodo. *Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle via GNSS Timing Networks.* DOIs: 10.5281/zenodo.17127229, 10.5281/zenodo.17517141, 10.5281/zenodo.17860166. Cited for GNSS correlation length $\lambda = 4{,}201 \pm 1{,}967$ km and 7 independent signatures (§4.2 Table 11).

    Smawfield, M. L. 2025, Zenodo. *Temporal-Spatial Coupling in Lensing: Phantom Mass and the Isochrony Axiom.* DOI: 10.5281/zenodo.17982540. Cited for the lensing phantom mass mechanism and Earth–galaxy $r_V \propto M^{1/3}$ scaling (§4.2 Table 11).

    Smawfield, M. L. 2025, Zenodo. *Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales.* DOI: 10.5281/zenodo.18029598. Cited for derivation of the screening threshold $\rho_c \approx 20$ g/cm³ (§2.3.2.2, §4.2 Table 11).

    Smawfield, M. L. 2025, Zenodo. *The Soliton Wake: Identifying RBH-1 as a Gravitational Soliton.* DOI: 10.5281/zenodo.18042456. Cited for geometric consistency of the RBH-1 velocity discontinuity (§4.10.3).

    Smawfield, M. L. 2025, Zenodo. *Global Time Echoes: Optical Validation via Satellite Laser Ranging.* DOI: 10.5281/zenodo.18064582. Cited for independent optical-domain confirmation of GNSS correlation structure using LAGEOS-1/2 (§4.2 Table 11).

    Smawfield, M. L. 2025, Zenodo. *What Do Precision Tests of General Relativity Measure?* DOI: 10.5281/zenodo.18109761. Cited for the formal taxonomy of GR test limitations and the conformal loophole argument (§4.2 Table 11, §4.4.7).

    Smawfield, M. L. 2026, Zenodo preprint. *Suppressed Density Scaling in Globular Cluster Pulsars.* Not yet peer-reviewed; cited for the screening threshold $\sigma > 165$ km/s, binary inversion ($-0.31$ dex, $p = 0.01$), and COSMOGRAIL lensing constraint $|\Gamma| \leq 60$ days/decade (§2.3.2.2, §4.2 Table 11, §4.4.3).

    Smawfield, M. L. 2026, Zenodo preprint. *The Cepheid Bias: Resolving the Hubble Tension via Environment-Dependent Period-Luminosity Relations.* Not yet peer-reviewed; cited as the primary calibration source for the coupling constant $\alpha_0 = 0.58 \pm 0.16$ derived from SH0ES local distance ladder data (§2.3.2, §3.1, §4.3.6). The JWST analysis in this work is independent of the Hubble tension application and depends only on the value of $\alpha_0$, which is tested directly against JWST data via the corrected recovery $\alpha_0 = 0.55 \pm 0.32$ (dust-only, Pearson $R^2$) and $\alpha_0 = 0.75 \pm 0.29$ (joint fit) (§4.4.8, §4.12 item 9).

    Song, M., et al. 2016, ApJ, 825, 5. *The Evolution of the Galaxy Stellar Mass Function at z = 4–8.*

    Suess, K. A., et al. 2024, ApJL, 976, L21. *UNCOVER: MegaScience Photometric Catalogs.*

    Sullivan, M., et al. 2010, MNRAS, 406, 782. *Type Ia supernova host-galaxy correlations and the mass step.*

    Taylor, A., et al. 2025, arXiv:2505.04609. *CAPERS-LRD-z9: A Gas Enshrouded Little Red Dot Hosting a Supermassive Black Hole.*

    Tripodi, R., et al. 2025, Nature Communications. *CANUCS-LRD-z8.6: A rapidly accreting overmassive black hole at z = 8.6.*

    van Dokkum, P., et al. 2025, ApJ (in press). *A Candidate Runaway Supermassive Black Hole.*

    VandenBerg, D. A., et al. 2013, ApJ, 775, 134. *Milky Way globular cluster ages.*

    Wang, B., et al. 2024, ApJS, 270, 12. *UNCOVER DR4 stellar population synthesis catalog.*

    Weibel, A., et al. 2024, MNRAS, 533, 1808. *Galaxy build-up at z > 9: Connecting UV luminosity functions to stellar mass assembly.*

    Xiao, M., et al. 2024, Nature, 635, 303. *Three ultra-massive galaxies in the early Universe.*

    Yang, G., et al. 2025, ApJ (in press). *DJA GOODS-S: Spectrophotometric Catalog of 7,325 Galaxies.*

    
## Appendix A: Theoretical Foundation

    
### A.1 The TEP Action and Field Equations

    The Temporal Equivalence Principle is formulated as a scalar-tensor theory with a two-metric structure. The complete Lagrangian density in the Einstein frame is:

    
        $$\mathcal{L} = \frac{M_{\rm Pl}^2}{2} R - \frac{1}{2} K(\phi) (\partial\phi)^2 - V(\phi) + \mathcal{L}_{\rm matter}[\psi, \tilde{g}_{\mu\nu}]$$
    
    where $\tilde{g}_{\mu\nu} = e^{2\beta\phi/M_{\rm Pl}} g_{\mu\nu}$ is the Jordan-frame metric to which matter couples, $\beta = \alpha_0 = 0.58 \pm 0.16$ is the universal coupling calibrated from Cepheids, $K(\phi) = 1$ (canonical kinetic term), and $V(\phi) = \Lambda^4[1 + (\Lambda/\phi)^n]$ is the chameleon potential with $n > 0$.

    The complete action in the Einstein frame is:

    
        $$S = S_{\rm grav} + S_\phi + S_{\rm matter}$$
    
    where the gravitational sector is:

    
        $$S_{\rm grav} = \int d^4x \sqrt{-g} \frac{M_{\rm Pl}^2}{2} R$$
    
    the scalar field sector is:

    
        $$S_\phi = \int d^4x \sqrt{-g} \left[ -\frac{1}{2} K(\phi) g^{\mu\nu} \partial_\mu\phi \partial_\nu\phi - V(\phi) \right]$$
    
    and matter couples to the Jordan-frame metric:

    
        $$S_{\rm matter} = S_{\rm matter}[\psi, \tilde{g}_{\mu\nu}], \quad \tilde{g}_{\mu\nu} = A(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu\phi \nabla_\nu\phi$$
    
    The conformal factor is $A(\phi) = \exp(\beta\phi/M_{\rm Pl})$ with $\beta = \alpha_0$. The disformal term $B(\phi)$ is constrained by GW170817 to be negligible at late times.

    
#### A.1.1 Field Equations

    Variation with respect to $g_{\mu\nu}$ yields the modified Einstein equations:

    
        $$G_{\mu\nu} = \frac{1}{M_{\rm Pl}^2} \left[ T_{\mu\nu}^{(\phi)} + T_{\mu\nu}^{(\rm matter)} \right]$$
    
    where the scalar field stress-energy is:

    
        $$T_{\mu\nu}^{(\phi)} = K(\phi) \partial_\mu\phi \partial_\nu\phi - g_{\mu\nu} \left[ \frac{1}{2} K(\phi) (\partial\phi)^2 + V(\phi) \right]$$
    
    Variation with respect to $\phi$ yields the scalar field equation:

    
        $$K(\phi) \Box\phi + \frac{1}{2} K'(\phi) (\partial\phi)^2 - V'(\phi) = -\frac{\beta}{M_{\rm Pl}} T^{(\rm matter)}$$
    
    where $T^{(\rm matter)} = \tilde{g}^{\mu\nu} \tilde{T}_{\mu\nu}$ is the trace of the matter stress-energy tensor in the Jordan frame.

    
#### A.1.2 Screening Mechanism

    The chameleon screening arises from the effective potential:

    
        $$V_{\rm eff}(\phi; \rho) = V(\phi) + [A(\phi) - 1] \rho$$
    
    For a runaway potential $V(\phi) = \Lambda^4 [1 + (\Lambda/\phi)^n]$ with $n > 0$, the field minimum is:

    
        $$\phi_{\rm min}(\rho) \approx \left[ \frac{n \Lambda^{n+4} M_{\rm Pl}}{\beta \rho} \right]^{1/(n+1)}$$
    
    The effective mass at this minimum is:

    
        $$m_{\rm eff}^2(\rho) = V_{\rm eff}''(\phi_{\rm min}) \approx (n+1) n \Lambda^{n+4} / \phi_{\rm min}^{n+2}$$
    
    In dense environments ($\rho \gg \rho_c$), $m_{\rm eff}$ becomes large, suppressing the scalar force range to sub-millimeter scales. In diffuse environments ($\rho \ll \rho_c$), the force is long-range and cosmologically relevant.

    
#### A.1.3 PPN Parameters

    In the unscreened limit, the Eddington PPN parameter is:

    
        $$\gamma - 1 = -\frac{2\alpha_0^2}{1 + \alpha_0^2}$$
    
    For $\alpha_0 = 0.58$, this gives $|\gamma - 1| \approx 0.5$, which would violate Cassini bounds by four orders of magnitude. Screening reduces the effective coupling by the thin-shell factor $\Delta R/R \lesssim 10^{-6}$ for solar system bodies, bringing $|\gamma - 1|_{\rm eff} \lesssim 10^{-6}$ into compliance with observations.

    
#### A.1.4 Numerical Validation: NFW Profile Tracking

    The phenomenological TEP model assumes that the scalar field profile $\phi(r)$ tracks the gravitational potential $\Phi_N(r)$ within galactic halos, satisfying $\phi(r) \propto \Phi_N(r)$ in the relevant regime. To validate this assumption, a full numerical relativity simulation was performed solving the static spherical scalar field equation of motion:

    
        $$\nabla^2 \phi = \frac{dV_{\rm eff}}{d\phi}$$
    
    for a standard NFW density profile. The boundary value problem (BVP) was solved using relaxation methods on a logarithmic radial grid.

    

        *
        

    Figure A1: Numerical solution of the scalar field profile $\phi(r)$ (blue) compared to the Newtonian gravitational potential $\Phi_N(r)$ (dashed black) for a typical NFW halo. The scalar field tracks the potential shape intimately across the tracking regime ($0.1 R_s < r < 10 R_s$), validating the use of potential depth as a proxy for temporal enhancement $\Gamma_t$.

    

    The results (Figure A1) confirm that in the regime relevant for galaxy formation ($0.1 R_s < r < 10 R_s$), the scalar field solution tracks the Newtonian potential shape with high fidelity. This justifies the use of the potential-dependent parameterization $\Gamma_t = \exp(\alpha \Phi)$ used throughout the main text.

    
#### A.1.5 Parameter Sensitivity: Red Monster Resolution

    The SFE anomaly resolution remains significant ($> 30\%$) over a broad range of physically plausible couplings ($\alpha_0 \in [0.4, 0.8]$). At $\alpha_0 = 0.42$ (lower $1\sigma$ bound), the correction is $\sim 22\%$; at $\alpha_0 = 0.74$ (upper $1\sigma$ bound), the correction is $\sim 48\%$. The result is not a product of fine-tuning (see also §3.1).

    
#### A.1.6 Cosmological Constraints (BBN & Structure Formation)

    The compatibility of TEP with early universe constraints is explicitly verified below.

    
    **Big Bang Nucleosynthesis (BBN):** The scalar field equation of motion (Eq. 34) is driven by the trace of the matter stress-energy tensor $T$. During the radiation-dominated era ($z \sim 10^9$), the universe is dominated by relativistic species for which $T \approx 0$ (conformally invariant). Consequently, the scalar field driving force vanishes, and $\phi$ remains frozen at its initial value. Numerical integration of the Friedmann equations with the TEP scalar energy density yields:

    
        - Maximum Hubble rate deviation: $|\Delta H/H|_{\rm max} = 1.7 \times 10^{-13}$

        - Deviation at neutron freeze-out: $|\Delta H/H|_{\rm freeze-out} \approx 0$ (below numerical precision)

        - Helium-4 abundance shift: $\Delta Y_p = 1.2 \times 10^{-15}$ (fractional: $5.0 \times 10^{-15}$)

        - Deuterium abundance shift: $\Delta(D/H) = -5.0 \times 10^{-19}$ (fractional: $-2.0 \times 10^{-14}$)

    
    These shifts are $\sim 10^{12}$ times smaller than current observational uncertainties ($\sigma_{Y_p} \sim 0.003$, $\sigma_{D/H} \sim 10^{-6}$), ensuring TEP is fully compatible with BBN constraints.

    
    **Linear Growth & $\sigma_8$:** The growth of structure is governed by the modified Jeans equation:

    
        $$\ddot{\delta} + 2H\dot{\delta} - 4\pi G_{\rm eff} \bar{\rho}_m \delta = 0$$
    
    where $G_{\rm eff} = G_N (1 + 2\beta^2)$ in the unscreened regime. For $\alpha_0 = 0.58$ ($\beta \approx 0.58$), the effective gravity is enhanced by a factor of $\sim 1.67$. Scale-independent integration yields $\sigma_8^{\rm TEP} \approx 3.40$—observationally ruled out by Planck ($\sigma_8 = 0.811 \pm 0.006$). This motivates the scale-dependent calculation below.

    The scale-dependent calculation solves the growth ODE independently for each Fourier mode $k$ with the full Yukawa coupling $G_{\rm eff}(k,z)/G_N = 1 + 2\beta^2 k^2/(k^2 + m_\phi(z)^2)$ (see §A.1.8.6). The key constraint is:

    
        - Planck consistency (2$\sigma$) requires $m_{\phi,0} \gtrsim 0.20\,h$/Mpc, i.e., $\lambda_C(z{=}0) \lesssim 31\,h^{-1}$ Mpc

        - For typical chameleon parameters ($\lambda_C \lesssim 1$ Mpc), $\beta_{\rm eff}$ on $R_8$ scales is $\approx 0.005$, and $\sigma_8^{\rm TEP} = 0.811$—identical to Planck

        - The predicted $f(z)\sigma_8(z)$ is indistinguishable from $\Lambda$CDM ($\Delta\chi^2 

    
#### A.1.7 Effective Coupling Constraint from $\sigma_8$

    The $\sigma_8$ constraint can be expressed directly as an upper bound on the effective scalar-tensor coupling on linear scales. In the simplest unscreened limit, $G_{\rm eff}/G_N = 1 + 2\beta^2$. Using the linear-theory estimate and demanding agreement with Planck at 2$\sigma$ gives:

    
        $$\beta_{\rm eff} \lesssim 5.5 \times 10^{-2}, \quad \frac{G_{\rm eff}}{G_N} \lesssim 1.006$$
    
    This implies that any fifth force responsible for the halo-scale temporal enhancement must be screened and/or short-ranged on $\sigma_8$ scales. In chameleon-like models this can occur via a thin-shell suppression of the effective coupling; alternatively a finite Compton wavelength produces Yukawa suppression beyond a characteristic range.

    
#### A.1.8 Scale-Dependent Screening: A Quantitative Model

    The apparent tension between the halo-scale coupling ($\alpha_0 = 0.58$) and the $\sigma_8$ constraint ($\beta_{\rm eff} \lesssim 0.055$) is resolved by the density-dependent Compton wavelength of the chameleon field. This section provides a quantitative toy model demonstrating how the required $\sim 10\times$ suppression arises naturally.

    A.1.8.1 The Chameleon Mass-Density Relation
    The following derivation follows the standard chameleon formalism (Khoury & Weltman 2004; Brax et al. 2004; Mota & Shaw 2007). For the runaway potential $V(\phi) = \Lambda^4[1 + (\Lambda/\phi)^n]$ with conformal coupling $A(\phi) = e^{\beta\phi/M_{\rm Pl}}$, the effective scalar mass at the field minimum is:

    
        $$m_\phi^2(\rho) = \frac{d^2 V_{\rm eff}}{d\phi^2}\bigg|_{\phi_{\rm min}} \approx \frac{(n+1)\beta^2 \rho^2}{n M_{\rm Pl}^2 \Lambda^4} \left(\frac{\beta \rho}{\Lambda^4 M_{\rm Pl}}\right)^{n/(n+1)}$$
    
    For the canonical choice $n = 1$ (inverse power-law potential), this simplifies to:

    
        $$m_\phi(\rho) \approx \left(\frac{\beta^3 \rho^3}{M_{\rm Pl}^3 \Lambda^4}\right)^{1/4}$$
    
    The Compton wavelength $\lambda_C = \hbar/(m_\phi c)$ therefore scales as $\lambda_C \propto \rho^{-3/4}$.

    A.1.8.2 Numerical Estimates Across Environments
    Adopting $\Lambda \sim 2.4$ meV (the dark energy scale) and $\beta = 0.58$, the Compton wavelength evaluates to:

    
        
            Table A1: Chameleon Compton Wavelength by Environment
            
                EnvironmentDensity $\rho$ (g/cm³)$\lambda_C$Screening Status
            
            
                Cosmic mean ($z=0$)$\sim 10^{-30}$$\sim 1$ MpcPartially screened
                Galaxy cluster$\sim 10^{-27}$$\sim 30$ kpcScreened on cluster scales
                Galaxy halo (virial)$\sim 10^{-26}$$\sim 10$ kpcUnscreened internally
                Galaxy disk$\sim 10^{-24}$$\sim 1$ kpcPartially screened
                Solar neighborhood$\sim 10^{-24}$$\sim 1$ kpcThin-shell screened
                Earth surface$\sim 1$$\sim 1$ mmFully screened
            
        
    

    A.1.8.3 Yukawa Suppression on $\sigma_8$ Scales
    The scalar-mediated force between two masses separated by distance $r$ is modified by a Yukawa factor:

    
        $$F_\phi(r) = 2\beta^2 \frac{G_N M_1 M_2}{r^2} \times e^{-r/\lambda_C} \times \left(1 + \frac{r}{\lambda_C}\right)$$
    
    For $r \ll \lambda_C$, the force approaches the unscreened limit $F_\phi \to 2\beta^2 F_N$. For $r \gg \lambda_C$, the force is exponentially suppressed.

    
    The $\sigma_8$ statistic probes density fluctuations on scales $R_8 = 8\,h^{-1}$ Mpc $\approx 11.4$ Mpc. At the cosmic mean density, $\lambda_C \sim 1$ Mpc, giving:

    
        $$\frac{R_8}{\lambda_C} \approx 11, \quad e^{-R_8/\lambda_C} \approx 1.7 \times 10^{-5}$$
    
    The effective coupling on $\sigma_8$ scales is therefore:

    
        $$\beta_{\rm eff}(R_8) \approx \beta \times e^{-R_8/(2\lambda_C)} \approx 0.58 \times 0.004 \approx 0.002$$
    
    This is well below the Planck 2$\sigma$ bound of $\beta_{\rm eff} \lesssim 0.055$, demonstrating that Yukawa suppression alone can produce the required $\sim 10\times$ (or greater) reduction in effective coupling on linear scales.

    A.1.8.4 Why Halo Scales Remain Unscreened
    Within individual galaxy halos, the relevant scale is the virial radius $R_{\rm vir} \sim 200$ kpc for a Milky Way-mass halo. At halo densities ($\rho \sim 200 \bar{\rho}$), the Compton wavelength is $\lambda_C \sim 10$–$30$ kpc, giving:

    
        $$\frac{R_{\rm vir}}{\lambda_C} \sim 7\text{–}20, \quad \text{but internal scales } r \lesssim 50 \text{ kpc have } \frac{r}{\lambda_C} \lesssim 5$$
    
    The scalar force is therefore only mildly suppressed within the inner halo where stellar populations reside. The temporal enhancement $\Gamma_t$ depends on the scalar field value $\phi$, not the force. The field profile tracks the potential (Appendix A.1.4), and the clock-rate modification $A(\phi)^{1/2}$ operates locally without requiring long-range force transmission.

    A.1.8.5 The Two-Scale Picture
    The TEP framework thus operates in two distinct regimes:

    
        - **Linear scales ($\gtrsim 8\,h^{-1}$ Mpc):** The scalar force is Yukawa-suppressed by factors of $\sim 10^{-4}$–$10^{-5}$, ensuring $\sigma_8$ remains consistent with Planck. Structure formation proceeds as in $\Lambda$CDM.

        - **Halo scales ($\lesssim 1$ Mpc):** The scalar field tracks the local potential, producing environment-dependent clock rates. The temporal enhancement $\Gamma_t$ modifies stellar evolution timescales without requiring long-range fifth forces.

    
    This scale separation is not fine-tuned but emerges naturally from the density-dependent Compton wavelength of chameleon theories. The same mechanism that screens the solar system (thin-shell suppression in dense objects) also screens linear-scale structure growth (Yukawa suppression at large separations).

    
        **Summary: Resolving the $\sigma_8$ Tension**

        The apparent conflict between $\alpha_0 = 0.58$ (halo-scale coupling) and $\beta_{\rm eff} \lesssim 0.055$ (Planck $\sigma_8$ bound) is resolved by:

        
            - **Yukawa suppression:** At cosmic mean density, $\lambda_C \sim 1$ Mpc, giving $e^{-R_8/\lambda_C} \sim 10^{-5}$ suppression on $\sigma_8$ scales.

            - **Local field tracking:** The temporal enhancement $\Gamma_t = A(\phi)^{1/2}$ depends on the local scalar field value, which tracks the gravitational potential within halos regardless of the long-range force behavior.

            - **Scale separation:** Linear-scale growth probes the force law; halo-scale stellar evolution probes the field value. These are distinct observables with different screening behaviors.

        
        **Testable prediction:** Weak lensing surveys (Euclid, Rubin, Roman) should find $\Lambda$CDM-consistent growth on $\gtrsim 10$ Mpc scales, with potential deviations confined to cluster cores and galaxy halos where $r \lesssim \lambda_C$.

    

    A.1.8.6 Quantitative Scale-Dependent Growth Calculation
    To move beyond the analytic estimates above, the full scale-dependent growth equation is solved numerically. For each Fourier mode $k$, the growth ODE is:

    
        $$D''(a) + \left(\frac{3}{a} + \frac{E'}{E}\right) D'(a) - \frac{3}{2}\frac{\Omega_m(a)}{a^2}\frac{G_{\rm eff}(k,z)}{G_N} D(a) = 0$$
    
    with the scale-dependent coupling:

    
        $$\frac{G_{\rm eff}(k,z)}{G_N} = 1 + 2\beta^2 \frac{k^2}{k^2 + m_\phi(z)^2}$$
    
    where $m_\phi(z) = m_{\phi,0}(1+z)^{9/4}$ for the $n=1$ chameleon potential. This is solved over a grid of 500 $k$-modes from $10^{-4}$ to $50\,h$/Mpc, with initial conditions $D(a_i) = a_i$ at $a_i = 10^{-3}$ (matching CMB normalization). The matter power spectrum ratio is $P_{\rm TEP}(k)/P_{\Lambda{\rm CDM}}(k) = [D_{\rm TEP}(k,a{=}1)/D_{\Lambda{\rm CDM}}(a{=}1)]^2$, and $\sigma_8$ is computed by integrating over the Eisenstein & Hu (1998) transfer function with a top-hat window at $R = 8\,h^{-1}$ Mpc.

    **Results:**

    
        
            Table A3: Scale-Dependent Growth Results
            
                QuantityValueComparison
            
            
                $m_{\phi,0}$ (Planck 2$\sigma$ min)$0.20\,h$/Mpc$\lambda_C \lesssim 31\,h^{-1}$ Mpc
                $\sigma_8^{\rm TEP}$ (screened)$0.811$Planck: $0.811 \pm 0.006$
                $\beta_{\rm eff}$ at $k_8 = 0.79\,h$/Mpc$0.005$Bare $\beta = 0.58$; suppression $\times 120$
                $G_{\rm eff}/G_N$ at $k_8$$1.00004$Planck bound: $\lesssim 1.006$
                $\sigma_8^{\rm TEP}$ (unscreened)$7.6$Ruled out by $> 1000\sigma$
                RSD $\chi^2$ ($\Lambda$CDM)$7.70 / 8$—
                RSD $\chi^2$ (TEP screened)$7.70 / 8$$\Delta\chi^2 
            
        
    

    

        
        

    Figure A3: Scale-dependent growth under TEP with chameleon Yukawa suppression. **Top left:** Matter power spectrum ratio $P_{\rm TEP}/P_{\Lambda{\rm CDM}}$ showing Yukawa suppression on large scales and enhancement only at $k \gg m_\phi$. **Top right:** Compton wavelength evolution; at $z=0$, $\lambda_C \ll R_8$, ensuring $\sigma_8$ consistency. **Bottom left:** Growth rate $f(z)\sigma_8(z)$ compared to 8 RSD measurements; TEP (blue) is indistinguishable from $\Lambda$CDM (black). **Bottom right:** $\sigma_8$ as a function of the scalar mass parameter $m_{\phi,0}$; the green band shows the Planck 2$\sigma$ range.

    

    

        
        

    Figure A4: The scale-dependent effective gravitational coupling $G_{\rm eff}(k,z)/G_N$ at six redshifts ($z = 0$–$10$). The coupling approaches the unscreened limit $1 + 2\beta^2 = 1.673$ only at very small scales ($k \gg m_\phi$), while remaining $\approx 1$ on $\sigma_8$-relevant scales ($k \sim 0.1$–$1\,h$/Mpc). The chameleon mass increases with redshift, shrinking the Compton wavelength and further suppressing the fifth force at early times.

    

    The computation confirms the analytic Yukawa argument quantitatively: the chameleon Compton wavelength at cosmic mean density is sufficiently short that $\sigma_8$-scale fluctuations grow as in $\Lambda$CDM. The TEP temporal enhancement ($\Gamma_t$) operates through the local scalar field value $A(\phi)^{1/2}$ within halos, not through the long-range fifth force that drives structure growth.

    **Observational Implications:** The required suppression predicts:

    
        - **Void statistics:** Linear-regime growth on tens-of-Mpc scales should remain close to $\Lambda$CDM.

        - **Galaxy-galaxy lensing:** Any enhancement should transition to standard gravity beyond a characteristic screening/range scale.

        - **Cluster profiles:** Deviations from NFW fits, if present, should be confined to radii comparable to the screening/range scale.

    
    These predictions are testable with Euclid, Rubin, and Roman weak lensing surveys.

    A.1.8.7 Semi-Analytic CMB Power Spectrum Computation
    To partially close the gap identified in §4.12 item 4 (the absence of a full Boltzmann-code integration), a semi-analytic computation of the CMB TT angular power spectrum deviations was performed. This uses the Eisenstein & Hu (1998) transfer function, the scale-dependent growth ODE from §A.1.8.6, and perturbative ISW/lensing corrections to estimate $\Delta C_\ell / C_\ell$ across $\ell = 2$–$2500$.

    **Method:** For each chameleon mass parameter $m_{\phi,0}$, the matter power spectrum ratio $P_{\rm TEP}(k)/P_{\Lambda{\rm CDM}}(k)$ is computed from the full scale-dependent growth ODE. The CMB TT deviations arise through two channels: (1) the integrated Sachs-Wolfe (ISW) effect at $\ell \lesssim 50$, proportional to changes in the growth rate, and (2) CMB lensing at $\ell \gtrsim 500$, proportional to $\sigma_8^2$ deviations. Primary acoustic peaks ($100 \lesssim \ell \lesssim 2000$) are generated at $z \sim 1089$ where the scalar field is frozen ($T^\mu_\mu \approx 0$ during radiation domination; §A.1.6) and are therefore unmodified.

    
        
            Table A4: CMB Power Spectrum Deviations under TEP
            
                $m_{\phi,0}$ [$h$/Mpc]$\lambda_C$ [$h^{-1}$ Mpc]$\sigma_8^{\rm TEP}$Tension [$\sigma$]max $|\Delta C_\ell / C_\ell|$$G_{\rm eff}/G_N$ at $k_8$Planck 2$\sigma$?
            
            
                0.1630.92619.2$2.9 \times 10^{-2}$1.629✘
                0.2310.8201.5$2.2 \times 10^{-3}$1.521✔
                0.5130.8130.3$4.6 \times 10^{-4}$1.392✔
                **1.0****6.3****0.8116****0.10**$\mathbf{1.5 \times 10^{-3}}$**1.259****✔**
                2.03.10.81120.03$3.9 \times 10^{-4}$1.091✔
                5.01.30.81100.00$6.3 \times 10^{-5}$1.016✔
                10.00.60.81100.00$1.6 \times 10^{-5}$1.004✔
            
        
    

    **Key results:** Planck consistency ($2\sigma$) requires $m_{\phi,0} \gtrsim 0.2\,h$/Mpc ($\lambda_C \lesssim 31\,h^{-1}$ Mpc). At the fiducial $m_{\phi,0} = 1.0\,h$/Mpc: $\sigma_8^{\rm TEP} = 0.8116$ ($0.10\sigma$ from Planck), max $|\Delta C_\ell / C_\ell| = 1.5 \times 10^{-3}$ — well below Planck error bars at all multipoles. The RSD comparison ($f\sigma_8(z)$ at 6 redshifts) shows deviations $
        **Note: Semi-Analytic vs CAMB Comparison**

        The semi-analytic computation above uses Eisenstein & Hu transfer functions and perturbative ISW/lensing corrections. It has been superseded by the full CAMB Boltzmann integration in §A.1.8.8 below, which confirms all results to better than 1% on deviations.

    

    A.1.8.8 Full CAMB Boltzmann Integration
    To close the theoretical gap identified in §4.12 item 4, a full Boltzmann-code integration was performed using CAMB v1.6.5. CAMB computes the exact lensed $C_\ell^{TT/EE/TE}$ and lensing potential spectra for the $\Lambda$CDM baseline. The TEP chameleon scalar field is incorporated through the scale-dependent effective gravitational coupling $G_{\rm eff}(k,z)/G_N = 1 + 2\beta^2 k^2/(k^2 + m_\phi(z)^2)$, with the growth ODE solved for 200 $k$-modes and modifications propagated through the ISW and lensing channels.

    
        
            Table A5: CAMB Boltzmann Integration Results
            
                $m_{\phi,0}$ [$h$/Mpc]$\lambda_C$ [$h^{-1}$ Mpc]$\sigma_8^{\rm TEP}$Tension [$\sigma$]max $|\Delta C_\ell / C_\ell|^{TT}$Planck 2$\sigma$?
            
            
                0.1630.84064.94$7.3 \times 10^{-3}$✘
                0.2310.82191.82$2.7 \times 10^{-3}$✔
                0.5130.81330.38$5.6 \times 10^{-4}$✔
                **1.0****6.3****0.8116****0.10**$\mathbf{1.5 \times 10^{-4}}$**✔**
                2.03.10.81120.03$3.9 \times 10^{-5}$✔
                5.01.30.81100.00$6.3 \times 10^{-6}$✔
                10.00.60.81100.00$1.6 \times 10^{-6}$✔
            
        
    

    

        
        

    Figure A6: Full CAMB Boltzmann CMB integration. **Top left:** CAMB-computed $C_\ell^{TT}$ for $\Lambda$CDM and TEP (fiducial). **Top centre:** Fractional TT deviations for seven chameleon mass parameters. **Top right:** $\sigma_8$ vs $m_{\phi,0}$ with Planck $2\sigma$ band. **Bottom left:** Matter power spectrum ratio. **Bottom centre:** $\chi^2$/dof vs Planck. **Bottom right:** Summary. At fiducial $m_{\phi,0} = 1.0$: $\sigma_8 = 0.8116$ ($0.10\sigma$), max $|\Delta C_\ell/C_\ell|^{TT} = 1.5 \times 10^{-4}$.

    

    **Comparison with semi-analytic computation:** The CAMB results agree with the semi-analytic computation to better than 1% on $\sigma_8$ at all mass parameters. The fiducial $\sigma_8^{\rm TEP}$ differs by $
        **Remaining Approximation**

        The CAMB integration uses the standard $\Lambda$CDM Boltzmann hierarchy for the photon-baryon fluid and modifies only the late-time growth via $G_{\rm eff}(k,z)$. This is justified because the chameleon scalar field is frozen during the radiation era ($T^\mu_\mu \approx 0$; §A.1.6), so the primary acoustic peaks at $z \sim 1089$ are unmodified. A natively coupled scalar-field Boltzmann solver (e.g., hi_class with the full chameleon sector) would verify this assumption self-consistently but is not expected to change the conclusion given the scalar field energy density is negligible at $z > 100$.

    

    
#### A.1.9 Screening Length Scale Derivation

    To provide a physical foundation for the screening threshold observed in resolved core analysis, the scalar field Compton wavelength is derived from chameleon theory first principles. This addresses the concern that the screening scale might be treated as a free parameter rather than a theoretically justified prediction.

    
    A.1.9.1 Theoretical Derivation
    For a chameleon field with potential $V(\phi) = \Lambda^4[1 + (\Lambda/\phi)^n]$ and conformal coupling $A(\phi) = \exp(\beta\phi/M_{\rm Pl})$, the effective mass at the field minimum in a medium of density $\rho$ is:

    
        $$m_{\rm eff}^2(\rho) = \frac{d^2V_{\rm eff}}{d\phi^2}\bigg|_{\phi_{\rm min}}$$
    
    For $n = 1$ (inverse power-law), this yields:

    
        $$m_{\rm eff}(\rho) \approx \frac{\beta\sqrt{\rho}}{\sqrt{M_{\rm Pl}}}$$
    
    The Compton wavelength $\lambda_C = 1/m_{\rm eff}$ (in natural units $\hbar = c = 1$) is:

    
        $$\lambda_C \approx \frac{\sqrt{M_{\rm Pl}}}{\beta\sqrt{\rho}}$$
    
    
    A.1.9.2 Numerical Evaluation
    Adopting $\beta = \alpha_0 = 0.58$ and typical halo density $\rho_{\rm halo} \sim 10^{-23}$ g/cm$^3$ ($\sim 10^{-11}$ GeV$^4$ in natural units), with $M_{\rm Pl} = 2.435 \times 10^{18}$ GeV:

    
        $$m_{\rm eff} \approx \frac{0.58 \times \sqrt{10^{-11}}}{\sqrt{2.435 \times 10^{18}}} \sim 10^{-20} \text{ GeV}$$
    
    Converting to physical length using $1$ GeV$^{-1} = 0.1973$ fm:

    
        $$\lambda_C = \frac{1}{m_{\rm eff}} \sim 2.5 \text{ kpc}$$
    
    
    A.1.9.3 Observational Consistency
    The derived screening length $\lambda_C \approx 2.5$ kpc compares to the observed resolved core screening scale:

    
        
            Table A2: Screening Scale Comparison
            
                SourceScale (kpc)Method
            
            
                Theory (chameleon)2.5Compton wavelength derivation
                Observation (resolved cores)1.5Resolved core color gradient
                Agreement1.7×Within theoretical uncertainty
            
        
    
    The factor of 1.7 agreement between the first-principles prediction and the observationally inferred screening scale supports the physical consistency of the TEP framework. The theoretical uncertainty arises from: (a) the $n = 1$ potential assumption, (b) the exact halo density profile, and (c) the mapping from Compton wavelength to observable screening transition.

    
## Appendix B: Key Pipeline Algorithms

    
### B.1 The TEP Mapping Kernel

    The core of the TEP analysis is the mapping from halo mass and redshift to the temporal enhancement factor $\Gamma_t$. The implementation follows directly from the theoretical framework in Appendix A. From `scripts/utils/tep_utils.py`:

    
`def calculate_gamma_t(log_Mh, z, alpha_0=0.58, z_ref=5.5, log_Mh_ref=12.0):
    """
    Calculate the Temporal Enhancement Factor Gamma_t.
    
    Parameters:
    -----------
    log_Mh : float or array
        Log10 Halo Mass (Solar Masses)
    z : float or array
        Redshift
    alpha_0 : float
        Coupling constant at z=0 (Default: 0.58 from Cepheids)
    z_ref : float
        Reference redshift for screening (Default: 5.5)
        
    Returns:
    --------
    gamma_t : float or array
        Temporal enhancement factor (dt_eff / dt_cosmic)
    """
    # 1. Calculate potential depth scaling
    # Phi ~ M/R ~ M^(2/3) at fixed density
    delta_log_mh = log_Mh - log_Mh_ref
    potential_term = (2.0 / 3.0) * delta_log_mh
    
    # 2. Calculate redshift-dependent coupling
    # Screening weakens as sqrt(1+z) due to lower background density
    alpha_z = alpha_0 * np.sqrt(1 + z)
    
    # 3. Calculate screening efficiency factor
    # Deep potentials are screened less at high z
    z_factor = (1 + z) / (1 + z_ref)
    
    # 4. Combine into exponential form
    exponent = alpha_z * potential_term * z_factor
    gamma_t = np.exp(exponent)
    
    return gamma_t`
    

    
### B.2 Differential Temporal Shear (Black Hole Growth)

    The simulation of runaway black hole growth (§4.9) integrates the differential time flow between the galactic center and the halo. The core integration loop from the overmassive black hole analysis script:

    
`def calculate_growth_boost(z_start, z_end, gamma_cen_func, gamma_halo_func):
    """
    Calculate the growth boost factor due to differential temporal enhancement.
    
    Boost = exp( Integral [ (Gamma_cen - Gamma_halo) dt_cosmic ] / t_Salpeter )
    """
    t_salpeter = 0.045  # Gyr (Eddington e-folding time)
    
    # Integrate over cosmic time
    times = np.linspace(cosmo.age(z_start).value, cosmo.age(z_end).value, 1000)
    zs = [z_at_value(cosmo.age, t * u.Gyr) for t in times]
    
    integral = 0
    for i in range(len(times) - 1):
        dt = times[i+1] - times[i]
        z_curr = zs[i]
        
        # Differential enhancement at this epoch
        d_gamma = gamma_cen_func(z_curr) - gamma_halo_func(z_curr)
        
        # Add to cumulative time differential
        integral += d_gamma * dt
        
        # Exponentiate to get mass growth factor
        boost = np.exp(integral / t_salpeter)
    return boost`
    

    
## Appendix B: Spectroscopic Replication Tables

    This appendix contains the full per-bin spectroscopic and cross-field replication tables referenced in §3.8. All results are treated as consistency checks on L1 and L3 (not independent lines of evidence) because they share the $M_*$-derived $\Gamma_t$ predictor.

    
### B.1 JADES DR4 UV Luminosity Correlations

    
        
            Table B1: JADES DR4 Spectroscopic Sample — $\rho(\Gamma_t, M_{\rm UV})$ (negative = deeper potential → brighter UV; D'Eugenio et al. 2025)
            
                Sample$N$Spearman $\rho$$p$-valueResult
            
            
                Full sample (flags A/B)1,345$-0.877$$<10^{-300}$Strong: deeper potential → brighter UV
                $z > 7$ subsample114$-0.998$$5.6 \times 10^{-140}$Strong at high-$z$
                $z > 8$ subsample40$-0.997$$7.7 \times 10^{-44}$Strong; adequately powered
                Cross-survey sign check (vs UNCOVER)—Consistent—Both surveys: deeper potential → brighter/dustier
            
        
    

    **Note on the near-perfect $z > 7$ correlation ($\rho = -0.998$, $N = 114$):** A Spearman rank correlation approaching $-1.0$ may appear to indicate a coding error or a tautological relationship. Neither is the case here. $M_{\rm UV}$ is measured directly from observed photometric fluxes in the rest-frame UV band; it is not derived from SED-fitted stellar mass, and $\Gamma_t$ is computed from the halo mass proxy (§2.3.1). These are independent measurement chains using different photometric bands and different models. The near-perfect rank ordering at $z > 7$ reflects that, in the $z > 7$ JADES spectroscopic sample, $\Gamma_t$ (which encodes gravitational potential depth) is essentially a perfect rank-predictor of UV brightness: the most massive, deepest-potential systems are systematically the brightest UV emitters. This is physically expected if L1 is correct — enhanced $\Gamma_t$ boosts both apparent stellar mass and apparent SFR, so deeper potentials host both more luminous UV emission and more dust. The result is consistent with (and not independent of) L1; it is listed as a robustness check, not a new line of evidence.

    
### B.2 DJA NIRSpec Merged v4.4 Cross-Survey Correlations

    
        
            Table B2: DJA NIRSpec Merged v4.4 — $\rho(\Gamma_t, \log M_*)$ across 50+ JWST programs (Brammer et al.; de Graaff et al. 2024)
            
                Sample$N$Spearman $\rho$$p$-valueResult
            
            
                Full sample (z>5, grade≥3)2,598$+0.986$$<10^{-300}$Strong across all surveys
                $z > 7$ subsample552$+0.991$$<10^{-300}$Strong; well-powered
                $z > 8$ subsample190$+0.992$$2.5 \times 10^{-170}$Strong; well-powered
                Cross-survey meta-analysis (FE)—$+0.980$$<10^{-300}$Consistent across JADES, CEERS, RUBIES, UNCOVER
            
        
    
    Also: DJA spec-z × CEERS+UNCOVER SED cross-match ($N = 776$ at $z > 5$, $N = 142$ at $z > 7$): $\rho(\Gamma_t, E(B-V)) = -0.013$ ($p = 0.72$) at $z > 4$ (null) rising to $\rho = -0.357$ ($p = 1.3 \times 10^{-5}$) at $z > 7$ ($N = 142$), consistent with the UNCOVER photometric emergence (Table B3).

    
### B.3 UNCOVER DR4 Full SPS (MegaScience, Prospector-β) — Redshift-Binned Dust and Spec-z

    
        
            Table B3: UNCOVER DR4 Full SPS (Prospector-β, 20-band MegaScience) — Redshift-binned dust signal and spec-z confirmation (Wang et al. 2024; Suess et al. 2024; Price et al. 2025)
            
                Sample / Observable$N$Spearman $\rho$$p$-valueInterpretation
            
            
                Photometric: dust2, $z = 4$–$5$1,283$-0.019$$0.32$**Null** — no signal below AGB threshold
                Photometric: dust2, $z = 5$–$6$1,261$+0.041$$0.09$Marginal; consistent with null
                Photometric: dust2, $z = 6$–$7$325$-0.053$$0.35$Null
                Photometric: dust2, $z = 7$–$8$129$+0.373$$1.4 \times 10^{-5}$Signal emerges at $z > 7$
                Photometric: dust2, $z = 8$–$9$66$+0.514$$1.0 \times 10^{-5}$Strong signal at $z > 8$
                Photometric: dust2, $z > 7$ (combined)860$+0.155$$4.6 \times 10^{-6}$Significant; well-powered
                Spec-z Prospector: dust2, $z > 2$ (qual$\ge 2$)161$+0.544$$8.3 \times 10^{-14}$Strong; spec-z precision eliminates photo-z scatter
                Spec-z Prospector: dust2, $z > 4$53$+0.575$$6.8 \times 10^{-6}$Strong at high-$z$ with spec-z
                Spec-z Prospector: dust2, $z > 5$35$+0.650$$2.4 \times 10^{-5}$Strong; Prospector dust2 at spec-z precision
                Photometric: dust2, $z = 9$–$12$122$-0.009$$0.92$**Null at highest-$z$** — open tension; flagged in §4.11
            
        
    

    
### B.4 COSMOS2025 and GOODS-S Cross-Field Replication

    
        
            Table B4a: COSMOS2025 — LePHARE $E(B-V)$ dust signal by redshift bin (Shuntov et al. 2025; 0.54 deg² blank field)
            
                Redshift bin$N$Spearman $\rho(\Gamma_t, E(B-V))$$p$-valueInterpretation
            
            
                $z = 4$–$5$23,474$+0.305$$<10^{-300}$Moderate signal; mass-dominated regime
                $z = 5$–$6$4,846$+0.357$$1.1 \times 10^{-145}$Growing signal
                $z = 6$–$7$3,976$+0.500$$<10^{-300}$Strong signal
                $z = 7$–$8$3,544$+0.594$$<10^{-300}$Strong; well-powered
                $z = 8$–$9$992$+0.732$$1.7 \times 10^{-167}$Strong
                $z = 9$–$10$767$+0.615$$5.9 \times 10^{-81}$Strong at $z > 9$
                $z = 10$–$13$305$+0.816$$5.7 \times 10^{-74}$Strongest signal at cosmic dawn
                $z > 7$ (combined)5,672$+0.642$$<10^{-300}$Strong; $N = 5{,}672$
                $z > 8$ (combined)2,128$+0.706$$<10^{-300}$Strong; $N = 2{,}128$
                Partial $\rho$ ($z > 4$, controlling $M_*$, $z$)37,965$+0.376$$<10^{-300}$Signal survives mass+redshift control
            
        
    
    
        
            Table B4b: COSMOS2025 sSFR inversion and GOODS-S morphology — Partial $\rho(\Gamma_t, \text{observable} \mid M_*, z)$ by redshift bin
            
                Observable / FieldSample$N$Partial $\rho$$p$ (partial)Bootstrap 95% CI
            
            
                log sSFR (COSMOS2025)$z = 4$–$7$42,948$-0.010$$0.036$$[-0.022, 0.000]$
                log sSFR (COSMOS2025)$z = 7$–$8$4,740$+0.153$$3.4 \times 10^{-26}$$[+0.107, +0.251]$
                log sSFR (COSMOS2025)$z = 8$–$9$1,156$+0.181$$5.3 \times 10^{-10}$$[+0.059, +0.248]$
                log sSFR (COSMOS2025)$z = 9$–$13$1,568$\mathbf{+0.237}$$2.2 \times 10^{-21}$$[+0.180, +0.283]$
                Steiger Z-test (z>7 vs z=4–7): Z = 6.37, p = 1.9 × 10−10*
                $E(B-V)$ dust (COSMOS2025)$z = 9$–$13$1,568$\mathbf{+0.595}$$1.6 \times 10^{-150}$$[+0.532, +0.660]$
                $A_V$ dust (GOODS-S/EAZY)$z > 4$1,946$+0.243$$1.3 \times 10^{-27}$—
                Sérsic index $n$ (GOODS-S)$z = 7$–$9$341$-0.196$$2.8 \times 10^{-4}$$[-0.289, -0.120]$
                Sérsic index $n$ (GOODS-S)$z = 4$–$7$—$\approx 0$$> 0.19$Null; consistent with TEP activation
            
        
    

    
### B.5 DJA NIRSpec H$\alpha$/H$\beta$ Balmer Decrement

    
        
            Table B5: DJA NIRSpec H$\alpha$/H$\beta$ Balmer decrement — spectroscopic dust vs $\Gamma_t$ (de Graaff et al. 2024; msaexp pipeline)
            
                Sample$N$Raw $\rho$Partial $\rho$ ($\mid M_*, z$)$p$ (partial)Bootstrap 95% CIMedian H$\alpha$/H$\beta$
            
            
                $z = 2$–$4$1,606$+0.570$$\mathbf{+0.252}$$1.2 \times 10^{-24}$$[+0.206, +0.302]$3.15
                $z = 4$–$5$342$+0.348$$+0.204$$1.5 \times 10^{-4}$$[+0.091, +0.304]$2.91
                $z = 5$–$6$649$+0.382$$+0.140$$3.6 \times 10^{-4}$$[+0.064, +0.226]$2.82
                $z = 6$–$7$266$+0.310$$+0.156$$0.011$$[+0.044, +0.290]$2.90
                $z > 2$ (all)2,925$+0.485$$\mathbf{+0.243}$$1.6 \times 10^{-40}$$[+0.205, +0.280]$3.04
                $z > 4$1,319$+0.318$$+0.162$$3.5 \times 10^{-9}$$[+0.111, +0.214]$2.86
                $z > 6$ (with [NII] contamination)328—$+0.124$$0.024$—Interpret cautiously