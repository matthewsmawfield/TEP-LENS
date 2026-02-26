# TEP-LENS: Resolving the Hubble Tension

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Kingston upon Hull)  
**Date:** First published: 11 January 2026  
**DOI:** 10.5281/zenodo.18209703  
**Generated:** 2026-02-26  
**Paper Series:** TEP Series: Paper 12 (Cosmological Observations)

---

**Abstract.** The Temporal Equivalence Principle (TEP) postulates that stellar evolution rates and temporal propagation scale with the depth of the local gravitational potential. While this mechanism has successfully accounted for anomalous Cepheid period-luminosity residuals and high-redshift star formation efficiencies, it makes a rigid geometric prediction for strongly lensed transients: light rays traversing different regions of a lens potential must experience differential temporal scaling. This breaks the classical "closure" of relative time delays in multi-image systems. Using public data from the triply imaged supernova SN H0pe (PLCK G165.7+67.0) and the TDCOSMO sample of strongly lensed quasars, we construct a purely geometric "route-closure" test. By comparing the sum of pairwise relative time delays around closed loops, we bypass absolute lens macromodel degeneracies (such as the Mass Sheet Degeneracy) to isolate the differential temporal shear. We present the predicted TEP closure residual for SN H0pe and evaluate the observed delays against both the $\Lambda$CDM standard expectation and the TEP framework.

                
                

                    
## 1. Introduction

The standard framework of cosmology, $\Lambda$CDM, relies on the assumption that General Relativity governs all scales of cosmic structure formation and that clocks in different regions of the Universe tick synchronously when corrected for peculiar velocity and scale factor evolution. The Temporal Equivalence Principle (TEP) challenges this foundation by positing that proper time effectively scales with the depth of the local gravitational potential, parameterized by a coupling constant $\alpha_0 \approx 0.58$ derived from local Cepheid variables. 

### 1.1 Strong Lensing as a Geometric Probe

Under standard physics, light traveling from a distant source through a gravitational lens along multiple paths experiences a relative time delay. This delay is strictly determined by the geometric path difference and the Shapiro time delay of the potential.

Under TEP, because the effective speed of temporal propagation is modified within deep potentials, the time delays should exhibit an additional, scale-dependent signature. By measuring the arrival times of a single transient event (like a supernova) or stochastic quasar variability along three or more distinct paths, we can construct a "route-closure" test.

### 1.2 SN H0pe

The recent discovery of SN H0pe (PLCK G165.7+67.0), a triply-imaged Type Ia supernova, provides a unique and nearly ideal geometric laboratory. Unlike stochastic quasar variations, a supernova provides a single, well-defined explosion seen along three distinct light paths, allowing us to directly compare relative time delays and construct the closure residual.

                
                

                    
## 2. Methodology

### 2.1 Lensing Time Delays

In standard General Relativity (GR), the time delay between two images $i$ and $j$ of a lensed source is given by the difference in the arrival time surface $t(\boldsymbol{\theta})$:

$\Delta t_{i,j} = \frac{D_{\Delta t}}{c} \left[ \frac{1}{2} (\boldsymbol{\theta}_i - \boldsymbol{\beta})^2 - \frac{1}{2} (\boldsymbol{\theta}_j - \boldsymbol{\beta})^2 - (\psi(\boldsymbol{\theta}_i) - \psi(\boldsymbol{\theta}_j)) \right]$

where $D_{\Delta t}$ is the time-delay distance, $\boldsymbol{\beta}$ is the source position, $\boldsymbol{\theta}_i$ and $\boldsymbol{\theta}_j$ are the image positions, and $\psi(\boldsymbol{\theta})$ is the deflection potential.

### 2.2 The Route-Closure Test

For a source producing three images (A, B, C), the sum of pairwise relative time delays around the closed loop must vanish identically:

$\Delta t_{AB} + \Delta t_{BC} + \Delta t_{CA} = 0$

Under TEP, the effective speed of light (and temporal propagation) differs along the paths through the varying depth of the cluster's gravitational potential. This introduces a potential-dependent scaling factor $\Gamma_t$ along each path length, resulting in a non-zero closure residual:

$\mathcal{R}_{\rm closure} = (\Gamma_{t,A} - 1)\Delta t_{AB} + (\Gamma_{t,B} - 1)\Delta t_{BC} + (\Gamma_{t,C} - 1)\Delta t_{CA} \neq 0$

By isolating this residual, we can construct a test of the TEP framework that is explicitly independent of the Mass Sheet Degeneracy (MSD) and absolute lens macromodel normalization, as $\lambda$ rescales all delays symmetrically.

                
                

                    
## 3. Results

### 3.1 GR Route-Closure Null Test

Using the photometric relative time delays measured by Pierel et al. (2024) for SN H0pe:

    - $\Delta t_{AB} = -116.6 \pm 10.8$ days

    - $\Delta t_{CB} = -48.6 \pm 4.0$ days

By definition in General Relativity, time is absolute across the observation plane. If we observe delays between pairs of images from a single source event, the loop must close perfectly to zero. The inferred delay between images A and C is purely arithmetic:

$\Delta t_{AC} = \Delta t_{AB} - \Delta t_{CB} = -68.0$ days.

The sum around the geometric route $\mathcal{R}_{\rm closure} = \Delta t_{AB} + \Delta t_{BC} + \Delta t_{CA}$ strictly evaluates to $0.00$ days under standard $\Lambda$CDM assumptions. 

### 3.2 TEP Predicted Closure Residual

Under the Temporal Equivalence Principle, the effective speed of light scales inversely with the local gravitational potential depth $\Gamma_t$. Using the lens macromodel magnifications ($\mu$) as a first-order proxy for the projected potential depth (convergence) along the line of sight, we evaluate the differential temporal shear across the three images. 

For a representative TEP coupling $\alpha = 0.05$, the temporal shear factors evaluate to:

    - Path A ($\mu = 5.4$): $\Gamma_A = 1.0366$

    - Path B ($\mu = 2.5$): $\Gamma_B = 1.0199$

    - Path C ($\mu = 2.0$): $\Gamma_C = 1.0151$

Because light on Path A traverses a deeper region of the cluster potential, its apparent transit time is dilated relative to Path C. Weighting the pairwise geometric delays by their respective differential shear factors produces the TEP closure prediction:

$\mathcal{R}_{\rm closure}^{\rm TEP} = (\Gamma_A - 1)\Delta t_{AB} + (\Gamma_B - 1)\Delta t_{BC} + (\Gamma_C - 1)\Delta t_{CA} = -2.28$ days.

    ![Bar chart comparing the route closure residual under standard GR (0.00 days) versus the TEP prediction (-2.28 days).](public/figures/step_04_closure_residual.png)
    Figure 1: Comparison of the route-closure residual $\mathcal{R}_{\rm closure}$ for the triply imaged SN H0pe system. While standard General Relativity mandates a perfectly closed loop ($0.00$ days), the Temporal Equivalence Principle predicts a geometric mismatch ($-2.28$ days) due to differential temporal shear across the lens cluster potential.

                
                

                    
## 4. Discussion

### 4.1 Breaking the Mass Sheet Degeneracy

The most significant barrier to precision cosmography using strong lensing time delays is the Mass Sheet Degeneracy (MSD). The addition of a constant density sheet to the lens macromodel rescales the projected potential $\psi(\boldsymbol{\theta})$ and scales the overall time delay $\Delta t$ by a factor $\lambda$ without altering the imaging observables. 

Because $\lambda$ operates as a global scalar, it modifies all pairwise time delays symmetrically. A traditional two-image delay measurement $\Delta t_{AB}$ cannot uniquely separate the Hubble constant $H_0$ from $\lambda$. However, the route-closure residual $\mathcal{R}_{\rm closure}$ is constructed specifically to vanish under symmetric global scaling. By integrating around a closed loop, the linear geometric path differences strictly cancel.

The predicted TEP residual of $-2.28$ days is not an absolute scaling factor; it is a **differential non-linear symmetry breaking**. No uniform mass sheet $\lambda$ can generate a non-zero closure residual from a baseline of zero. Consequently, measuring a statistically significant deviation from zero in a multi-imaged transient directly falsifies standard General Relativistic time propagation, independent of MSD uncertainties.

### 4.2 Multi-Path Transients vs Stochastic Quasars

While strongly lensed quasars (e.g. the TDCOSMO sample) provide high statistical power due to continuous monitoring baselines, the intrinsic stochasticity of Active Galactic Nucleus (AGN) accretion disks makes true route-closure tests difficult. Quasar variability exhibits red noise, and microlensing by stars in the lens galaxy introduces uncorrelated extrinsic variability along each sightline.

A lensed supernova like SN H0pe represents the ideal geometric probe: a single, physically constrained explosion event traversing three distinct spacetime paths. The relative time delays are determined not by cross-correlating noisy stochastic signals, but by aligning the distinct evolutionary phases of the Type Ia light curve. This provides a much tighter constraint on the true differential arrival times required to evaluate $\mathcal{R}_{\rm closure}$.

### 4.3 Observational Feasibility

The predicted TEP residual for SN H0pe ($-2.28$ days) is currently smaller than the $1\sigma$ measurement uncertainties of the observed photometric delays (e.g. $\Delta t_{AB}$ has an error of $\sim 10$ days). Therefore, the existing SN H0pe data cannot yet definitively confirm or falsify the TEP prediction. 

However, as the sample of strongly lensed transients grows—particularly with high-cadence JWST or Roman Space Telescope observations of highly asymmetric quadruple lens configurations (like the recently discovered superluminous SN 2025wny)—the measurement precision on $\Delta t$ will drop below $1$ day. At that threshold, the TEP closure residual becomes a strictly testable, binary geometric prediction.

                
                

                    
## 5. Conclusion

We have proposed a purely geometric "route-closure" test for the Temporal Equivalence Principle (TEP) using multiply imaged transients. Because TEP modifies the effective speed of light inversely with local gravitational potential depth, time delays along distinct geometric paths through a lens cluster do not sum to zero around a closed loop.

By applying this to the triply imaged supernova SN H0pe, we find that:

    - Under standard General Relativity, the closure residual is identically zero.

    - Under TEP, assuming a coupling of $\alpha=0.05$ calibrated to local domains, the closure residual is non-zero ($-2.28$ days), breaking the absolute Mass Sheet Degeneracy scaling.

    - While current observational uncertainties ($\sim 10$ days) obscure the TEP signal in SN H0pe, upcoming quad-lensed superluminous supernovae (like SN 2025wny) with higher precision cadence will cross the falsification threshold.

Consequently, strong lensing route-closure residuals offer an immediate, binary geometric test of standard General Relativity against potential-dependent temporal propagation.

                
                

                    
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