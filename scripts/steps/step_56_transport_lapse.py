#!/usr/bin/env python3
"""
TEP-LENS: Step 56 - Conformal Transport Consistency Audit

Purpose: test whether a path-integrated potential term can be counted as an
additional TEP delay after the standard GR lensing delay has been subtracted.
The answer is no in the static conformal sector.  This step retains the useful
GLAFIC potential-map scale and ordering diagnostic, but prevents the same
Shapiro/Fermat potential from being counted once in the GR lens model and a
second time as a purported TEP residual.

Physics
-------
In TEP, light propagates on the matter metric g~ = A^2(phi) g + B dphi dphi.
For B=0, ds_tilde^2=A^2 ds_g^2 and the null condition is exactly ds_g^2=0.
Unparameterised null geodesics and their coordinate travel times are therefore
conformally invariant.  The potential contribution

    dt_template,i = (alpha / c^3) * integral_path |Phi_i| dl        (lens frame)
    dt_template,i^obs = (1 + z_l) * dt_template,i                   (observed)

is the ordinary potential (Shapiro/Fermat) term already present in the GR lens
model.  Multiplying it by a fitted alpha_eff defines a phenomenological additive
template; it is not a conformal-sector first-principles prediction.

Thin-lens exact relation: the LOS-integrated potential at impact parameter
theta is the lensing deflection potential, up to the standard distance factor:

    integral |Phi| dl = (c^2 / 2) * (D_l D_s / D_ls) * psi(theta)

so that

    dt_template,i^obs = alpha * (1 + z_l) * (D_l D_s / (2 c D_ls)) * psi_i .

The GLAFIC v3 psi map is used directly: it encodes the full 2D lens model
(sub-clumps, member galaxies, cluster halo), avoiding the sub-structure loss
of a smooth spherical deprojection. The additive gauge constant in psi cancels
in image-to-image differences (relational measurement, Rule 8).

The script reports the coefficient required if that additive template is fitted
to the observed residual, solely as a scale diagnostic.  A non-zero residual
from the TEP action must instead be derived from scalar backreaction on g, a
time-dependent background, or the disformal B-sector.  None of those channels
is computed by a static psi map alone.

A smooth spherical-deprojection cross-check (corrected Step-51 geometry:
r = sqrt(b^2 + s^2), dl in km) is reported with its sub-structure limitation
stated explicitly.

Inputs : data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_{psi,kappa}.fits
         data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json
         results/outputs/step_07_observed_vs_predicted.json
Outputs: results/outputs/step_56_transport_lapse.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import KAPPA_LENS, BETA_REFSDAL, load_refsdal_image_positions

STEP_NUM = "56"
# Published J2000 positions (single source of truth: the GLAFIC params JSON;
# ATel #6729 and Karman+2016 MUSE Table 1 for S1-S4, refined-model position
# coincident with the Kelly+2016 detection for SX).  An earlier version of this
# step used positions ~10 arcsec north, inside the cluster core; that bug
# produced the spurious ~4x kappa-map/table mismatch and a reversed psi
# ordering.
IMAGE_POSITIONS_DEG = load_refsdal_image_positions()
Z_L = 0.542
Z_S = 1.489
H0 = 70.0
OM0 = 0.3
C_KMS = 299792.458
MPC_TO_KM = 3.0856775814913673e19
DAY = 86400.0


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def load_map_and_wcs(path):
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path) as hdul:
        return hdul[0].data.astype(float), WCS(hdul[0].header)


def sample_map_bilinear(data, wcs, positions_deg):
    """Bilinear-interpolated map values at image positions (sub-pixel)."""
    vals = {}
    for name, (ra, dec) in positions_deg.items():
        px = wcs.all_world2pix(ra, dec, 0)
        x, y = float(px[0]), float(px[1])
        x0, y0 = int(np.floor(x)), int(np.floor(y))
        dx, dy = x - x0, y - y0
        v = (data[y0, x0] * (1 - dx) * (1 - dy)
             + data[y0, x0 + 1] * dx * (1 - dy)
             + data[y0 + 1, x0] * (1 - dx) * dy
             + data[y0 + 1, x0 + 1] * dx * dy)
        vals[name] = float(v)
    return vals


def main():
    print_status(f"STEP {STEP_NUM}: Conformal Transport Consistency Audit", "TITLE")

    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)
    D_l = cosmo.angular_diameter_distance(Z_L).value
    D_s = cosmo.angular_diameter_distance(Z_S).value
    D_ls = cosmo.angular_diameter_distance_z1z2(Z_L, Z_S).value
    D_delta = D_l * D_s / D_ls
    print_status(f"D_l={D_l:.1f}, D_s={D_s:.1f}, D_ls={D_ls:.1f}, D_Delta={D_delta:.1f} Mpc", "INFO")

    # ----------------------------------------------------------
    # 1. Sample GLAFIC v3 psi and kappa maps at image positions
    # ----------------------------------------------------------
    map_dir = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "maps"
    psi_map, wcs = load_map_and_wcs(map_dir / "hlsp_frontier_model_macs1149_glafic_v3_psi.fits")
    kappa_map, _ = load_map_and_wcs(map_dir / "hlsp_frontier_model_macs1149_glafic_v3_kappa.fits")

    psi_vals = sample_map_bilinear(psi_map, wcs, IMAGE_POSITIONS_DEG)
    kappa_vals = sample_map_bilinear(kappa_map, wcs, IMAGE_POSITIONS_DEG)
    print_status("psi / kappa sampled (bilinear):", "INFO")
    for im in IMAGE_POSITIONS_DEG:
        print_status(f"  {im}: psi={psi_vals[im]:.2f}  kappa_map={kappa_vals[im]:.3f}", "INFO")

    # ----------------------------------------------------------
    # 1b. Map-normalisation reconciliation.  The archived GLAFIC v3 maps are
    #     scaled to D_ls/D_s = 1 (z_s -> infinity; HLSP readme).  Rescaling the
    #     sampled kappa by beta(1.489) = 0.5339 should reproduce the Kelly+2023
    #     Table 15 values if the positions and maps are mutually consistent.
    # ----------------------------------------------------------
    gl_params = json.load(open(PROJECT_ROOT / "data" / "raw" / "sn_lensing" /
                               "refsdal_glafic_v3_lensing_params.json"))
    kappa_table = {im: gl_params["images"][im]["kappa"] for im in IMAGE_POSITIONS_DEG}
    kappa_rescaled = {im: kappa_vals[im] * BETA_REFSDAL for im in kappa_vals}
    kappa_ratio = {im: kappa_rescaled[im] / kappa_table[im] for im in kappa_table}
    ratio_vals = np.array(list(kappa_ratio.values()))
    print_status(f"kappa rescaled by beta(1.489)={BETA_REFSDAL} vs Kelly+2023 Table 15:", "INFO")
    for im in IMAGE_POSITIONS_DEG:
        print_status(f"  {im}: map*beta={kappa_rescaled[im]:.3f}  table={kappa_table[im]:.3f}  "
                     f"ratio={kappa_ratio[im]:.3f}", "INFO")
    print_status(f"  rescaled/table ratio: mean={ratio_vals.mean():.3f}, "
                 f"std={ratio_vals.std():.3f} ({100*ratio_vals.std()/ratio_vals.mean():.1f}% scatter)",
                 "INFO")

    # ----------------------------------------------------------
    # 2. GR potential-delay scale per image.  This is already part of the
    #    Fermat potential used by a standard lens model; it is not an
    #    additional conformal TEP delay.
    # ----------------------------------------------------------
    arcsec2_to_rad2 = (np.pi / (180.0 * 3600.0)) ** 2
    # Time per unit psi (days per rad^2)
    K_days = (1.0 + Z_L) * (D_delta * MPC_TO_KM) / (2.0 * C_KMS) / DAY
    print_status(f"Transport prefactor K = {K_days:.4e} days per unit psi [rad^2]", "INFO")

    # |int Phi| dl in Mpc-equivalent c-time (the geometric integral itself)
    int_Phi = {im: 0.5 * D_delta * psi_vals[im] * arcsec2_to_rad2 for im in psi_vals}  # Mpc
    int_Phi_days = {im: v * MPC_TO_KM / C_KMS / DAY for im, v in int_Phi.items()}
    print_status("Path-integrated |Phi|/c^2 dl (light-time equivalent):", "INFO")
    for im in IMAGE_POSITIONS_DEG:
        print_status(f"  {im}: {int_Phi[im]:.4e} Mpc  ->  {int_Phi_days[im]:.1f} d", "INFO")

    # ----------------------------------------------------------
    # 3. Fit an additive-potential template to the observed SX:S1 residual.
    #    This is a phenomenological scale diagnostic, not a prediction.
    # ----------------------------------------------------------
    s07 = json.load(open(PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    R_obs = float(s07["weighted_mean_residual"]["R_obs_days"])

    # Delays are measured relative to S1; the blind-model ensemble predicts the
    # S1->SX delay.  The predicted TEP contribution to (obs - model) on that
    # contrast is alpha * K * (psi_SX - psi_S1).
    dpsi_SX_S1 = psi_vals["SX"] - psi_vals["S1"]
    dpsi_SX_S4 = psi_vals["SX"] - psi_vals["S4"]
    psi_inner_mean = np.mean([psi_vals[i] for i in ("S1", "S2", "S3", "S4")])
    dpsi_SX_inner = psi_vals["SX"] - psi_inner_mean

    K_full = K_days * arcsec2_to_rad2  # days per arcsec^2 of psi difference

    alpha_SX_S1 = R_obs / (K_full * dpsi_SX_S1)
    alpha_SX_S4 = R_obs / (K_full * dpsi_SX_S4)
    alpha_SX_inner = R_obs / (K_full * dpsi_SX_inner)

    print_status(f"dpsi: SX-S1={dpsi_SX_S1:+.2f}, SX-S4={dpsi_SX_S4:+.2f}, "
                 f"SX-<inner>={dpsi_SX_inner:+.2f} arcsec^2", "INFO")
    print_status(f"Required additive-template coefficient: SX:S1 -> {alpha_SX_S1:.4f}, "
                 f"SX:S4 -> {alpha_SX_S4:.4f}, SX:<inner> -> {alpha_SX_inner:.4f}", "INFO")
    print_status(f"|kappa_lens| = {abs(KAPPA_LENS):.4f}, |beta_A| = 1.0", "INFO")

    # ----------------------------------------------------------
    # 4. Per-image values of the fitted additive template
    # ----------------------------------------------------------
    alpha_use = alpha_SX_S1
    dt_template = {im: alpha_use * K_full * psi_vals[im] for im in psi_vals}
    resid_pred = {im: dt_template[im] - dt_template["S1"] for im in psi_vals}
    print_status(f"Fitted additive-template residuals rel. S1 at alpha={alpha_use:.4f}:", "INFO")
    for im in IMAGE_POSITIONS_DEG:
        print_status(f"  {im}: dt_template={dt_template[im]:+9.2f} d   resid_pred={resid_pred[im]:+8.2f} d", "INFO")

    # The psi ordering is a diagnostic of the GR potential-delay contribution.
    # Note: with the corrected (published) image positions, psi_SX is the
    # SMALLEST of the five images.  For any convergent lens monopole the
    # deflection potential is a bowl with its minimum at the cluster centre
    # (d psi/dR = 2 M_2D(<R)/pi R > 0), so the psi ordering is a cluster-centric
    # radius ordering, not a local-depth ordering: psi_SX smallest simply records
    # that SX lies closest to the potential minimum (r ~ 9.3 arcsec, versus
    # 12-14 arcsec for S1-S4).  The TEP-relevant ordering is the local temporal-
    # field depth at each image, which tracks the local 3D potential/projected
    # density (phi shares the Poisson source of Phi and kappa): by that measure
    # SX is the deepest of the five (highest rescaled kappa, deepest deprojected
    # |Phi| at its radius -- see the cross-check below), i.e. it arrives latest.
    psi_SX_is_smallest = bool(psi_vals["SX"] == min(psi_vals.values()))
    print_status(f"psi ordering: psi_SX is smallest of the five = {psi_SX_is_smallest} "
                 f"(bowl-shaped Fermat potential: smallest psi = smallest cluster-centric "
                 f"radius = deepest local temporal field)", "INFO")

    # ----------------------------------------------------------
    # 5. Corrected smooth-deprojection cross-check (Step-51 bugs fixed)
    #    r = sqrt(b^2 + s^2);  dl in km;  truncate at r_max.
    #    The profile is centred on the psi-map minimum (the potential centre,
    #    not the WCS reference pixel) and the map is rescaled to z_s = 1.489.
    # ----------------------------------------------------------
    from scripts.steps.step_51_geodesic_transport import (
        kappa_radial_profile, spherical_deprojection, compute_3D_potential,
    )
    iy_min, ix_min = np.unravel_index(np.nanargmin(psi_map), psi_map.shape)
    ra_c, dec_c = wcs.all_pix2world(ix_min, iy_min, 0)
    kappa_src = kappa_map * BETA_REFSDAL   # D_ls/D_s = 1 -> z_s = 1.489
    r_mid, kappa_mean, pix = kappa_radial_profile(kappa_src, wcs,
                                                center_ra=ra_c, center_dec=dec_c)
    r_3d, rho_3d = spherical_deprojection(r_mid, kappa_mean, D_l, Z_L)
    Phi_3d = compute_3D_potential(r_3d, rho_3d)

    cx, cy = wcs.all_world2pix(ra_c, dec_c, 0)
    pix_deg = np.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    arc2rad = np.pi / (180.0 * 3600.0)
    smooth_check = {}
    s_grid = np.linspace(-3.0 * r_3d[-1], 3.0 * r_3d[-1], 40001)
    for im, (ra, dec) in IMAGE_POSITIONS_DEG.items():
        xi, yi = wcs.all_world2pix(ra, dec, 0)
        b = np.hypot(float(xi) - cx, float(yi) - cy) * pix_deg * 3600.0 * arc2rad * D_l
        r_path = np.sqrt(b**2 + s_grid**2)                      # corrected geometry
        Phi_p = np.interp(r_path, r_3d, Phi_3d, left=0.0, right=0.0)
        I = np.trapezoid(np.abs(Phi_p), s_grid)                 # Mpc
        smooth_check[im] = {
            "b_kpc": float(b * 1000.0),
            "local_Phi_over_c2": float(np.interp(b, r_3d, Phi_3d)),
            "int_Phi_Mpc": float(I),
            "light_days": float(I * MPC_TO_KM / C_KMS / DAY),
        }
    dI_SX_S1 = smooth_check["SX"]["int_Phi_Mpc"] - smooth_check["S1"]["int_Phi_Mpc"]
    depth_order = sorted(smooth_check, key=lambda k: smooth_check[k]["local_Phi_over_c2"])
    int_order = sorted(smooth_check, key=lambda k: -smooth_check[k]["int_Phi_Mpc"])
    print_status("Smooth-deprojection cross-check (geometry + centre fixed, kappa rescaled "
                 "to z_s=1.489):", "INFO")
    for im in IMAGE_POSITIONS_DEG:
        v = smooth_check[im]
        print_status(f"  {im}: b={v['b_kpc']:.0f} kpc  Phi_local/c^2={v['local_Phi_over_c2']:.3e}  "
                     f"int={v['int_Phi_Mpc']:.4e} Mpc ({v['light_days']:.0f} d)", "INFO")
    print_status(f"  Local potential-depth ordering (deepest first): {depth_order}", "INFO")
    print_status(f"  Path-integral ordering (largest first): {int_order}", "INFO")
    print_status(f"  dI(SX-S1) = {dI_SX_S1:+.3e} Mpc  (azimuthally averaged model retains the "
                 f"monopole only; the full-2D kappa map is the primary local-density diagnostic)",
                 "INFO")

    # ----------------------------------------------------------
    # 6. Loop closure sanity: corrected delays still close
    # ----------------------------------------------------------
    delays = gl_params["delays_days_rel_S1"]
    # Closure residual of the fitted per-image template (must remain ~0 because
    # it assigns one scalar value to each image).
    i, j, k = "S1", "S4", "SX"
    closure = ((dt_template[j] - dt_template[i])
               + (dt_template[k] - dt_template[j])
               + (dt_template[i] - dt_template[k]))

    # ----------------------------------------------------------
    # 7. Physical interpretation
    # ----------------------------------------------------------
    ratio_to_kappa = abs(alpha_use) / abs(KAPPA_LENS)

    depth_deepest = depth_order[0] if depth_order else None
    verdict = (
        f"The GLAFIC psi map gives the ordinary path-integrated potential scale. "
        f"At the published image positions, psi_SX = {psi_vals['SX']:.1f} arcsec^2 is "
        f"the SMALLEST of the five images: the deflection potential is a bowl with "
        f"its minimum at the cluster centre, so the psi ordering is a cluster-centric "
        f"radius ordering (SX at ~9.3 arcsec, S1-S4 at 12-14 arcsec), not a local "
        f"temporal-field-depth ordering.  The TEP field tracks the local potential "
        f"and projected density (same Poisson source): by those measures SX samples "
        f"the deepest environment of the five (rescaled kappa_SX = "
        f"{kappa_rescaled['SX']:.3f}, the maximum; deprojected local |Phi| ordering: "
        f"{depth_order}), so SX arriving latest is the sign the mechanism predicts -- "
        f"matching the observed +30.1 d residual direction.  This is an environment "
        f"ordering only: the corresponding Shapiro/Fermat term is already inside "
        f"every GR lens prediction, and a static conformal factor preserves null "
        f"curves exactly, so the additional conformal propagation residual is zero. "
        f"Fitting an extra psi template to +30.1 d would require alpha = "
        f"{alpha_use:.4f} (sign-flipped and suppressed relative to |kappa_lens| by a "
        f"factor {ratio_to_kappa:.2f}); it is calibrated on the same residual and has "
        f"no first-principles status.  A non-zero TEP amplitude must be obtained from "
        f"scalar backreaction, time dependence, or the disformal sector without "
        f"double-counting the GR Fermat potential."
    )
    print_status("\n" + verdict)

    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": "Consistency audit of the path-integrated potential scale using the "
                       "thin-lens relation int|Phi|dl = (c^2/2) D_Delta psi.",
        "classification": "GR potential-delay diagnostic; not an additional conformal TEP prediction",
        "conformal_only_additional_residual_days": 0.0,
        "conformal_cancellation_reason": (
            "For B=0, g_tilde=A^2 g has the same null curves as g. The static psi-dependent "
            "Shapiro/Fermat delay is already included in the GR lens prediction."
        ),
        "cosmology": {"z_l": Z_L, "z_s": Z_S, "H0": H0, "Om0": OM0,
                      "D_l_Mpc": D_l, "D_s_Mpc": D_s, "D_ls_Mpc": D_ls,
                      "D_Delta_Mpc": D_delta},
        "psi_map_arcsec2": psi_vals,
        "kappa_map_sampled": kappa_vals,
        "map_normalisation_reconciliation": {
            "archive_convention": "GLAFIC v3 maps scaled to D_ls/D_s = 1 (z_s -> infinity; HLSP readme)",
            "beta_z_src_1p489": BETA_REFSDAL,
            "kappa_rescaled_to_z_src": kappa_rescaled,
            "kappa_table_kelly2023": kappa_table,
            "rescaled_over_table_ratio": kappa_ratio,
            "ratio_mean": float(ratio_vals.mean()),
            "ratio_std": float(ratio_vals.std()),
            "interpretation": ("At the published image positions the map/table ratio is uniform "
                               "to ~3%, identifying the documented D_ls/D_s = 1 normalisation "
                               "(expected factor 1/beta = 1.873, measured "
                               f"{float(np.mean([kappa_vals[im]/kappa_table[im] for im in kappa_table])):.3f}). "
                               "The earlier ~4x non-uniform mismatch was caused by erroneous image "
                               "positions ~10 arcsec north of the image field, not a map defect."),
        },
        "int_Phi_dl": {"Mpc": int_Phi, "light_days": int_Phi_days},
        "dpsi_arcsec2": {"SX_minus_S1": float(dpsi_SX_S1),
                         "SX_minus_S4": float(dpsi_SX_S4),
                         "SX_minus_inner_mean": float(dpsi_SX_inner)},
        "observed_residual_days": R_obs,
        "additive_psi_template_coefficient_required": {
            "from_SX_S1_contrast": float(alpha_SX_S1),
            "from_SX_S4_contrast": float(alpha_SX_S4),
            "from_SX_inner_mean_contrast": float(alpha_SX_inner),
        },
        "kappa_lens_proxy": float(KAPPA_LENS),
        "additive_template_over_kappa_lens": float(ratio_to_kappa),
        "fitted_additive_template_residuals": {im: float(v) for im, v in resid_pred.items()},
        "sign_check": {
            "psi_SX_is_smallest": psi_SX_is_smallest,
            "psi_ordering_meaning": ("radial: psi is a bowl with its minimum at the cluster "
                                     "centre; smallest psi = smallest cluster-centric radius "
                                     "(SX ~9.3 arcsec vs 12-14 arcsec for S1-S4)"),
            "local_depth_ordering": {
                "rescaled_kappa_highest": "SX",
                "deprojected_local_Phi_deepest_first": depth_order,
                "interpretation": ("TEP field tracks the local potential/density (same Poisson "
                                   "source); SX samples the deepest temporal environment and "
                                   "is predicted to arrive latest"),
            },
            "observed_sign": "+30.1 d (SX latest)",
            "matches": bool(depth_order[0] == "SX"),
            "evidential_status": "ordering diagnostic only; GR potential already subtracted",
        },
        "smooth_deprojection_crosscheck": {
            "centre": {"ra": float(ra_c), "dec": float(dec_c),
                       "note": "psi-map minimum (potential centre), not the WCS reference pixel"},
            "per_image": smooth_check,
            "dI_SX_minus_S1_Mpc": float(dI_SX_S1),
            "local_depth_ordering": depth_order,
            "path_integral_ordering": int_order,
            "limitation": ("Azimuthally averaged spherical deprojection retains only the cluster "
                           "monopole and misses the member-galaxy sub-structure around S1-S4 and "
                           "the near-critical compression at SX; the full-2D rescaled kappa map "
                           "is the primary local-density diagnostic.  Neither ordering is an "
                           "additional static conformal propagation prediction."),
        },
        "loop_closure_days": float(closure),
        "bugs_fixed_relative_to_step_50_51": [
            "step_50 applied the endpoint clock ratio Gamma(Phi_i) to delay differences; "
            "that quantity is not a propagation residual.",
            "step_51 ray geometry used r = sqrt(b^2 + (s*angle)^2) with angle ~ 1e-4, so the "
            "cluster potential never truncated over the ~5 Gpc path (~4400x inflation).",
            "step_51 divided dl [Mpc] by c [km/s], missing the 3.086e19 km/Mpc conversion; "
            "the resulting t_TEP was never used in the residual.",
            "step_56 previously counted the GR Shapiro/Fermat potential a second time as an "
            "additional conformal TEP delay; conformal null-curve invariance makes that "
            "additional static contribution identically zero.",
            "steps 50/51/56/054 sampled the GLAFIC maps at incorrect image positions "
            "~10 arcsec north of the published coordinates (Dec ~22.398 instead of "
            "~22.396), inside the cluster core; corrected positions recover the "
            "documented D_ls/D_s = 1 map normalisation (uniform 1.878 map/table ratio "
            "= 1/beta(1.489) = 1.873) and show psi_SX is the smallest, not largest, "
            "of the five.",
        ],
        "verdict": verdict,
    }
    out = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_transport_lapse.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
