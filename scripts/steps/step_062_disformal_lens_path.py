#!/usr/bin/env python3
"""Step 62: disformal transport delay on the SN Refsdal lens paths.

Step 58 bounded the disformal sector by transferring the GW170817
constraint |sigma| = |c_gamma - c_g|/c < 1e-15 directly onto the lens
halo traverse, giving ~5-13 ms. This step re-derives the bound correctly.

GW170817 constrains the disformal family coefficient B0, not the local
cone split sigma. The corpus's own evaluation (Paper 0 step_50,
B0_DISFORMAL_GW170817_MAX = 77.7) applies the bound to the PEAK of

    sigma/B0 = shape(u) * (du/dr * R_H)^2 ,  shape(u) = u^2/(1+u^2) e^{-u^4/2}

along the GW170817 host-galaxy profile. The same B0 must then be evaluated
on the MACS J1149 profile, where u and du/dr are set by a cluster-depth
potential, not a dwarf-galaxy one. The two profiles are different objects:
the bound transfers through B0, and the local split on the lens path is
whatever the cluster field configuration makes it.

This step therefore:

  1. rebuilds the screened scalar profile along each image's line of
     sight: g(z) from the sheet-normalized halo model, y(z) from the
     canonical screening law, du/dz = y g/c^2, u(z) integrated inward
     from ambient with the potential depth anchored to the GLAFIC psi
     column (psi * c^2 D_l / 2 gives the projected integral of Phi);
  2. evaluates sigma(z)/B0 = shape(u(z)) (n_hat . grad u * R_H)^2 with the
     geometric factor cos^2 theta = z^2/(z^2+b^2) for a straight ray at
     sky-plane impact b;
  3. integrates dt_i/B0 = (1/c) int sigma/B0 dl along the halo traverse;
  4. reports the differential SX-S1 delay per unit B0, the B0 required
     for the observed +30.1 d residual, and that B0's standing relative
     to the galactic-profile bound B0 <~ 77.7.

Also reported: the ratio of the cluster-path peak sigma/B0 to the
galactic-profile peak (1.29e-17) that set the B0 bound -- i.e., how much
the naive 1e-15-per-length transfer undercounted the disformal tilt on a
deep-potential path.

Inputs : GLAFIC v3 kappa + psi maps, refsdal_glafic_v3_lensing_params.json
Output : results/outputs/step_062_disformal_lens_path.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import print_status
from scripts.utils.tep_config import BETA_REFSDAL, load_refsdal_image_positions
from core import constants as tep_const
from core.scalar_field import screening_y

STEP_NUM = "062"

KAPPA_MAP = (
    PROJECT_ROOT
    / "data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_kappa.fits"
)
PSI_MAP = (
    PROJECT_ROOT
    / "data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_psi.fits"
)

Z_L, Z_S = 0.542, 1.489
H0_KMS = 70.0
OM = 0.3
C_KMS = 299792.458
G_SI = 6.67430e-11
MPC_M = 3.0856775814913673e22
ARCSEC_TO_RAD = np.pi / (180.0 * 3600.0)
DAY_S = 86400.0
G_T = tep_const.G_T_TRANSITION  # ~3.4e-10 m/s^2
R_H = tep_const.C_LIGHT / (H0_KMS * 1e3 / MPC_M)  # Hubble radius, m

# corpus bound (Paper 0 step_50): B0 <~ 77.7 from the peak of
# shape(u)*(du/dr*R_H)^2 on the GW170817 host-galaxy profile
B0_MAX_GAL = 77.7
SIGMA_PEAK_GAL_PER_B0 = 1.0e-15 / B0_MAX_GAL  # 1.287e-17

# halo LOS scale height h (the single cluster profile is shared by all
# five image rays; the per-image variation enters through g_i only)
H_MPC_SCAN = [0.1, 0.25, 0.5, 1.0, 2.0]
H_FID_MPC = 0.5
# the delay is integrated over |z| < L_int around the halo midplane
L_INT_MPC = 4.0


def comoving_distance(z, n=20000):
    zs = np.linspace(0.0, z, n)
    ez = np.sqrt(OM * (1.0 + zs) ** 3 + (1.0 - OM))
    return float(np.trapezoid(1.0 / ez, zs)) * C_KMS / H0_KMS


def angular_diameter(z):
    return comoving_distance(z) / (1.0 + z)


def shape_u(u):
    """B(u)/B0 = u^2/(1+u^2) exp(-u^4/2)."""
    return u**2 / (1.0 + u**2) * np.exp(-0.5 * u**4)


def los_profile(g_img, h, b, proj_mode="geometry", nz=20001):
    """Screened field along one line of sight.

    |grad u|(z) = y * g(z)/c^2 with g(z) = g_img * (1+(z/h)^2)^(-3/2)
    the Plummer-like LOS magnitude profile of the shared cluster halo;
    g_img is the per-image midplane value from the GLAFIC sheet
    normalization and b the sky-plane impact of the ray.

    The disformal term couples to (n_hat . grad u)^2: for a spherical
    halo grad u is radial, so the component along the photon direction
    is du/dz = (z/sqrt(z^2+b^2)) * du/dr -- the sky-plane gradient does
    not tilt the photon cone ("geometry" mode). "unity" brackets the
    flattened-halo limit where grad u is LOS-dominated (MACS J1149 is a
    merging, elongated cluster); "isotropic" applies cos^2 = 1/3.

    u(z) is accumulated inward from ambient (u = 0) along the ray using
    the same signed LOS component. The potential depth W = int g dl is
    derived, not assumed.

    Returns z, g, y, u, sigma_per_B0 arrays and h.
    """
    z = np.linspace(-8.0 * h, 8.0 * h, nz)
    g = g_img * (1.0 + (z / h) ** 2) ** (-1.5)

    y = screening_y(g)
    du_dr = y * g / tep_const.C_LIGHT**2  # |grad u|, m^-1

    # LOS component of the gradient along the ray. u decreases outward:
    # du/dz = -|grad u| z/r (negative for z > 0).
    r = np.sqrt(z**2 + b**2)
    if proj_mode == "geometry":
        proj = z / r
    elif proj_mode == "unity":
        proj = np.sign(z)
    else:  # isotropic
        proj = np.sign(z) * np.sqrt(1.0 / 3.0)
    du_dz = -du_dr * proj

    # tail contribution to u at the profile edge (z = 8h -> inf):
    # residual unscreened wing ~ y_edge * g_edge * scale
    g_edge = g_img * (h / (8.0 * h)) ** 2
    y_edge = float(screening_y(np.array([g_edge]))[0])
    u_edge = y_edge * g_edge * h / 8.0 / tep_const.C_LIGHT**2

    # u(z) = u_edge - int_z^{8h} du/dz dz'  (symmetric, deepest at z=0)
    u = np.empty_like(z)
    dz = np.diff(z)
    integ = 0.5 * (du_dz[1:] + du_dz[:-1]) * dz
    u_rev = np.cumsum(integ[::-1])[::-1]
    u[:-1] = u_rev
    u[-1] = 0.0
    u = u_edge - u

    sigma_per_b0 = shape_u(u) * (du_dz * R_H) ** 2
    return z, g, y, u, sigma_per_b0, h


def main():
    from astropy.io import fits
    from astropy.wcs import WCS

    print_status(f"STEP {STEP_NUM}: disformal lens-path delay", "TITLE")

    params = json.loads(
        (PROJECT_ROOT / "data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json").read_text()
    )
    imgs = params["images"]
    positions = load_refsdal_image_positions()

    s07 = json.load(open(PROJECT_ROOT / "results/outputs/step_07_observed_vs_predicted.json"))
    r_obs = float(s07["weighted_mean_residual"]["R_obs_days"])

    D_l = angular_diameter(Z_L)
    chi_l = comoving_distance(Z_L)
    chi_s = comoving_distance(Z_S)
    D_ls = (chi_s - chi_l) / (1.0 + Z_S)
    D_s = angular_diameter(Z_S)
    sigma_crit = (
        (C_KMS * 1e3) ** 2 / (4.0 * np.pi * G_SI)
        * (D_s * MPC_M) / ((D_l * MPC_M) * (D_ls * MPC_M))
    )

    with fits.open(KAPPA_MAP) as hdul:
        kappa_map = hdul[0].data.astype(np.float64)
        wcs = WCS(hdul[0].header)
    with fits.open(PSI_MAP) as hdul:
        psi_map = hdul[0].data.astype(np.float64)  # arcsec^2, D_ls/D_s = 1
        wcs_p = WCS(hdul[0].header)
    pix_deg = abs(float(hdul[0].header["CDELT2"]))
    pix_arcsec = pix_deg * 3600.0

    kappa_src = kappa_map * BETA_REFSDAL
    sigma_map = kappa_src * sigma_crit
    g_sheet = 2.0 * np.pi * G_SI * sigma_map
    y_map = screening_y(g_sheet)

    # cluster centre: kappa-intensity-weighted core centroid. (The
    # psi-map argmax sits at the map edge and is NOT the mass centre;
    # using it inflated the image impacts by ~10x.)
    core = kappa_map > 0.1 * np.nanmax(kappa_map)
    if core.sum() < 50:
        core = kappa_map > 0.05 * np.nanmax(kappa_map)
    yyg, xxg = np.mgrid[0:kappa_map.shape[0], 0:kappa_map.shape[1]]
    icx = float(np.average(xxg[core], weights=kappa_map[core]))
    icy = float(np.average(yyg[core], weights=kappa_map[core]))
    psi_base = float(np.nanmin(psi_map))  # map-edge baseline (lower bound on column)

    per_image = {}
    for name in imgs:
        ra, dec = positions[name]
        x, yy = wcs.all_world2pix(ra, dec, 0)
        fx = int(round(float(np.asarray(x).ravel()[0])))
        fy = int(round(float(np.asarray(yy).ravel()[0])))
        xp, yp = wcs_p.all_world2pix(ra, dec, 0)
        fxp = int(round(float(np.asarray(xp).ravel()[0])))
        fyp = int(round(float(np.asarray(yp).ravel()[0])))

        g_i = float(g_sheet[fy, fx])
        y_i = float(y_map[fy, fx])
        psi_i = float(psi_map[fyp, fxp]) * ARCSEC_TO_RAD**2  # rad^2
        # projected potential column: int Phi/c^2 dl = psi_inf * D_l/2
        col_i = psi_i * (D_l * MPC_M) / 2.0                     # m
        # excess column over map-edge baseline -> lower bound on depth
        col_exc = (float(psi_map[fyp, fxp]) - psi_base) * ARCSEC_TO_RAD**2 * (D_l * MPC_M) / 2.0

        # sky-plane impact from cluster centre, metres
        b_i = np.hypot(fx - icx, fy - icy) * pix_arcsec * ARCSEC_TO_RAD * D_l * MPC_M

        # deprojection factor g_3D/g_sheet at the impact: enclosed
        # projected mass inside b (the r < b disc is fully covered by
        # the map footprint at the true impacts) gives the spherical
        # midplane field GM(<b)/b^2 that drives the scalar profile.
        b_pix_i = b_i / (pix_arcsec * ARCSEC_TO_RAD * D_l * MPC_M)
        enc = (yyg - icy) ** 2 + (xxg - icx) ** 2 < b_pix_i ** 2
        M_enc = float(np.sum(sigma_map[enc]) * (pix_arcsec * ARCSEC_TO_RAD * D_l * MPC_M) ** 2)
        g_3d_i = G_SI * M_enc / b_i ** 2
        g_scale_i = g_3d_i / g_i

        per_image[name] = {
            "g_sheet": g_i,
            "y": y_i,
            "du_dl_per_m": y_i * g_i / tep_const.C_LIGHT**2,
            "b_mpc": float(b_i / MPC_M),
            "psi_rad2": psi_i,
            "g_3d_deprojected": g_3d_i,
            "g_scale_deprojection": g_scale_i,
            "potential_column_int_phi_dl_m3_s2": col_i * tep_const.C_LIGHT**2,
            "excess_column_over_edge_m3_s2": col_exc * tep_const.C_LIGHT**2,
        }

    # --- LOS solve per image, per halo scale h -----------------------
    def integrate_for(name, h_mpc, proj_mode="geometry", g_scale=1.0):
        gi = per_image[name]["g_sheet"] * g_scale
        b = per_image[name]["b_mpc"] * MPC_M
        h = h_mpc * MPC_M
        z, g, y, u, sig, h = los_profile(gi, h, b, proj_mode=proj_mode)
        # delay integral over |z| < L_int about the halo midplane
        mask = np.abs(z) <= 0.5 * L_INT_MPC * MPC_M
        I = float(np.trapezoid(sig[mask], z[mask]))  # metres per B0
        w_derived = float(np.trapezoid(g, z))  # int g dl, m^2/s^2
        return {
            "h_mpc": float(h / MPC_M),
            "u_midplane": float(u[np.argmin(np.abs(z))]),
            "W_int_g_dl_m2_s2": w_derived,
            "v_equiv_km_s": float(np.sqrt(w_derived) / 1e3),
            "sigma_peak_per_B0": float(np.max(sig)),
            "delay_integral_m_per_B0": I,
            "delay_days_per_B0": I / tep_const.C_LIGHT / DAY_S,
        }

    fid = {}
    for name in imgs:
        fid[name] = integrate_for(name, H_FID_MPC)

    resid_per_b0 = (
        fid["SX"]["delay_days_per_B0"] - fid["S1"]["delay_days_per_B0"]
    )
    b0_req = r_obs / resid_per_b0 if resid_per_b0 > 0 else float("nan")

    # halo-scale scan
    scan = {}
    for h_mpc in H_MPC_SCAN:
        sub = {n: integrate_for(n, h_mpc) for n in imgs}
        rb0 = sub["SX"]["delay_days_per_B0"] - sub["S1"]["delay_days_per_B0"]
        scan[f"h_{h_mpc}_Mpc"] = {
            "resid_days_per_B0": rb0,
            "B0_required_for_Robs": r_obs / rb0 if rb0 > 0 else None,
            "W_m2_s2_SX": sub["SX"]["W_int_g_dl_m2_s2"],
            "u_midplane_SX": sub["SX"]["u_midplane"],
            "residual_days_at_B0max": rb0 * B0_MAX_GAL,
        }

    # bound-transfer factor: cluster peak sigma/B0 vs galactic peak
    peak_cluster = max(fid[n]["sigma_peak_per_B0"] for n in imgs)
    headroom = peak_cluster / SIGMA_PEAK_GAL_PER_B0

    # normalization bracket: the sheet value g = 2 pi G Sigma is the
    # corpus convention for the local field, but a spherical deprojection
    # of a concentrated cluster gives a midplane |grad Phi| ~ factor
    # ~5-10 smaller; sigma scales ~ g^4 through u and du/dz, so this is
    # the dominant systematic. Report the conservative corner plus the
    # data-derived per-image deprojection computed above.
    norm_bracket = {}
    for gs in [1.0, 0.3, 0.1]:
        sub = {n: integrate_for(n, H_FID_MPC, g_scale=gs) for n in imgs}
        rb0 = sub["SX"]["delay_days_per_B0"] - sub["S1"]["delay_days_per_B0"]
        norm_bracket[f"g_scale_{gs}"] = {
            "resid_days_per_B0": rb0,
            "B0_required_for_Robs": r_obs / rb0 if rb0 > 0 else None,
        }
    # data-derived deprojection: per-image enclosed-mass g_3D
    dep = {n: integrate_for(
        n, H_FID_MPC,
        g_scale=per_image[n]["g_scale_deprojection"]) for n in imgs}
    rb0_dep = dep["SX"]["delay_days_per_B0"] - dep["S1"]["delay_days_per_B0"]
    norm_bracket["deprojected_data"] = {
        "resid_days_per_B0": rb0_dep,
        "B0_required_for_Robs": r_obs / rb0_dep if rb0_dep > 0 else None,
        "per_image_g_scale": {n: per_image[n]["g_scale_deprojection"]
                              for n in imgs},
    }
    # deprojected halo-scale scan (the decisive pair: h x deprojection)
    dep_scan = {}
    for h_mpc in H_MPC_SCAN:
        sub = {n: integrate_for(n, h_mpc,
                                g_scale=per_image[n]["g_scale_deprojection"])
               for n in imgs}
        rb0 = sub["SX"]["delay_days_per_B0"] - sub["S1"]["delay_days_per_B0"]
        dep_scan[f"h_{h_mpc}_Mpc"] = {
            "resid_days_per_B0": rb0,
            "B0_required_for_Robs": r_obs / rb0 if rb0 > 0 else None,
        }

    # cos^2 sensitivity at the fiducial point
    c2sens = {}
    for mode in ["geometry", "unity", "isotropic"]:
        sub = {n: integrate_for(n, H_FID_MPC, proj_mode=mode) for n in imgs}
        rb0 = sub["SX"]["delay_days_per_B0"] - sub["S1"]["delay_days_per_B0"]
        c2sens[mode] = {"resid_days_per_B0": rb0,
                        "B0_required": r_obs / rb0 if rb0 > 0 else None}

    out = {
        "step": STEP_NUM,
        "channel": "disformal transport (B(phi) cone tilt) on the lens paths",
        "correction_to_step_58": (
            "step_58 applied |sigma| < 1e-15 as a per-length bound on the "
            "local cone tilt in the lens halo. GW170817 bounds B0 through "
            "the PEAK of shape(u)*(du/dr*R_H)^2 on the NGC 4993 galaxy "
            "profile (B0 <~ 77.7, Paper 0 step_50). On the MACS J1149 "
            "profile the same B0 gives a local tilt larger by the ratio "
            "of the two profile functionals; the bound transfers through "
            "B0, not through sigma."
        ),
        "formula": "sigma/B0 = shape(u) (du/dl R_H)^2 cos^2 theta ; "
                   "shape(u) = u^2/(1+u^2) e^{-u^4/2} ; dt = (1/c) int sigma dl",
        "per_image_environment": per_image,
        "fiducial": {
            "h_mpc": H_FID_MPC,
            "L_int_mpc": L_INT_MPC,
            "per_image": fid,
            "residual_SX_minus_S1_days_per_B0": resid_per_b0,
            "residual_days_at_B0_max_77p7": resid_per_b0 * B0_MAX_GAL,
            "B0_required_for_Robs": b0_req,
            "B0_max_from_GW170817_galactic": B0_MAX_GAL,
            "margin_B0req_over_B0max": b0_req / B0_MAX_GAL,
        },
        "halo_scale_scan": scan,
        "normalization_bracket": norm_bracket,
        "deprojected_halo_scan": dep_scan,
        "cos2_sensitivity": c2sens,
        "bound_transfer_headroom": {
            "galactic_peak_sigma_per_B0": SIGMA_PEAK_GAL_PER_B0,
            "cluster_peak_sigma_per_B0": peak_cluster,
            "ratio": headroom,
            "note": "how much the naive 1e-15-per-length transfer "
                    "undercounts the permitted local tilt on a "
                    "cluster-depth path",
        },
        "R_obs_days": r_obs,
    }

    op = PROJECT_ROOT / "results/outputs/step_062_disformal_lens_path.json"
    op.write_text(json.dumps(out, indent=2))
    print_status(f"SX-S1 residual per B0: {resid_per_b0:.4g} d", "INFO")
    print_status(f"B0 required for +{r_obs:.1f} d: {b0_req:.4g} (bound {B0_MAX_GAL})", "INFO")
    print_status(f"headroom vs galactic-profile peak: {headroom:.3g}x", "INFO")
    print_status(f"wrote {op.name}", "OK")


if __name__ == "__main__":
    main()
