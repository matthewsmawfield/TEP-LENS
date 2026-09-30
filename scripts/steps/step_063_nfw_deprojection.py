#!/usr/bin/env python3
"""Step 63: decisive deprojection for the disformal lens-path delay.

Step 62 established that the disformal channel's viability is decided by
the halo's three-dimensional field profile. This step replaces the
Plummer + sheet-normalization shortcut with a data-determined deprojection
of the GLAFIC v3 convergence map:

  1. Cluster-scale centroid from the heavily smoothed kappa map
     (sigma = 10 arcsec), which washes out the host-galaxy substructure
     spike that contaminates the unsmoothed kappa-weighted centroid.
  2. Spherical NFW + constant-background fit to the azimuthal surface
     density profile about that centroid.
  3. Host-galaxy component at the S1-S4 cross centroid: compact-mass
     excess inside a 60 kpc aperture after subtracting the fitted smooth
     halo, modelled as a spherical Plummer halo of scale a_gal.
  4. Per image: the true 3D field along the ray,
        g_z(z) = sum over components of g_mag(r) z/r,
        |g_tot|  = |g_vec| (LOS + sky components added vectorially)
     du/dz = y(|g_tot|) g_z / c^2 accumulated inward from ambient,
     sigma(z)/B0 = shape(u) (du/dz R_H)^2, delay = (1/c) int sigma dl.

Validation: the model's projected convergence at the five image positions
is compared with the tabulated GLAFIC kappa values; the SX-deepest
ordering must survive for the model to be admissible.

Inputs : GLAFIC v3 kappa map, refsdal_glafic_v3_lensing_params.json
Output : results/outputs/step_063_nfw_deprojection.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.logger import print_status
from scripts.utils.tep_config import BETA_REFSDAL, load_refsdal_image_positions
from core import constants as tep_const
from core.scalar_field import screening_y

STEP_NUM = "063"

KAPPA_MAP = (
    PROJECT_ROOT
    / "data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_kappa.fits"
)

Z_L, Z_S = 0.542, 1.489
H0_KMS = 70.0
OM = 0.3
C_KMS = 299792.458
G_SI = 6.67430e-11
MPC_M = 3.0856775814913673e22
KPC_M = 3.0856775814913673e19
MSUN = 1.98847e30
ARCSEC_TO_RAD = np.pi / (180.0 * 3600.0)
DAY_S = 86400.0
G_T = tep_const.G_T_TRANSITION
R_H = tep_const.C_LIGHT / (H0_KMS * 1e3 / MPC_M)

B0_MAX_GAL = 77.7

# galaxy-component geometry (cross host; Plummer scale of the compact halo)
A_GAL_KPC = 15.0
APERTURE_GAL_KPC = 25.0   # tight aperture: captures the compact core spike,
                          # excludes the halo-dominated outer annulus
SMOOTH_SIGMA_ARCSEC = 10.0


def comoving_distance(z, n=20000):
    zs = np.linspace(0.0, z, n)
    ez = np.sqrt(OM * (1.0 + zs) ** 3 + (1.0 - OM))
    return float(np.trapezoid(1.0 / ez, zs)) * C_KMS / H0_KMS


def angular_diameter(z):
    return comoving_distance(z) / (1.0 + z)


def shape_u(u):
    return u**2 / (1.0 + u**2) * np.exp(-0.5 * u**4)


def nfw_F(x):
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    gt = x > 1.0
    lt = x < 1.0
    out[gt] = 2.0 * np.arctan(np.sqrt((x[gt] - 1) / (x[gt] + 1))) / np.sqrt(x[gt] ** 2 - 1)
    out[lt] = 2.0 * np.arctanh(np.sqrt((1 - x[lt]) / (x[lt] + 1))) / np.sqrt(1 - x[lt] ** 2)
    return out


def sigma_nfw(R, sig_s, r_s):
    """Sigma(R) = sig_s * (1 - F(x))/(x^2 - 1), sig_s = 2 rho_s r_s.
    The removable x = 1 pole is set to its l'Hopital limit 1/3."""
    x = np.asarray(R, dtype=float) / r_s
    out = np.empty_like(x)
    near = np.abs(x - 1.0) < 1e-4
    out[near] = 1.0 / 3.0
    far = ~near
    out[far] = (1.0 - nfw_F(x[far])) / (x[far] ** 2 - 1.0)
    return sig_s * out


def nfw_menc(r, rho_s, r_s):
    c = r / r_s
    return 4.0 * np.pi * rho_s * r_s**3 * (np.log(1.0 + c) - c / (1.0 + c))


def plum_mag(r, m_gal, a):
    """Plummer field magnitude GM r / (r^2 + a^2)^{3/2}."""
    return G_SI * m_gal * r / (r**2 + a**2) ** 1.5


def main():
    from astropy.io import fits
    from astropy.wcs import WCS

    print_status(f"STEP {STEP_NUM}: NFW+galaxy deprojection, disformal delay", "TITLE")

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
        pix_deg = abs(float(hdul[0].header["CDELT2"]))
    pix_arcsec = pix_deg * 3600.0
    m_per_pix = pix_arcsec * ARCSEC_TO_RAD * D_l * MPC_M
    pix_area_m2 = m_per_pix**2

    kappa_src = kappa_map * BETA_REFSDAL
    sigma_map = kappa_src * sigma_crit  # kg/m^2

    # --- cluster centroid: peak of the smoothed map -----------------------
    sm = gaussian_filter(kappa_map, SMOOTH_SIGMA_ARCSEC / pix_arcsec)
    icy_c, icx_c = np.unravel_index(np.argmax(sm), sm.shape)
    ra_c, dec_c = wcs.all_pix2world(icx_c, icy_c, 0)

    yyg, xxg = np.mgrid[0:kappa_map.shape[0], 0:kappa_map.shape[1]]
    r_m = np.hypot(xxg - icx_c, yyg - icy_c) * m_per_pix

    # --- azimuthal NFW fit -----------------------------------------------
    edges = np.array([40, 60, 80, 110, 150, 200, 270, 350, 450, 550]) * KPC_M
    rb, sb, nb = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (r_m >= lo) & (r_m < hi)
        rb.append(0.5 * (lo + hi))
        sb.append(float(np.nanmean(sigma_map[m])))
        nb.append(int(m.sum()))
    rb, sb, nb = np.array(rb), np.array(sb), np.array(nb)

    def model(R, sig_s, r_s, bg):
        return sigma_nfw(R, sig_s, r_s) + bg

    def fit_nfw(r_s_fixed=None):
        if r_s_fixed is None:
            popt, _ = curve_fit(model, rb, sb, p0=[3.0, 3.0e22, 0.3 * float(np.min(sb))],
                                bounds=([0, 1e21, -1.0], [100.0, 1e24, 5.0]), maxfev=20000)
            return popt
        popt, _ = curve_fit(lambda R, s, b: model(R, s, r_s_fixed, b), rb, sb,
                            p0=[3.0, 0.3 * float(np.min(sb))],
                            bounds=([0, -1.0], [500.0, 5.0]), maxfev=20000)
        return np.array([popt[0], r_s_fixed, popt[1]])

    popt = fit_nfw()
    sig_s, r_s_m, bg = popt
    rho_s = sig_s / (2.0 * r_s_m)
    resid_frac = float(np.std((model(rb, *popt) - sb) / sb))

    # --- host-galaxy mass excess at the cross centroid --------------------
    pos_arr = {n: np.array(wcs.all_world2pix(*positions[n], 0)).ravel() for n in imgs}
    cross_c = np.mean([pos_arr[n] for n in ["S1", "S2", "S3", "S4"]], axis=0)
    ap_pix = APERTURE_GAL_KPC * KPC_M / m_per_pix
    ap_mask = (xxg - cross_c[0]) ** 2 + (yyg - cross_c[1]) ** 2 < ap_pix**2
    m_ap = float(np.sum(sigma_map[ap_mask]) * pix_area_m2 / MSUN)
    m_nfw_ap = float(np.sum(model(r_m[ap_mask], *popt)) * pix_area_m2 / MSUN)
    m_gal = max(m_ap - m_nfw_ap, 0.0) * MSUN
    a_gal = A_GAL_KPC * KPC_M

    b_cl = {n: float(np.hypot(*(pos_arr[n] - np.array([icx_c, icy_c]))) * m_per_pix)
            for n in imgs}
    b_gal = {n: float(np.hypot(*(pos_arr[n] - cross_c)) * m_per_pix) for n in imgs}

    # model kappa at image positions (validation vs tabulated values)
    def sigma_plummer(R, m_g, a):
        return (m_g / np.pi) * a**2 / (R**2 + a**2) ** 2

    kappa_model = {}
    for n in imgs:
        sig_tot = sigma_nfw(b_cl[n], sig_s, r_s_m) + bg + sigma_plummer(b_gal[n], m_gal, a_gal)
        kappa_model[n] = float(sig_tot / sigma_crit)  # physical kappa at z_s

    # --- per-image LOS solve ---------------------------------------------
    def sky_unit(p):
        e_cl = p - np.array([icx_c, icy_c]); e_cl /= np.linalg.norm(e_cl)
        e_g = p - cross_c; e_g /= np.linalg.norm(e_g)
        return e_cl, e_g

    NZ = 40001

    def solve_ray(name, rho_s_v, r_s_v, m_gal_v, u_amb=0.0):
        b_c, b_g = b_cl[name], b_gal[name]
        e_cl, e_g = sky_unit(pos_arr[name])
        hmax = max(5.0 * r_s_v, 8.0 * (b_c + b_g))
        z = np.linspace(-hmax, hmax, NZ)
        r_cl = np.sqrt(z**2 + b_c**2)
        r_g = np.sqrt(z**2 + b_g**2)
        gm_cl = G_SI * nfw_menc(r_cl, rho_s_v, r_s_v) / r_cl**2
        gm_g = plum_mag(r_g, m_gal_v, a_gal)
        g_z = gm_cl * (z / r_cl) + gm_g * (z / r_g)
        sky2 = ((gm_cl * b_c / r_cl)[None, :] * e_cl[:, None]
                + (gm_g * b_g / r_g)[None, :] * e_g[:, None])
        g_tot = np.sqrt(g_z**2 + np.linalg.norm(sky2, axis=0) ** 2)
        y = screening_y(g_tot)
        du_dz = -y * g_z / tep_const.C_LIGHT**2
        dz = np.diff(z)
        integ = 0.5 * (du_dz[1:] + du_dz[:-1]) * dz
        u_rev = np.cumsum(integ[::-1])[::-1]          # int_z^top du_dz
        u = np.empty_like(z)
        u[:-1] = -u_rev                              # halo perturbation
        u[-1] = 0.0
        u_tot = u + u_amb                            # + cosmic ambient u_bar(z_l)
        sig = shape_u(u_tot) * (np.abs(du_dz) * R_H) ** 2
        I = float(np.trapezoid(sig, z))
        return {
            "b_cluster_kpc": b_c / KPC_M,
            "b_galaxy_kpc": b_g / KPC_M,
            "u_midplane": float(u[np.argmin(np.abs(z))]),
            "sigma_peak_per_B0": float(np.max(sig)),
            "delay_days_per_B0": I / tep_const.C_LIGHT / DAY_S,
        }

    def run_component_model(rho_s_v, r_s_v, m_gal_v, u_amb=0.0):
        pi = {n: solve_ray(n, rho_s_v, r_s_v, m_gal_v, u_amb) for n in imgs}
        rp = pi["SX"]["delay_days_per_B0"] - pi["S1"]["delay_days_per_B0"]
        return pi, rp, (r_obs / rp if rp > 0 else float("nan"))

    per_img, resid_pb, b0_req = run_component_model(rho_s, r_s_m, m_gal)
    _, resid_pb_ng, b0_ng = run_component_model(rho_s, r_s_m, 0.0)

    # Ambient-consistent variant: B(u) is a function of the TOTAL field.
    # Corpus convention (Rule 4): phi_bar(z) = ln(1+z) M_Pl, so the
    # ambient field at the lens plane is u_bar(z_l) = 0.433, and
    # shape(u) on the path is ~0.155, not ~1e-8. The same convention
    # applied to the GW170817 path adds the cosmic-drift term
    # (du_bar/dx)R_H = 1 that the step_50 galactic-profile bound omits:
    # |Dc/c| ~ B0 * int shape(u_bar) du_bar over the path ~ B0*3e-5
    # -> B0_max_drift ~ 3e-11. Recorded for both-side consistency.
    u_amb_l = float(np.log(1.0 + Z_L))
    _, resid_pb_amb, b0_amb = run_component_model(rho_s, r_s_m, m_gal,
                                                u_amb=u_amb_l)
    u_gw = float(np.log(1.0 + 0.0097))
    i_drift = float(
        np.trapezoid(shape_u(np.linspace(0, u_gw, 4000)),
                     np.linspace(0, u_gw, 4000)))
    frac_drift = i_drift * tep_const.C_LIGHT / (
        (H0_KMS * 1e3 / MPC_M) * 40e6 * 3.0856775814913673e16)
    b0_max_drift = 1.0e-15 / frac_drift

    # ---- admissible-form variant: B_eff = B(phi_tot) * G(delta_u) --
    # The corpus's only admissible disformal form (registry, 29-4;
    # Paper 29 gate10b): the gate is closed on the low-X ambient and
    # saturates inside nonlinear wells above the DLA depth scale
    # (~1e-7), with population consistency bounding G(1e-5) <~
    # few x 1e-3/u_bar. The cluster well (delta u ~ 1e-4) is deep in
    # the saturated regime, and B(phi) reads the total field, so
    # shape(u_bar + delta_u) ~ 0.155 applies inside the halo while
    # the ambient/drift bounds are evaded by gate closure.
    # Gate model: G(u) = g_sat * min(1, u/u_act) in the local
    # perturbation u (well depth above ambient); bracket the corpus's
    # inferred saturation level.
    # Activation depth: gate10b requires the DLA wells (delta_phi ~
    # 2e-8--1.5e-7) to sit inside the activation transition to supply
    # the measured dG/dln(dphi) ~ 4--7e-4 slopes; U_ACT ~ 1e-7 places
    # them there. (1e-6 would leave the DLAs only ~2--15% open, short
    # of the absorber contrast.) The cluster (1e-4) and the GW host
    # (5e-7) are both saturated at this scale.
    U_ACT = 1.0e-7

    def solve_ray_gated(name, rho_s_v, r_s_v, m_gal_v, g_sat,
                        u_act=U_ACT):
        b_c, b_g = b_cl[name], b_gal[name]
        e_cl, e_g = sky_unit(pos_arr[name])
        hmax = max(5.0 * r_s_v, 8.0 * (b_c + b_g))
        z = np.linspace(-hmax, hmax, NZ)
        r_cl = np.sqrt(z**2 + b_c**2)
        r_g = np.sqrt(z**2 + b_g**2)
        gm_cl = G_SI * nfw_menc(r_cl, rho_s_v, r_s_v) / r_cl**2
        gm_g = plum_mag(r_g, m_gal_v, a_gal)
        g_z = gm_cl * (z / r_cl) + gm_g * (z / r_g)
        sky2 = ((gm_cl * b_c / r_cl)[None, :] * e_cl[:, None]
                + (gm_g * b_g / r_g)[None, :] * e_g[:, None])
        g_tot = np.sqrt(g_z**2 + np.linalg.norm(sky2, axis=0) ** 2)
        y = screening_y(g_tot)
        du_dz = -y * g_z / tep_const.C_LIGHT**2
        dz = np.diff(z)
        integ = 0.5 * (du_dz[1:] + du_dz[:-1]) * dz
        u_rev = np.cumsum(integ[::-1])[::-1]
        u = np.empty_like(z)
        u[:-1] = -u_rev
        u[-1] = 0.0
        gate = g_sat * np.minimum(1.0, u / u_act)
        sig = (shape_u(u + u_amb_l) * gate
               * (np.abs(du_dz) * R_H) ** 2)
        return float(np.trapezoid(sig, z))

    # matching bound: same gated form on the step_50 GW170817 host
    # profile (uniform sphere M=1e11 Msun, R=30 kpc; host well depth
    # ~5e-7 -- borderline saturated)
    KPC_M_ = 3.0856775814913673e19
    rgw0 = 30.0 * KPC_M_          # 30 kpc host radius (step_50 convention)
    r_gw = np.geomspace(1e-3 * rgw0, 50.0 * rgw0, 6000)
    mgw = 1.0e11 * MSUN
    ggw = np.where(r_gw < rgw0, G_SI * mgw * r_gw / rgw0**3,
                   G_SI * mgw / r_gw**2)
    ygw = screening_y(ggw)
    dudl = ygw * ggw / tep_const.C_LIGHT**2
    drg = np.diff(r_gw)
    ugw = np.zeros_like(r_gw)
    ugw[:-1] = np.cumsum((0.5 * (dudl[1:] + dudl[:-1]) * drg)[::-1])[::-1]

    gated_rows = {}
    for g_sat in (1.0e-3, 3.0e-3, 7.0e-3):
        ig = {n: solve_ray_gated(n, rho_s, r_s_m, m_gal, g_sat)
              for n in imgs}
        rp_g = ig["SX"] - ig["S1"]
        b0_req_g = r_obs / (rp_g / tep_const.C_LIGHT / DAY_S)
        gate_gw = g_sat * np.minimum(1.0, ugw / U_ACT)
        dbar_gw = shape_u(ugw + u_gw) * gate_gw * (dudl * R_H) ** 2
        b0_max_gw = 1.0e-15 / float(np.max(dbar_gw))
        gated_rows[f"{g_sat:.0e}"] = {
            "G_sat": g_sat,
            "resid_days_per_B0": rp_g / tep_const.C_LIGHT / DAY_S,
            "B0_required": b0_req_g,
            "B0_max_GW_same_convention": b0_max_gw,
            "margin": b0_max_gw / b0_req_g,
        }

    # activation-scale robustness: margin = B0_max/B0_req scales with
    # gate_cluster/gate_host. Scan U_ACT to show no admissible value
    # rescues the channel -- the g_sat factor cancels identically.
    u_act_scan = {}
    for u_a in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
        ig_u = {n: solve_ray_gated(n, rho_s, r_s_m, m_gal, 1.0e-3, u_a)
                for n in imgs}
        b0_req_u = r_obs / ((ig_u["SX"] - ig_u["S1"]) / tep_const.C_LIGHT / DAY_S)
        gate_gw_u = 1.0e-3 * np.minimum(1.0, ugw / u_a)
        b0_max_u = 1.0e-15 / float(
            np.max(shape_u(ugw + u_gw) * gate_gw_u * (dudl * R_H) ** 2))
        u_act_scan[f"{u_a:.0e}"] = {
            "U_ACT": u_a,
            "B0_required": b0_req_u,
            "B0_max_GW": b0_max_u,
            "margin": b0_max_u / b0_req_u,
            "dla_gate_open_fraction_at_5e-8": float(min(1.0, 5e-8 / u_a)),
        }

    # halo-concentration sensitivity: refit (sig_s, bg) at fixed r_s and
    # re-solve -- the u^2 gate makes the answer concentration-sensitive
    rs_scan = {}
    for rs_kpc in [200.0, 300.0, 500.0, float(r_s_m / KPC_M), 1000.0, 1500.0,
                   2000.0, 3000.0, 5000.0]:
        pp = fit_nfw(rs_kpc * KPC_M)
        rho_v = pp[0] / (2.0 * pp[1])
        # refit galaxy excess under this smooth model
        m_nfw_ap_v = float(np.sum(model(r_m[ap_mask], *pp)) * pix_area_m2 / MSUN)
        mg_v = max(m_ap - m_nfw_ap_v, 0.0) * MSUN
        _, rp_v, b0_v = run_component_model(rho_v, pp[1], mg_v)
        rs_scan[f"{rs_kpc:.0f}"] = {
            "sig_s": float(pp[0]), "bg": float(pp[2]), "m_gal_msun": float(mg_v / MSUN),
            "resid_days_per_B0": rp_v, "B0_required": b0_v,
            "scatter": float(np.std((model(rb, *pp) - sb) / sb)),
        }

    out = {
        "step": STEP_NUM,
        "channel": "disformal delay on deprojected NFW + host-galaxy field",
        "centroid": {
            "method": f"argmax of gaussian-smoothed kappa (sigma={SMOOTH_SIGMA_ARCSEC} arcsec)",
            "ra_deg": float(ra_c), "dec_deg": float(dec_c),
        },
        "nfw_fit": {
            "sig_s_kg_m2": float(sig_s), "r_s_kpc": float(r_s_m / KPC_M),
            "rho_s_kg_m3": float(rho_s), "background_kg_m2": float(bg),
            "profile_frac_scatter": resid_frac,
            "radii_kpc": (rb / KPC_M).tolist(),
            "sigma_kg_m2": sb.tolist(), "model": model(rb, *popt).tolist(),
        },
        "galaxy": {
            "aperture_kpc": APERTURE_GAL_KPC, "plummer_a_kpc": A_GAL_KPC,
            "m_aperture_msun": m_ap, "m_smooth_in_aperture_msun": m_nfw_ap,
            "m_gal_excess_msun": float(m_gal / MSUN),
        },
        "kappa_validation": {
            n: {"model": kappa_model[n], "tabulated": params["images"][n]["kappa"]}
            for n in imgs
        },
        "per_image": per_img,
        "residual_SX_minus_S1_days_per_B0": resid_pb,
        "residual_days_at_B0max": resid_pb * B0_MAX_GAL,
        "B0_required_for_Robs": b0_req,
        "B0_max": B0_MAX_GAL,
        "no_galaxy_variant": {
            "residual_days_per_B0": resid_pb_ng,
            "B0_required": b0_ng,
        },
        "r_s_scan_kpc": rs_scan,
        "ambient_consistent_variant": {
            "u_bar_z_l": u_amb_l,
            "shape_at_ambient": float(shape_u(u_amb_l)),
            "resid_days_per_B0": resid_pb_amb,
            "B0_required": b0_amb,
            "B0_max_drift_term_GW170817": b0_max_drift,
            "signature_wall_B0_at_u_amb": float(1.0 / shape_u(u_amb_l)),
            "note": "total-field convention: B(u) evaluated at u_bar(z_l)+u_halo; "
                    "drift term (du_bar/dx)R_H=1 tightens the multimessenger bound "
                    "to ~3e-11; signature wall at z_l ambient is B0 <~ 6.5",
        },
        "admissible_gated_form": {
            "model": "B_eff = B(phi_tot)*G(delta_u); G closed on ambient, "
                     "saturating above the DLA depth scale (gate10b). "
                     "G(u) = g_sat*min(1, u/U_ACT), U_ACT=1e-7",
            "u_act": U_ACT,
            "G_sat_bracket": gated_rows,
            "u_act_scan": u_act_scan,
            "note": "corpus admissible coupling form (registry 29-4); the "
                    "uniform-B0 exclusion above does not apply to this form "
                    "because the ambient/drift bounds are evaded by gate "
                    "closure, while the GW bound under the same gated "
                    "convention is evaluated on the step_50 host profile. "
                    "Corrected margin: the same-convention bound is "
                    "B0_max ~ 2e-5 (host well saturated at the "
                    "gate10b-consistent U_ACT=1e-7), so the gated channel "
                    "is also excluded by ~10^3; the u_act_scan shows no "
                    "activation scale rescues it (best-case margin ~0.5 "
                    "requires U_ACT ~ 1e-4, which closes the DLA wells and "
                    "destroys the absorber-sector fit).",
        },
        "R_obs_days": r_obs,
        "verdict_basis": "B0_required vs 77.7 galactic-profile bound",
    }
    op = PROJECT_ROOT / "results/outputs/step_063_nfw_deprojection.json"
    op.write_text(json.dumps(out, indent=2))
    print_status(
        f"NFW: r_s={r_s_m/KPC_M:.0f} kpc, rho_s={rho_s:.2e}, bg={bg:.3f}, scatter={resid_frac:.2f}",
        "INFO",
    )
    print_status(f"M_gal excess = {m_gal/MSUN:.3e} Msun", "INFO")
    for n in imgs:
        print_status(
            f"{n}: b_cl={b_cl[n]/KPC_M:.1f} kpc b_gal={b_gal[n]/KPC_M:.1f} kpc "
            f"k_mod={kappa_model[n]:.2f} k_tab={params['images'][n]['kappa']:.2f} "
            f"u_mid={per_img[n]['u_midplane']:.2e} dt/B0={per_img[n]['delay_days_per_B0']:.3f} d",
            "INFO",
        )
    print_status(
        f"SX-S1 residual per B0: {resid_pb:.4g} d -> B0_req = {b0_req:.4g} (bound {B0_MAX_GAL})",
        "INFO",
    )
    print_status(f"wrote {op.name}", "OK")


if __name__ == "__main__":
    main()
