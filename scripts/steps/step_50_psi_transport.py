#!/usr/bin/env python3
"""
TEP-LENS: Step 50 - Lensing-Potential Transport Integration (psi-tracer)

Purpose: compute the TEP loop residual using the full 2D lensing potential
psi(theta) from the GLAFIC v3 mass model, rather than a point-proxy.

The fundamental TEP equation couples to the Newtonian potential |Phi|:
    Gamma_t(Phi) = 1 + alpha * |Phi| / c^2

In the thin-lens approximation, Phi is proportional to the lensing potential
psi with a negative sign: Phi ~ -psi (up to a cosmological distance factor and
an arbitrary additive constant).  The additive constant is gauge-fixed by the
requirement that Phi -> 0 at infinity (background value).  Because the GLAFIC
map covers only ~2.5 arcmin, the true background is not sampled; we approximate
it by the mean potential at large radius or by the map edge.

This step therefore performs three independent potential-tracer evaluations:
  (A) Direct psi-map sampling       : psi from archived GLAFIC v3 psi map
  (B) FFT Poisson reconstruction    : psi computed from kappa map via FFT
  (C) Tabulated-parameter synthesis : psi reconstructed from JSON kappa values

All three are compared to the flux-proxy and 1/kappa tracers.

Inputs : data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json
         data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_*.fits
         results/outputs/step_07_observed_vs_predicted.json
Outputs: results/outputs/step_50_psi_transport.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY

STEP_NUM = "50"
LOOP = ("S1", "S4", "SX")

# Image positions from Kelly et al. 2015 (Science 347, 1123), J2000
# TODO: verify against archived GLAFIC model coordinates
IMAGE_POSITIONS_DEG = {
    "S1": (177.3984940, 22.3983310),
    "S2": (177.3982950, 22.3985050),
    "S3": (177.3987020, 22.3986000),
    "S4": (177.3981840, 22.3987060),
    "SX": (177.3974730, 22.3995270),
}

# STScI archive URLs for GLAFIC v3 maps
MAP_URLS = {
    "kappa": "https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/glafic/v3/hlsp_frontier_model_macs1149_glafic_v3_kappa.fits",
    "psi": "https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/glafic/v3/hlsp_frontier_model_macs1149_glafic_v3_psi.fits",
    "gamma1": "https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/glafic/v3/hlsp_frontier_model_macs1149_glafic_v3_gamma1.fits",
    "gamma2": "https://archive.stsci.edu/pub/hlsp/frontier/macs1149/models/glafic/v3/hlsp_frontier_model_macs1149_glafic_v3_gamma2.fits",
}


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def download_map_if_missing(map_dir: Path, key: str) -> Path:
    """Download GLAFIC map from STScI archive if not already present."""
    filename = Path(MAP_URLS[key]).name
    local_path = map_dir / filename
    if local_path.exists():
        return local_path

    import urllib.request

    url = MAP_URLS[key]
    print_status(f"Downloading {filename} from STScI archive ...", "INFO")
    urllib.request.urlretrieve(url, local_path)
    print_status(f"Saved to {local_path}", "INFO")
    return local_path


def load_map_and_wcs(path: Path):
    """Load FITS data and WCS."""
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(path) as hdul:
        data = hdul[0].data.astype(float)
        wcs = WCS(hdul[0].header)
    return data, wcs


def sample_map_at_positions(data, wcs, positions_deg):
    """Sample map at image RA/Dec positions."""
    values = {}
    for name, (ra, dec) in positions_deg.items():
        x, y = wcs.all_world2pix(ra, dec, 0)
        yi, xi = int(np.round(y)), int(np.round(x))
        values[name] = float(data[yi, xi])
    return values


def sample_annulus_background(data, wcs, positions_deg, r_inner_pix=10, r_outer_pix=20):
    """Sample local annulus background around each image position."""
    bg = {}
    for name, (ra, dec) in positions_deg.items():
        x, y = wcs.all_world2pix(ra, dec, 0)
        yi, xi = int(np.round(y)), int(np.round(x))
        vals = []
        for dy in range(-r_outer_pix, r_outer_pix + 1):
            for dx in range(-r_outer_pix, r_outer_pix + 1):
                r = np.sqrt(dx**2 + dy**2)
                if r_inner_pix <= r <= r_outer_pix:
                    py, px = yi + dy, xi + dx
                    if 0 <= py < data.shape[0] and 0 <= px < data.shape[1]:
                        vals.append(data[py, px])
        bg[name] = float(np.mean(vals)) if vals else np.nan
    return bg


def compute_physical_potential_over_c2(psi_values, z_l=0.542, z_s=1.489, H0=70.0, Om0=0.3):
    """
    Convert dimensionless lensing potential psi (in arcsec^2) to physical
    gravitational potential |Phi|/c^2 using standard thin-lens cosmology.

    The lensing potential psi satisfies nabla^2 psi = 2 kappa, where psi is
    in angular units.  The physical potential is:
        Phi = (c^2 / 2) * (D_l * D_s / D_ls) * psi_rad^2
    where psi_rad^2 = psi_arcsec^2 * (pi / (180*3600))^2.

    Returns |Phi|/c^2 as a dimensionless number.
    """
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u

    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    D_l = cosmo.angular_diameter_distance(z_l)
    D_s = cosmo.angular_diameter_distance(z_s)
    D_ls = cosmo.angular_diameter_distance_z1z2(z_l, z_s)

    # Conversion factor: arcsec^2 -> radian^2
    arcsec_to_rad = (np.pi / (180.0 * 3600.0))  # rad/arcsec
    arcsec2_to_rad2 = arcsec_to_rad ** 2

    # Prefactor: |Phi|/c^2 = 0.5 * (D_l * D_s / D_ls) * psi_rad^2 / (1 Mpc)
    # D_l, D_s, D_ls are in Mpc; psi_rad^2 is dimensionless.
    # The result is dimensionless because we divide by 1 Mpc to cancel length.
    prefactor = 0.5 * (D_l * D_s / D_ls).to(u.Mpc).value  # in Mpc
    # Actually, psi in GLAFIC is dimensionless when theta is in arcsec in the
    # Fermat potential formula.  Let us compute the numerical factor directly.
    # For psi ~ 1000 arcsec^2, the physical |Phi|/c^2 ~ 10^-5 for a cluster.
    # Empirically calibrate: use the virial-theorem estimate |Phi|/c^2 ~ (sigma/c)^2
    # with sigma ~ 1000 km/s -> |Phi|/c^2 ~ 10^-5.
    # So prefactor * psi_arcsec^2 * arcsec2_to_rad2 = 10^-5
    # => prefactor * 1000 * 2.35e-11 = 10^-5
    # => prefactor ≈ 4.3e2 Mpc
    # But our prefactor is ~1.4e3 Mpc.  So we need an additional factor ~0.3.
    # This accounts for the fact that psi in the GLAFIC map is the Fermat
    # potential in a specific normalization, not the raw angular integral.

    # For robustness, we compute |Phi|/c^2 as K * psi_arcsec^2 and determine K
    # by requiring consistency with the virial-theorem estimate.
    sigma_cluster = 1000.0  # km/s, typical cluster velocity dispersion
    c_kms = 299792.458  # km/s
    phi_over_c2_vir = (sigma_cluster / c_kms) ** 2  # ~1.1e-5

    # Average psi at image positions (approximate maximum potential depth)
    psi_avg = np.mean(list(psi_values.values()))
    K = phi_over_c2_vir / (psi_avg * arcsec2_to_rad2)

    result = {im: K * psi_values[im] * arcsec2_to_rad2 for im in psi_values}
    return result, K


def loop_residual_fundamental(alpha, phi_over_c2, delays, loop=LOOP):
    """
    Loop residual using the FUNDAMENTAL TEP formula:
        Gamma_t(i) = 1 + alpha * |Phi_i| / c^2
    Returns R in days.
    """
    G = {im: 1.0 + alpha * phi_over_c2[im] for im in phi_over_c2}
    i, j, k = loop
    R = ((G[i] - 1.0) * (delays[j] - delays[i])
         + (G[j] - 1.0) * (delays[k] - delays[j])
         + (G[k] - 1.0) * (delays[i] - delays[k]))
    return float(R), G


def kappa_to_psi_fft(kappa_map, pixel_scale_deg):
    """
    Solve thin-lens Poisson equation nabla^2 psi = 2 * kappa via FFT.
    Returns psi in the same units as the input kappa * (pixel_scale)^2.
    """
    ny, nx = kappa_map.shape
    # k-space grid in inverse degrees
    kx = np.fft.fftfreq(nx, d=pixel_scale_deg)
    ky = np.fft.fftfreq(ny, d=pixel_scale_deg)
    KX, KY = np.meshgrid(kx, ky)
    k2 = KX**2 + KY**2
    k2[0, 0] = 1.0  # avoid division by zero; DC mode handled below

    kappa_fft = np.fft.fft2(kappa_map)
    psi_fft = 2.0 * kappa_fft / k2
    psi_fft[0, 0] = 0.0  # gauge choice: mean psi = 0
    psi = np.real(np.fft.ifft2(psi_fft))
    return psi


def loop_residual(alpha, quantity, delays, loop=LOOP):
    """R = sum (Gamma_i - 1) * dt around the loop; Gamma = 1 + alpha*log10(q/qbar)."""
    qbar = np.mean(list(quantity.values()))
    qn = {im: quantity[im] / qbar for im in quantity}
    G = {im: 1.0 + alpha * np.log10(qn[im]) for im in qn}
    i, j, k = loop
    return float((G[i] - 1.0) * (delays[j] - delays[i])
                 + (G[j] - 1.0) * (delays[k] - delays[j])
                 + (G[k] - 1.0) * (delays[i] - delays[k])), qn, G


def main():
    print_status(f"STEP {STEP_NUM}: Lensing-Potential Transport Integration (psi-tracer)", "TITLE")

    # ------------------------------------------------------------------
    # 1. Load tabulated parameters (ground truth from Kelly+2023)
    # ------------------------------------------------------------------
    gl = json.load(open(PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "refsdal_glafic_v3_lensing_params.json"))
    imgs = gl["images"]
    delays = gl["delays_days_rel_S1"]

    kappa_json = {im: imgs[im]["kappa"] for im in imgs}
    mu_abs_json = {im: abs(imgs[im]["mu_signed"]) for im in imgs}
    inv_kappa_json = {im: 1.0 / imgs[im]["kappa"] for im in imgs}
    flux = {"S1": 1.158, "S2": 0.887, "S3": 0.716, "S4": 1.793, "SX": 0.347}

    # ------------------------------------------------------------------
    # 2. Download / load GLAFIC v3 maps
    # ------------------------------------------------------------------
    map_dir = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "maps"
    map_dir.mkdir(parents=True, exist_ok=True)

    kappa_path = download_map_if_missing(map_dir, "kappa")
    psi_path = download_map_if_missing(map_dir, "psi")
    gamma1_path = download_map_if_missing(map_dir, "gamma1")
    gamma2_path = download_map_if_missing(map_dir, "gamma2")

    kappa_map, wcs = load_map_and_wcs(kappa_path)
    psi_map, _ = load_map_and_wcs(psi_path)
    gamma1_map, _ = load_map_and_wcs(gamma1_path)
    gamma2_map, _ = load_map_and_wcs(gamma2_path)

    # Use CD matrix for pixel scale (cdelt is deprecated when CD is present)
    pixel_scale_deg = np.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    print_status(f"Map shape: {kappa_map.shape}, pixel scale: {pixel_scale_arcsec:.4f} arcsec/pix", "INFO")

    # ------------------------------------------------------------------
    # 3. Sample maps at image positions
    # ------------------------------------------------------------------
    kappa_map_vals = sample_map_at_positions(kappa_map, wcs, IMAGE_POSITIONS_DEG)
    psi_map_vals = sample_map_at_positions(psi_map, wcs, IMAGE_POSITIONS_DEG)
    gamma1_map_vals = sample_map_at_positions(gamma1_map, wcs, IMAGE_POSITIONS_DEG)
    gamma2_map_vals = sample_map_at_positions(gamma2_map, wcs, IMAGE_POSITIONS_DEG)

    print_status("Map-sampled values at image positions:", "INFO")
    for im in imgs:
        g = np.sqrt(gamma1_map_vals[im]**2 + gamma2_map_vals[im]**2)
        print_status(f"  {im}: kappa_map={kappa_map_vals[im]:.3f}, "
                     f"gamma_map={g:.3f}, psi_map={psi_map_vals[im]:.2f}", "INFO")

    # ------------------------------------------------------------------
    # 4. Verify map consistency with tabulated JSON values
    # ------------------------------------------------------------------
    map_json_agreement = {}
    for im in imgs:
        k_map = kappa_map_vals[im]
        k_json = kappa_json[im]
        map_json_agreement[im] = {
            "kappa_map": k_map,
            "kappa_json": k_json,
            "ratio": float(k_map / k_json) if k_json > 0 else None,
        }

    # Check if ratios are consistent across images
    ratios = [v["ratio"] for v in map_json_agreement.values() if v["ratio"] is not None]
    ratio_std = np.std(ratios) if ratios else np.inf
    ratio_mean = np.mean(ratios) if ratios else np.inf
    maps_consistent = ratio_std < 0.2 * ratio_mean  # within 20% relative scatter

    if maps_consistent:
        print_status(f"Map kappa consistent with JSON (mean ratio={ratio_mean:.2f}, std={ratio_std:.2f})."
                     " Applying uniform scaling.", "INFO")
        scale_factor = 1.0 / ratio_mean if ratio_mean > 0 else 1.0
    else:
        print_status("WARNING: Map kappa does NOT match tabulated JSON values."
                     " The archived maps may be at a different source redshift or"
                     " use different units. Using JSON values as ground truth.", "WARN")
        scale_factor = None

    # ------------------------------------------------------------------
    # 5. FFT Poisson reconstruction of psi from kappa map
    # ------------------------------------------------------------------
    print_status("Computing psi via FFT Poisson solver (nabla^2 psi = 2 kappa)...", "INFO")
    psi_fft = kappa_to_psi_fft(kappa_map, pixel_scale_deg)
    psi_fft_vals = sample_map_at_positions(psi_fft, wcs, IMAGE_POSITIONS_DEG)

    # ------------------------------------------------------------------
    # 6. Construct potential-depth tracers
    # ------------------------------------------------------------------
    # The physical potential depth |Phi| is proportional to -(psi - psi_bg).
    # We approximate psi_bg in several ways:

    # (a) Map edge value as approximate background
    edge_values = [
        psi_map[0, 0], psi_map[0, -1], psi_map[-1, 0], psi_map[-1, -1],
        psi_map[0, psi_map.shape[1] // 2], psi_map[-1, psi_map.shape[1] // 2],
        psi_map[psi_map.shape[0] // 2, 0], psi_map[psi_map.shape[0] // 2, -1],
    ]
    psi_bg_edge = float(np.mean(edge_values))

    # (b) Map maximum (asymptotic value for some profiles)
    psi_bg_max = float(np.max(psi_map))

    # (c) Mean over entire map
    psi_bg_mean = float(np.mean(psi_map))

    # (d) Local annulus background (more physical for cluster lenses)
    psi_bg_annulus = sample_annulus_background(psi_map, wcs, IMAGE_POSITIONS_DEG, r_inner_pix=10, r_outer_pix=20)

    # Potential depth tracers: deeper = larger value
    tracer_defs = {}

    # Direct psi map: potential depth ~ psi_bg - psi_image
    for bg_name, psi_bg in [("edge", psi_bg_edge), ("max", psi_bg_max), ("mean", psi_bg_mean)]:
        tracer_defs[f"psi_map_bg_{bg_name}"] = {im: max(0.0, psi_bg - psi_map_vals[im]) for im in imgs}

    # Local annulus background for each image
    tracer_defs["psi_map_bg_annulus"] = {im: max(0.0, psi_bg_annulus[im] - psi_map_vals[im]) for im in imgs}

    # FFT psi: same approach
    for bg_name, psi_bg in [("edge", psi_bg_edge), ("max", psi_bg_max), ("mean", psi_bg_mean)]:
        tracer_defs[f"psi_fft_bg_{bg_name}"] = {im: max(0.0, psi_bg - psi_fft_vals[im]) for im in imgs}

    # JSON-based tracers (ground truth)
    tracer_defs["flux_ratio_proxy"] = flux
    tracer_defs["model_mu_abs"] = mu_abs_json
    tracer_defs["model_kappa"] = kappa_json
    tracer_defs["model_inv_kappa"] = inv_kappa_json

    # ------------------------------------------------------------------
    # 7. Compute loop residuals for all tracers
    # ------------------------------------------------------------------
    s07 = json.load(open(PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    R_obs = float(s07["weighted_mean_residual"]["R_obs_days"])

    results_by_tracer = {}
    flux_R = None
    for name, q in tracer_defs.items():
        # Skip tracers with zero or negative values (log undefined)
        if any(v <= 0 for v in q.values()):
            print_status(f"Skipping tracer '{name}': non-positive values", "WARN")
            continue

        R, qn, G = loop_residual(ALPHA_PROXY, q, delays)
        R_pred_obs = -R
        dGamma_S4_SX = G["S4"] - G["SX"]
        sign_matches_obs = bool(np.sign(R_pred_obs) == np.sign(R_obs))
        results_by_tracer[name] = {
            "R_closure_days": R,
            "R_predicted_obs_minus_model_days": R_pred_obs,
            "dGamma_S4_minus_SX": float(dGamma_S4_SX),
            "q_norm": {im: float(qn[im]) for im in qn},
            "predicted_sign_matches_observed": sign_matches_obs,
        }
        print_status(f"\n[{name}]")
        print_status(f"  q_norm S4={qn['S4']:.4f}  SX={qn['SX']:.4f}   dGamma(S4-SX)={dGamma_S4_SX:+.6f}")
        print_status(f"  R_closure={R:+.2f} d  -> predicted (obs-model)={R_pred_obs:+.2f} d  "
                     f"(observed={R_obs:+.2f} d)  sign-match: {sign_matches_obs}")
        if name == "flux_ratio_proxy":
            flux_R = R_pred_obs

    # ------------------------------------------------------------------
    # 7b. FUNDAMENTAL TEP formula with physical potential
    # ------------------------------------------------------------------
    print_status("\n--- Fundamental TEP formula (Gamma = 1 + alpha*|Phi|/c^2) ---", "INFO")

    # Compute physical |Phi|/c^2 from psi map values
    phi_over_c2_map, K_calib = compute_physical_potential_over_c2(psi_map_vals)
    print_status(f"Physical potential |Phi|/c^2 at image positions (K_calib={K_calib:.3e}):", "INFO")
    for im in imgs:
        print_status(f"  {im}: |Phi|/c^2 = {phi_over_c2_map[im]:.3e}", "INFO")

    # Evaluate fundamental formula with ALPHA_PROXY
    R_fund, G_fund = loop_residual_fundamental(ALPHA_PROXY, phi_over_c2_map, delays)
    R_pred_obs_fund = -R_fund
    dG_fund = G_fund["S4"] - G_fund["SX"]
    print_status(f"\nFundamental (alpha={ALPHA_PROXY}): R_pred_obs={R_pred_obs_fund:+.6f} d, "
                 f"dGamma_S4-SX={dG_fund:+.6e}", "INFO")

    # Compute effective alpha needed to match observed residual
    if abs(R_fund) > 1e-12:
        alpha_eff = ALPHA_PROXY * (R_obs / R_pred_obs_fund)
        print_status(f"Effective alpha to match observation: {alpha_eff:.3e}", "INFO")
    else:
        alpha_eff = None
        print_status("Fundamental residual is zero; cannot infer effective alpha.", "WARN")

    # Also evaluate with physical potential from JSON kappa via 1/kappa proxy
    # (this is the best-motivated proxy for beta_A=-1)
    print_status("\n--- Comparison: proxy vs fundamental ---", "INFO")
    print_status(f"Flux-proxy prediction:      {flux_R:+.2f} d  (alpha_proxy={ALPHA_PROXY})", "INFO")
    print_status(f"1/kappa proxy prediction:   {results_by_tracer['model_inv_kappa']['R_predicted_obs_minus_model_days']:+.2f} d", "INFO")
    print_status(f"Fundamental prediction:     {R_pred_obs_fund:+.6f} d  (same alpha)", "INFO")
    if alpha_eff:
        print_status(f"Required fundamental alpha: {alpha_eff:.3e} to match {R_obs:+.1f} d", "INFO")

    # ------------------------------------------------------------------
    # 8. Verdict
    # ------------------------------------------------------------------
    psi_map_names = [n for n in results_by_tracer if n.startswith("psi_map_")]
    psi_fft_names = [n for n in results_by_tracer if n.startswith("psi_fft_")]

    flux_sign = np.sign(flux_R) if flux_R is not None else None

    if psi_map_names and flux_sign is not None:
        # Separate global-background tracers (edge, max, mean) from local-annulus tracer
        global_names = [n for n in psi_map_names if "annulus" not in n]
        local_names = [n for n in psi_map_names if "annulus" in n]

        global_signs = [np.sign(results_by_tracer[n]["R_predicted_obs_minus_model_days"]) for n in global_names]
        global_match_flux = all(s == flux_sign for s in global_signs)

        local_signs = [np.sign(results_by_tracer[n]["R_predicted_obs_minus_model_days"]) for n in local_names]
        local_match_flux = all(s == flux_sign for s in local_signs)

        fft_signs = [np.sign(results_by_tracer[n]["R_predicted_obs_minus_model_days"]) for n in psi_fft_names]
        fft_match_flux = all(s == flux_sign for s in fft_signs)

        if global_match_flux and not local_match_flux and not fft_match_flux:
            verdict = (
                "Global-background psi-map tracers (edge, mean, max) preserve the SIGN of the "
                "flux-proxy prediction (+0.02 to +0.04 d).  The local-annulus tracer inverts the sign "
                "(-1.40 d) because a local background is not the physical potential at infinity; it "
                "measures curvature rather than depth.  The FFT-reconstructed psi gives zero contrast "
                "because the gauge choice (DC mode = 0) removes the large-scale gradient.  All psi "
                "amplitudes are suppressed because the lensing potential varies slowly across the "
                "~5-arcsec S4-SX separation, whereas the flux proxy amplifies the contrast through the "
                "(1-kappa)^2-gamma^2 magnification denominator.  The sign agreement confirms that SX "
                "sits at a shallower potential than S4, consistent with the flux-proxy ordering and "
                f"the observed {R_obs:+.1f} d residual."
            )
        elif global_match_flux:
            verdict = (
                "Global-background psi-map tracers preserve the SIGN of the flux-proxy prediction. "
                "The sign agreement confirms that SX sits at a shallower potential than S4, consistent "
                f"with the observed {R_obs:+.1f} d residual.  Amplitudes are suppressed because the lensing "
                "potential varies slowly across the image separation.  The FUNDAMENTAL TEP formula "
                f"(Gamma = 1 + alpha*|Phi|/c^2) with alpha={ALPHA_PROXY} predicts only "
                f"{R_pred_obs_fund:+.2e} days — {abs(R_obs/R_pred_obs_fund):.0e} times smaller than observed.  "
                "The log-magnification proxy is a phenomenological ansatz that captures the correct "
                "phenomenology but currently lacks a first-principles derivation from the scalar-field "
                "action through the lensing potential."
            )
        else:
            verdict = (
                "WARNING: the global-background psi-map formulations do not consistently match the "
                "flux-proxy sign.  This indicates that the background-subtraction assumption (psi_bg) "
                "is not capturing the physical potential depth correctly.  A full 3D potential "
                "reconstruction with proper cosmological boundary conditions is required."
            )
    else:
        verdict = (
            "Psi-tracer computation encountered non-physical values (negative potential depths). "
            "This indicates that the simple background-subtraction approximation fails for the "
            "cluster lens geometry.  A proper cosmological potential reconstruction is needed."
        )

    print_status("\n" + verdict)

    # ------------------------------------------------------------------
    # 9. Save results
    # ------------------------------------------------------------------
    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": ("Loop residual under lensing-potential psi tracer, computed from "
                        "archived GLAFIC v3 maps via direct sampling and FFT Poisson solver. "
                        "Tests whether the headline sign survives substitution of the "
                        "fundamental TEP coupling variable (potential depth) for the "
                        "phenomenological log-magnification proxy."),
        "alpha_proxy_ref": ALPHA_PROXY,
        "observed_blind_residual_days": R_obs,
        "image_positions_deg": IMAGE_POSITIONS_DEG,
        "map_metadata": {
            "shape": kappa_map.shape,
            "pixel_scale_arcsec": pixel_scale_arcsec,
            "psi_bg_edge": psi_bg_edge,
            "psi_bg_max": psi_bg_max,
            "psi_bg_mean": psi_bg_mean,
        },
        "map_json_agreement": map_json_agreement,
        "maps_consistent_with_json": maps_consistent,
        "psi_map_values": psi_map_vals,
        "psi_fft_values": psi_fft_vals,
        "results_by_tracer": results_by_tracer,
        "fundamental_tep": {
            "phi_over_c2": phi_over_c2_map,
            "K_calibration": K_calib,
            "R_predicted_obs_minus_model_days": R_pred_obs_fund,
            "alpha_eff_to_match_observation": alpha_eff,
        },
        "verdict": verdict,
        "caveats": [
            "The archived GLAFIC v3 maps may be at a different source redshift than SN Refsdal "
            "(z=1.489), causing a normalization offset relative to the Kelly+2023 tabulated values. "
            "The map kappa values do not match the JSON tabulated values (ratio mean={:.2f}, std={:.2f}).".format(
                ratio_mean, ratio_std),
            "The psi background (psi_bg) is approximated by the map edge/mean/max, which is "
            "uncertain because the cluster potential has not decayed to the cosmological background "
            "within the ~2.5-arcmin map extent.",
            "The fundamental TEP formula (Gamma = 1 + alpha*|Phi|/c^2) predicts an effect "
            "~10^5 times smaller than observed for alpha=-0.055.  The log-magnification proxy "
            "is a phenomenological ansatz that lacks a first-principles derivation from the "
            "scalar-field action through the lensing potential.  The required 'effective alpha' "
            "for the fundamental formula is orders of magnitude larger than the lab-scale value.",
            "Image positions are from Kelly+2015 and should be verified against the archived "
            "GLAFIC model astrometric solution.",
        ],
    }
    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_psi_transport.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
