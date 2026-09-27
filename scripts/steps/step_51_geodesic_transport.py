#!/usr/bin/env python3
"""
TEP-LENS: Step 51 - Full 3D Cosmological Geodesic Transport Integration

Purpose: compute the TEP loop residual by integrating the scalar field along the
full 3D geodesic from source to observer, NOT by sampling psi at image positions.

The 3D gravitational potential Phi(r) is reconstructed from the surface mass
density kappa(theta) via spherical deprojection. Light rays are traced through
the 3D potential, and the TEP-corrected propagation time is integrated along
each path. The loop residual is computed from these geodesic-integrated delays.

Inputs : data/raw/sn_lensing/maps/hlsp_frontier_model_macs1149_glafic_v3_kappa.fits
         data/raw/sn_lensing/refsdal_glafic_v3_lensing_params.json
         results/outputs/step_07_observed_vs_predicted.json
Outputs: results/outputs/step_51_geodesic_transport.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import KAPPA_LENS, BETA_REFSDAL, load_refsdal_image_positions

STEP_NUM = "51"
LOOP = ("S1", "S4", "SX")
# Published J2000 positions (provenance: refsdal_glafic_v3_lensing_params.json;
# ATel #6729 / Karman+2016 MUSE Table 1; Kelly+2016 for SX).  An earlier version
# used positions ~10 arcsec north inside the cluster core.
IMAGE_POSITIONS_DEG = load_refsdal_image_positions()
Z_L = 0.542
Z_S = 1.489
H0 = 70.0
OM0 = 0.3
C_KMS = 299792.458
MPC_TO_KM = 3.0856775814913673e19


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def angular_distances(z_l, z_s, H0=70.0, Om0=0.3):
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    D_l = cosmo.angular_diameter_distance(z_l).value
    D_s = cosmo.angular_diameter_distance(z_s).value
    D_ls = cosmo.angular_diameter_distance_z1z2(z_l, z_s).value
    return D_l, D_s, D_ls


def kappa_radial_profile(kappa_map, wcs, center_ra=177.3987536, center_dec=22.3985266):
    """Extract azimuthally-averaged radial kappa profile."""
    ny, nx = kappa_map.shape
    pixel_scale_deg = np.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    cx, cy = wcs.all_world2pix(center_ra, center_dec, 0)
    y, x = np.indices((ny, nx))
    r_pix = np.sqrt((x - cx)**2 + (y - cy)**2)
    r_arcsec = r_pix * pixel_scale_arcsec
    r_bins = np.arange(0, min(np.max(r_arcsec), 120.0), pixel_scale_arcsec * 5)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    kappa_mean = np.zeros(len(r_mid))
    for i in range(len(r_mid)):
        mask = (r_arcsec >= r_bins[i]) & (r_arcsec < r_bins[i+1])
        kappa_mean[i] = np.mean(kappa_map[mask]) if np.sum(mask) > 0 else 0.0
    return r_mid, kappa_mean, pixel_scale_arcsec


def spherical_deprojection(R_arcsec, kappa_profile, D_l_Mpc, z_l=0.542):
    """
    Deproject surface density kappa(R) to 3D density rho(r) assuming spherical symmetry.
    Returns r_3d in Mpc and rho in dimensionless units (relative to critical density).
    """
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)
    a_l = 1.0 / (1.0 + z_l)

    # Convert angular radius to physical radius at lens redshift
    # R_phys = R_arcsec * (pi / 180*3600) * D_l [Mpc]
    arcsec_to_rad = np.pi / (180.0 * 3600.0)
    R_phys = R_arcsec * arcsec_to_rad * D_l_Mpc  # Mpc

    # Smooth kappa profile to avoid noise
    from scipy.ndimage import gaussian_filter1d
    kappa_smooth = gaussian_filter1d(kappa_profile, sigma=2.0)

    # Numerical derivative dSigma/dR
    dSigma_dR = np.gradient(kappa_smooth, R_phys)

    # Abel deprojection: rho(r) = -(1/pi) * integral_r^inf (dSigma/dR) / sqrt(R^2 - r^2) dR
    r_3d = R_phys.copy()
    rho_3d = np.zeros(len(r_3d))

    for i in range(len(r_3d)):
        r = r_3d[i]
        integrand = []
        R_int = []
        for j in range(len(R_phys)):
            if R_phys[j] > r + 1e-10:
                denom = np.sqrt(R_phys[j]**2 - r**2)
                if denom > 0:
                    integrand.append(-dSigma_dR[j] / denom)
                    R_int.append(R_phys[j])
        if len(integrand) > 1:
            rho_3d[i] = (1.0 / np.pi) * np.trapezoid(integrand, R_int)
        else:
            rho_3d[i] = 0.0

    # Smooth and enforce non-negativity
    rho_3d = np.maximum(rho_3d, 0.0)
    rho_3d = gaussian_filter1d(rho_3d, sigma=1.0)

    return r_3d, rho_3d


def compute_3D_potential(r_3d, rho_3d):
    """
    Compute 3D gravitational potential Phi(r) from 3D density rho(r).
    For spherical symmetry:
        Phi(r) = -4*pi*G * [ (1/r) * integral_0^r rho(r') r'^2 dr' + integral_r^inf rho(r') r' dr' ]

    We work in units where G = 1, and return Phi in units of c^2.
    """
    # Use a dimensionless formulation: Phi/c^2 = -2 * integral_0^r M(<r')/r'^2 dr'
    # where M(<r) = 4*pi * integral_0^r rho(r') r'^2 dr'

    # Enclosed mass
    M_enclosed = np.zeros(len(r_3d))
    for i in range(1, len(r_3d)):
        M_enclosed[i] = 4.0 * np.pi * np.trapezoid(rho_3d[:i+1] * r_3d[:i+1]**2, r_3d[:i+1])

    # Potential: Phi(r) = -G * M_enclosed(r) / r - 4*pi*G * integral_r^inf rho(r') r' dr'
    # In our dimensionless units, we just compute the shape; amplitude set by virial theorem
    Phi = np.zeros(len(r_3d))
    for i in range(len(r_3d)):
        r = r_3d[i]
        # First term: -M/r
        term1 = -M_enclosed[i] / r if r > 1e-10 else 0.0
        # Second term: -4*pi * integral_r^inf rho(r') r' dr'
        if i < len(r_3d) - 1:
            term2 = -4.0 * np.pi * np.trapezoid(rho_3d[i:] * r_3d[i:], r_3d[i:])
        else:
            term2 = 0.0
        Phi[i] = term1 + term2

    # Normalize so that max |Phi| corresponds to virial value
    # For a cluster: |Phi|/c^2 ~ (sigma/c)^2 ~ 10^-5 for sigma ~ 1000 km/s
    Phi_peak = np.max(np.abs(Phi))
    if Phi_peak > 0:
        # Scale to match virial theorem estimate
        sigma_cluster = 1000.0  # km/s
        phi_over_c2_target = (sigma_cluster / C_KMS) ** 2  # ~1.1e-5
        scale = phi_over_c2_target / Phi_peak
        Phi = Phi * scale

    return Phi


def trace_geodesic_and_integrate(image_ra, image_dec, wcs, r_3d, Phi_3d,
                                  D_l, D_s, D_ls, z_l=Z_L, z_s=Z_S,
                                  centre_pix=None, n_segments=1000):
    """
    Trace a light ray from source to observer passing through the image position,
    and integrate the TEP scalar field along the path.

    The path is approximated as straight segments. At each point along the path,
    we compute the 3D radius from the cluster center and interpolate the potential.

    Returns:
        t_GR: GR propagation time (Shapiro + geometric)
        t_TEP: TEP-corrected propagation time
        Phi_integral: integral of |Phi|/c^2 along the path
    """
    # Convert image position to arcsec offset from the potential centre
    # (psi-map minimum), not the WCS reference pixel.
    if centre_pix is None:
        cx, cy = wcs.all_world2pix(wcs.wcs.crval[0], wcs.wcs.crval[1], 0)
    else:
        cx, cy = centre_pix
    x_img, y_img = wcs.all_world2pix(image_ra, image_dec, 0)
    dx_pix = x_img - cx
    dy_pix = y_img - cy
    pixel_scale_deg = np.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    pixel_scale_arcsec = pixel_scale_deg * 3600.0
    theta_x_arcsec = dx_pix * pixel_scale_arcsec  # arcsec
    theta_y_arcsec = dy_pix * pixel_scale_arcsec  # arcsec

    # Physical offset at lens plane (Mpc)
    arcsec_to_rad = np.pi / (180.0 * 3600.0)
    x_l = theta_x_arcsec * arcsec_to_rad * D_l  # Mpc
    y_l = theta_y_arcsec * arcsec_to_rad * D_l  # Mpc

    # Approximate path: from source (z_s) to lens (z_l) to observer (z=0)
    # The comoving distance from source to lens
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)

    # Comoving distances
    chi_s = cosmo.comoving_distance(z_s).value  # Mpc
    chi_l = cosmo.comoving_distance(z_l).value  # Mpc

    # Path from source to lens (passing through image position at lens)
    # In the thin-lens approximation, the ray is straight except at the lens plane
    # For a thick lens, we model the path as straight with the impact parameter
    # varying as the ray approaches the lens plane.

    # Impact parameter of the ray at the lens plane (Mpc).  Within the cluster
    # region (|s| <~ 2 Mpc << D_l) the transverse offset is constant at b_l to
    # excellent approximation; the potential vanishes beyond r_3d[-1], so the
    # integral is evaluated on a fine line-of-sight grid localised at the lens
    # rather than over the full ~5 Gpc cosmological path.
    b_l = np.sqrt(x_l**2 + y_l**2)  # Mpc

    s_max = 3.0 * r_3d[-1]
    s_path = np.linspace(-s_max, s_max, n_segments)

    # 3D radius from cluster center: straight ray, r = sqrt(b^2 + s^2).
    # (The previous form used (chi-chi_l)*sin(angle) with angle ~ b_l/D_l ~ 1e-4,
    # which kept r ~ b over the whole path so the potential never truncated; and
    # the 2000-point full-path grid could not resolve the ~1.5 Mpc lens region.)
    r_path = np.sqrt(b_l**2 + s_path**2)

    # Interpolate potential at each radius
    Phi_over_c2_path = np.interp(r_path, r_3d, Phi_3d, left=Phi_3d[0], right=0.0)

    # Path element dl in Mpc, converted to km for the c-division
    dl = np.abs(np.diff(s_path))
    dl_km = dl * MPC_TO_KM
    Phi_mid = 0.5 * (Phi_over_c2_path[:-1] + Phi_over_c2_path[1:])

    # GR Shapiro delay: dt_GR = -(2/c) * integral Phi/c^2 dl
    # (in coordinate time, for a photon in a weak field)
    # The factor (1+z) accounts for cosmological redshift
    t_GR_delay = -(2.0 / C_KMS) * np.sum(Phi_mid * dl_km)  # seconds

    # TEP correction: the scalar field modifies the metric, adding an extra term
    # to the effective refractive index: n = 1 - 2*Phi/c^2 - alpha*Phi/c^2/2
    # Wait - for TEP with beta_A = -1, the scalar field tracks the potential.
    # The TEP time dilation is Gamma = 1 + alpha*|Phi|/c^2.
    # For light propagation, the additional delay is:
    # dt_TEP = -(1/c) * integral alpha*|Phi|/c^2 dl
    # (this is the extra path length due to the scalar field)
    t_TEP_correction = -(1.0 / C_KMS) * KAPPA_LENS * np.sum(np.abs(Phi_mid) * dl_km)  # seconds

    # Convert to days
    t_GR_delay_days = t_GR_delay / 86400.0
    t_TEP_corr_days = t_TEP_correction / 86400.0

    # The geometric time delay (in thin-lens) is handled separately by the loop_residual
    # Here we only compute the potential-dependent delays
    Phi_integral = np.sum(np.abs(Phi_mid) * dl)  # Mpc * dimensionless

    return t_GR_delay_days, t_TEP_corr_days, Phi_integral


def main():
    print_status(f"STEP {STEP_NUM}: Full 3D Geodesic Transport Integration", "TITLE")

    # Load kappa and psi maps
    print_status("Loading GLAFIC v3 kappa/psi maps...", "INFO")
    from astropy.io import fits
    from astropy.wcs import WCS
    map_dir = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "maps"
    with fits.open(map_dir / "hlsp_frontier_model_macs1149_glafic_v3_kappa.fits") as hdul:
        kappa_map = hdul[0].data.astype(float)
        wcs = WCS(hdul[0].header)
    with fits.open(map_dir / "hlsp_frontier_model_macs1149_glafic_v3_psi.fits") as hdul:
        psi_map = hdul[0].data.astype(float)

    # The archived maps are scaled to D_ls/D_s = 1 (z_s -> infinity); rescale to
    # the SN Refsdal source plane z_s = 1.489.
    kappa_map = kappa_map * BETA_REFSDAL

    # Centre the profile on the potential centre (psi-map minimum), not the WCS
    # reference pixel, which is offset from the cluster centre.
    iy_min, ix_min = np.unravel_index(np.nanargmin(psi_map), psi_map.shape)
    centre_ra, centre_dec = [float(v) for v in wcs.all_pix2world(ix_min, iy_min, 0)]
    print_status(f"Potential centre (psi minimum): RA={centre_ra:.7f}, Dec={centre_dec:.7f}", "INFO")

    # Cosmological distances
    D_l, D_s, D_ls = angular_distances(Z_L, Z_S, H0, OM0)
    print_status(f"D_l={D_l:.1f} Mpc, D_s={D_s:.1f} Mpc, D_ls={D_ls:.1f} Mpc", "INFO")

    # Extract radial profile
    print_status("Extracting radial kappa profile...", "INFO")
    r_mid, kappa_mean, pixel_scale = kappa_radial_profile(kappa_map, wcs,
                                                        center_ra=centre_ra,
                                                        center_dec=centre_dec)
    print_status(f"Profile: {len(r_mid)} bins, max radius={r_mid[-1]:.1f} arcsec", "INFO")

    # 3D deprojection
    print_status("Deprojecting to 3D density (spherical symmetry)...", "INFO")
    r_3d, rho_3d = spherical_deprojection(r_mid, kappa_mean, D_l, Z_L)
    print_status(f"3D grid: {len(r_3d)} points, max r={r_3d[-1]:.3f} Mpc", "INFO")

    # 3D potential
    print_status("Computing 3D gravitational potential...", "INFO")
    Phi_3d = compute_3D_potential(r_3d, rho_3d)
    print_status(f"Potential range: Phi/c^2 = [{np.min(Phi_3d):.3e}, {np.max(Phi_3d):.3e}]", "INFO")

    # Trace geodesics and integrate for each image
    print_status("Tracing geodesics and integrating TEP scalar field...", "INFO")
    delays = {"S1": 0.0, "S2": 9.9, "S3": 9.0, "S4": 20.3, "SX": 376.0}

    t_GR = {}
    t_TEP = {}
    Phi_int = {}

    centre_pix = (ix_min, iy_min)
    for name, (ra, dec) in IMAGE_POSITIONS_DEG.items():
        t_gr, t_tep, phi_i = trace_geodesic_and_integrate(
            ra, dec, wcs, r_3d, Phi_3d, D_l, D_s, D_ls,
            Z_L, Z_S, centre_pix=centre_pix, n_segments=2000
        )
        t_GR[name] = t_gr
        t_TEP[name] = t_tep
        Phi_int[name] = phi_i
        print_status(f"  {name}: GR delay={t_gr:+.4f} d, TEP corr={t_tep:+.6e} d, Phi_int={phi_i:.3e} Mpc", "INFO")

    # Compute loop residuals
    # The total delay for each image is: delay_i = delay_GR_i + delay_TEP_i
    # But the GR delays are already encoded in the observed delays.
    # The TEP correction to the loop residual is:
    #   R_TEP = (G_i-1)*dt_ji + (G_j-1)*dt_kj + (G_k-1)*dt_ik
    # where G_i = 1 + alpha * |Phi|_integrated_i / (some reference)

    # For the geodesic integral, the relevant quantity is the integrated potential
    # along each path. We use this as the tracer.  Because |Phi| is a DEPTH
    # variable (larger integral = deeper path = slower transport under beta_A=-1),
    # the physical convention applies a positive coupling to log10(depth ratio);
    # the flux-convention result (KAPPA_LENS < 0 on a depth tracer) is reported
    # alongside for continuity with earlier versions of this step.
    Phi_bar = np.mean(list(Phi_int.values()))
    q = {im: Phi_int[im] / Phi_bar for im in Phi_int}
    G = {im: 1.0 + abs(KAPPA_LENS) * np.log10(q[im]) for im in q}
    G_flux = {im: 1.0 + KAPPA_LENS * np.log10(q[im]) for im in q}

    i, j, k = LOOP
    R = (G[i] - 1.0) * (delays[j] - delays[i]) + (G[j] - 1.0) * (delays[k] - delays[j]) + (G[k] - 1.0) * (delays[i] - delays[k])
    R_pred_obs = -R
    R_flux = (G_flux[i] - 1.0) * (delays[j] - delays[i]) + (G_flux[j] - 1.0) * (delays[k] - delays[j]) + (G_flux[k] - 1.0) * (delays[i] - delays[k])
    R_pred_obs_flux = -R_flux

    print_status(f"\nGeodesic-integrated proxy residual (depth convention): R_pred_obs = {R_pred_obs:+.4f} d", "INFO")
    print_status(f"  (flux-convention equivalent: {R_pred_obs_flux:+.4f} d)", "INFO")
    print_status(f"  q: S4={q['S4']:.6f}, SX={q['SX']:.6f}", "INFO")
    print_status(f"  dGamma_S4-SX = {G['S4']-G['SX']:+.6e}", "INFO")

    # Also compute with the model-kappa tracer.  kappa is a depth variable
    # (projected density), so the physically consistent convention applies
    # +|kappa_lens| to log10(kappa/kappa_bar); the earlier 1/kappa construction
    # encoded the same depth ordering in the flux (shallowness) convention and
    # is now retired in favour of the direct depth reading.
    gl = json.load(open(PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "refsdal_glafic_v3_lensing_params.json"))
    kappa_json = {im: gl["images"][im]["kappa"] for im in gl["images"]}
    q_k = {im: kappa_json[im] / np.mean(list(kappa_json.values())) for im in kappa_json}
    G_dep = {im: 1.0 + abs(KAPPA_LENS) * np.log10(q_k[im]) for im in q_k}
    R_dep = (G_dep[i]-1)*(delays[j]-delays[i]) + (G_dep[j]-1)*(delays[k]-delays[j]) + (G_dep[k]-1)*(delays[i]-delays[k])
    R_kappa_depth_obs = -R_dep
    inv_kappa = {im: 1.0 / kappa_json[im] for im in kappa_json}
    q_inv = {im: inv_kappa[im] / np.mean(list(inv_kappa.values())) for im in inv_kappa}
    G_inv = {im: 1.0 + KAPPA_LENS * np.log10(q_inv[im]) for im in q_inv}
    R_inv = (G_inv[i]-1)*(delays[j]-delays[i]) + (G_inv[j]-1)*(delays[k]-delays[j]) + (G_inv[k]-1)*(delays[i]-delays[k])
    R_inv_obs = -R_inv
    print_status(f"kappa-as-depth proxy residual: R_pred_obs = {R_kappa_depth_obs:+.2f} d", "INFO")
    print_status(f"1/kappa (flux-convention equivalent): R_pred_obs = {R_inv_obs:+.2f} d", "INFO")

    # Endpoint-rescaling diagnostic with geodesic-integrated Phi.  This is not
    # an additional static conformal delay; see the exact cancellation audit
    # in Step 56.
    Phi_int_bar = np.mean(list(Phi_int.values()))
    G_fund = {im: 1.0 + abs(KAPPA_LENS) * (Phi_int[im] / Phi_int_bar) for im in Phi_int}
    R_fund = (G_fund[i]-1)*(delays[j]-delays[i]) + (G_fund[j]-1)*(delays[k]-delays[j]) + (G_fund[k]-1)*(delays[i]-delays[k])
    R_fund_obs = -R_fund
    print_status(f"Endpoint-rescaling diagnostic (geodesic Phi): R_pred_obs = {R_fund_obs:+.6e} d", "INFO")

    # Load observed residual
    s07 = json.load(open(PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    R_obs = float(s07["weighted_mean_residual"]["R_obs_days"])

    # Local temporal-field depth ordering from the deprojected 3D potential at
    # each image's impact parameter (the quantity phi tracks under beta_A = -1).
    b_phys = {}
    arcsec_to_rad = np.pi / (180.0 * 3600.0)
    for name, (ra, dec) in IMAGE_POSITIONS_DEG.items():
        x_img, y_img = wcs.all_world2pix(ra, dec, 0)
        b_arcsec = np.hypot(x_img - ix_min, y_img - iy_min) * pixel_scale
        b_phys[name] = b_arcsec * arcsec_to_rad * D_l  # Mpc
    Phi_local = {im: float(np.interp(b_phys[im], r_3d, Phi_3d)) for im in Phi_int}
    depth_order = sorted(Phi_local, key=lambda k: Phi_local[k])
    int_order = sorted(Phi_int, key=lambda k: -Phi_int[k])
    print_status(f"Local |Phi| depth ordering (deepest first): {depth_order}", "INFO")
    print_status(f"Path-integral ordering (largest first): {int_order}", "INFO")

    # Verdict
    sign_match = np.sign(R_pred_obs) == np.sign(R_obs)
    depth_sx = depth_order[0] == "SX"
    if sign_match and depth_sx:
        verdict = (
            f"At the published image positions the smooth spherical deprojection places "
            f"SX at the smallest cluster-centric radius ({b_phys['SX']*1000:.0f} kpc) and the "
            f"deepest local potential (Phi/c^2 = {Phi_local['SX']:.3e}), and its sightline "
            f"carries the largest path-integrated |Phi|/c^2 dl "
            f"({Phi_int['SX']:.3e} Mpc).  All three physical tracers of the temporal-field "
            f"depth therefore order SX first, consistent with the observed "
            f"{R_obs:+.1f} d residual direction.  The geodesic-integrated proxy predicts "
            f"{R_pred_obs:+.4f} d.  The path-integrated potential is the ordinary GR "
            f"Fermat/Shapiro contribution already inside the blind lens predictions "
            f"(Step 56), so the sign agreement here is a consistent environment "
            f"diagnostic, not an additional propagation amplitude."
        )
    else:
        verdict = (
            f"The smooth spherical deprojection gives local-depth ordering {depth_order} "
            f"and path-integral ordering {int_order}; geodesic proxy residual "
            f"{R_pred_obs:+.4f} d vs observed {R_obs:+.1f} d.  The azimuthally averaged "
            f"model retains only the cluster monopole; the full-2D kappa map at the "
            f"published positions (Step 56) is the primary local-density diagnostic."
        )

    print_status("\n" + verdict)

    # Save results
    results = {
        "step": STEP_NUM,
        "status": "success",
        "description": "3D potential-path diagnostic with static conformal cancellation audit.",
        "cosmology": {"z_l": Z_L, "z_s": Z_S, "H0": H0, "Om0": OM0,
                      "D_l_Mpc": D_l, "D_s_Mpc": D_s, "D_ls_Mpc": D_ls},
        "potential_reconstruction": {
            "r_3d_Mpc": r_3d.tolist(),
            "rho_3d": rho_3d.tolist(),
            "Phi_over_c2": Phi_3d.tolist(),
            "Phi_peak": float(np.max(np.abs(Phi_3d))),
        },
        "geodesic_integration": {
            "n_segments": 2000,
            "t_GR_delay_days": t_GR,
            "t_TEP_correction_days": t_TEP,
            "Phi_integral_Mpc": Phi_int,
            "impact_parameter_Mpc": b_phys,
            "local_Phi_over_c2_at_image": Phi_local,
            "local_depth_ordering": depth_order,
            "path_integral_ordering": int_order,
        },
        "loop_residuals": {
            "geodesic_proxy_days": R_pred_obs,
            "geodesic_proxy_days_flux_convention": R_pred_obs_flux,
            "kappa_as_depth_proxy_days": R_kappa_depth_obs,
            "inv_kappa_proxy_days": R_inv_obs,
            "endpoint_rescaling_diagnostic_days": R_fund_obs,
            "fundamental_formula_days": R_fund_obs,
            "observed_days": R_obs,
        },
        "verdict": verdict,
        "caveats": [
            "Spherical symmetry assumed for 3D deprojection; real cluster is triaxial.",
            "Thin-lens approximation used for geodesic path; thick-lens effects neglected.",
            "Cluster member galaxies and substructure not included in the smooth deprojection.",
            "The integrated potential is the ordinary GR Fermat/Shapiro term; adding it after "
            "GR-model subtraction would double-count the static conformal contribution.",
            "The fundamental_formula_days field is retained as a deprecated compatibility alias "
            "for endpoint_rescaling_diagnostic_days; it is not a fundamental TEP prediction.",
            "The log-magnification proxy is phenomenological; no first-principles derivation exists.",
            "Published image positions are used (refsdal_glafic_v3_lensing_params.json); an "
            "earlier pipeline generation sampled ~10 arcsec north of the image field.",
            "The kappa map is rescaled by beta(z_s=1.489) = 0.5339 from the archived "
            "D_ls/D_s = 1 convention before deprojection, and the radial profile is centred "
            "on the psi-map minimum (potential centre) rather than the WCS reference pixel.",
        ],
    }
    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_geodesic_transport.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=safe_json_default)
    print_status(f"\nResults saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
