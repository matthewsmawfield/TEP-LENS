#!/usr/bin/env python3
"""
TEP-LENS: Step 36b - SN 2025wny post-measurement update (measured A-D baseline)

Step 36 is the prospective forecast and is left unchanged. It assumed the post-hoc
Witt-Wynne geometric baseline (~175 d, modelled as U(140, 210)). The SN 2025wny
time delays are now published (August 2026):
  - Townsend et al. 2026, arXiv:2608.28427 (photometric) and Johansson et al. 2026,
    arXiv:2608.28416 (spectroscopic)
  - Lens modelling and cosmography (Paper V), arXiv:2608.28430, Table 6:
      dt_AB = -10.4 +/- 1.7 d, dt_AC = 0.8 +/- 2.1 d, dt_AD = -65.7 +/- 3.5 d,
      H0 = 66.7 (+7.6/-6.3) km/s/Mpc (fiducial isothermal, filter-marginalized)
This step re-evaluates the same proxy-model Monte Carlo with the measured A-D baseline
(all other priors and the Step 07 calibrated coupling unchanged), and compares the
predicted residual with the current lens-model precision on that baseline.
Output: results/outputs/step_36b_sn2025wny_measured_baseline_update.json
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "outputs" / "step_36b_sn2025wny_measured_baseline_update.json"

def main():
    rng = np.random.default_rng(42)
    n = 100000
    s07 = json.load(open(ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"))
    bai = s07["bootstrap_kappa_inference"]
    k_mean = float(bai["kappa_mean"]); k_sig = float(bai.get("sigma_kappa_analytical", bai["kappa_std"]))
    mu_b = rng.uniform(15.0, 80.0, n); mu_f = rng.uniform(2.0, 15.0, n)
    mu_m1 = mu_f + rng.uniform(0, 1, n) * (mu_b - mu_f); mu_m2 = mu_f + rng.uniform(0, 1, n) * (mu_b - mu_f)
    mu_mean = (mu_b + mu_m1 + mu_m2 + mu_f) / 4.0
    dt = np.abs(rng.normal(65.7, 3.5, n))            # measured |dt_AD| (Paper V, Table 6)
    kap = rng.normal(k_mean, k_sig, n)
    dgam = kap * (np.log10(mu_b / mu_mean) - np.log10(mu_f / mu_mean))
    R = np.abs(dgam * dt * rng.uniform(0.6, 1.0, n))
    p = {f"p{q}": float(np.percentile(R, q)) for q in (16, 50, 84, 95)}
    # current lens-model precision on the A-D baseline from the H0 posterior width
    frac = 0.5 * (7.6 + 6.3) / 66.7
    sigma_model = frac * 65.7

    # Cosmographic H0 direction & shift under the proxy delay residual:
    # In strong-lensing cosmography H0_inferred ~ 1 / dt_obs; a positive TEP
    # delay residual delta_t_TEP > 0 lengthens the observed baseline and biases
    # a standard GR lens model toward a lower inferred H0 (just as in SN Refsdal:
    # 66.6 -> 69.3 km/s/Mpc). Mörtsell et al. (2026) obtain H0_GR = 66.7 (+7.6/-6.3),
    # sitting below the local Cepheid ladder (73.0) and within 0.21 sigma of TEP's
    # zero-parameter cosmic root H0 = 68.17 km/s/Mpc (Paper 18 Step 21); correcting
    # the 65.7 d baseline by the median +4.37 d proxy residual shifts H0 to 71.5.
    h0_gr = 66.7
    h0_gr_sigma_sym = 0.5 * (7.6 + 6.3)
    h0_tep_root = 68.167
    h0_tep_corrected_p50 = float(h0_gr * 65.7 / (65.7 - p["p50"]))
    h0_tep_corrected_p16 = float(h0_gr * 65.7 / (65.7 - p["p16"]))
    h0_tep_corrected_p84 = float(h0_gr * 65.7 / (65.7 - p["p84"]))

    # Action-level two-branch screening & phantom-mass fraction at the galaxy-scale
    # Einstein ring of SN 2025wny (z_l = 0.375, z_s = 2.011, theta_E ~ 0.85 arcsec):
    G_si = 6.67430e-11
    c_si = 299792458.0
    M_sun = 1.98847e30
    kpc_m = 3.085677581491367e19
    Mpc_m = 1e3 * kpc_m
    H0_si = h0_tep_root * 1e3 / Mpc_m
    g_t_si = c_si * H0_si / 2.0
    k_tb = 16.03
    Om0 = 0.3026

    def _dc(z):
        zz = np.linspace(0.0, z, 10000)
        return (c_si / H0_si) * float(np.trapezoid(1.0 / np.sqrt(Om0 * (1.0 + zz) ** 3 + (1.0 - Om0)), zz))

    Dl = _dc(0.375) / 1.375
    Ds = _dc(2.011) / 3.011
    Dls = (_dc(2.011) - _dc(0.375)) / 3.011
    sigma_crit = (c_si ** 2 / (4.0 * np.pi * G_si)) * Ds / (Dl * Dls)
    theta_E_rad = 0.85 * (np.pi / (180.0 * 3600.0))
    R_E_m = Dl * theta_E_rad
    M_E_msun = float(np.pi * R_E_m ** 2 * sigma_crit / M_sun)

    def _solve_y_tb(q_val):
        y_val = 1.0
        c1 = k_tb / np.sqrt(2.0)
        for _ in range(30):
            f_v = c1 * q_val * y_val ** 2 + (q_val ** 2) * (y_val ** 3) - 1.0
            fp_v = 2.0 * c1 * q_val * y_val + 3.0 * (q_val ** 2) * (y_val ** 2)
            y_val -= f_v / fp_v
        return float(y_val)

    einstein_ring_benchmarks = {}
    for fb_label, fb_val in [("fb_0p25_early_type", 0.25), ("fb_0p35_early_type", 0.35)]:
        sigma_b = fb_val * 0.5 * sigma_crit
        g_b = 2.0 * np.pi * G_si * sigma_b
        q_val = float(g_b / g_t_si)
        y_val = _solve_y_tb(q_val)
        m_tot_over_mb = float(1.0 + 2.0 * y_val)
        f_phantom = float((2.0 * y_val) / (1.0 + 2.0 * y_val))
        einstein_ring_benchmarks[fb_label] = {
            "f_b_enclosed": fb_val,
            "q_gb_over_gt": q_val,
            "y_two_branch": y_val,
            "M_tot_over_M_b_2x": m_tot_over_mb,
            "apparent_dark_matter_fraction_inside_R_E": f_phantom,
        }

    out = {
        "step": "36b",
        "inputs": {"dt_AD_days": [-65.7, 3.5], "dt_AB_days": [-10.4, 1.7], "dt_AC_days": [0.8, 2.1],
                   "dt_AE_days_spec": [3.7, 8.8],
                   "resolved_images_space_data": ["A", "B", "C", "D", "E"],
                   "z_sn_space": 2.0151,
                   "z_deflectors": [0.3755, 0.3766],
                   "theta_E_two_deflectors_arcsec": [1.6, 0.75],
                   "H0_lens_model_isothermal": [66.7, 7.6, 6.3],
                   "H0_lens_model_free_slope": [70.8, 8.2, 6.1],
                   "H0_lens_model": [66.7, 7.6, 6.3], "kappa_calibrated_step07": [k_mean, k_sig],
                   "sources": ["arXiv:2608.28413", "arXiv:2608.28416", "arXiv:2608.28427", "arXiv:2608.28430"]},
        "forecast_step36_premise_baseline_days": [140.0, 210.0],
        "predicted_residual_days_measured_baseline": p,
        "current_lens_model_fractional_precision": frac,
        "current_sigma_model_days_on_AD": sigma_model,
        "required_sigma_model_days": {"3sigma_median": p["p50"] / 3, "5sigma_median": p["p50"] / 5},
        "cosmographic_h0_implication": {
            "H0_GR_inferred_km_s_Mpc": h0_gr,
            "H0_TEP_zero_param_cosmic_root_km_s_Mpc": h0_tep_root,
            "tension_vs_TEP_cosmic_root_sigma": float(abs(h0_gr - h0_tep_root) / h0_gr_sigma_sym),
            "H0_TEP_proxy_corrected_p50_km_s_Mpc": h0_tep_corrected_p50,
            "H0_TEP_proxy_corrected_p16_p84_km_s_Mpc": [h0_tep_corrected_p16, h0_tep_corrected_p84],
        },
        "two_branch_einstein_ring_screening": {
            "theta_E_arcsec": 0.85,
            "R_E_kpc": float(R_E_m / kpc_m),
            "M_E_enclosed_Msun": M_E_msun,
            "benchmarks": einstein_ring_benchmarks,
        },
        "status": ("open: predicted residual is below the current lens-model precision on the "
                   "measured A-D baseline; neither confirmed nor excluded"),
    }
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
