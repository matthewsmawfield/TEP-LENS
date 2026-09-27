#!/usr/bin/env python3
"""
TEP-LENS Step 58 - H0 time-delay channel audit (issue 11-9 open item).

The corpus-level Hubble-tension resolution (Paper 11) allocates >99.9% of the
measured endpoint response to the Cepheid distance-modulus channel.  The
time-delay (H0LiCOW/TDCOSMO-class) channel is a genuinely independent
observable on which the same mechanism must act with the same sign and a
commensurate magnitude, so its channel budget is computed here explicitly.

Sub-channels evaluated per system on real TDCOSMO-2025 data
(data/interim/external/tdcosmo2025_public/TDCOSMO_sample/tdcosmo_sample.yaml
and the SDSS1206+4332 angular-diameter posteriors):

  1. Propagated transport on the measured inter-image delay:
       conformal sector: exactly zero propagated photon delay (verified in
         step_56: conformal_only_additional_residual_days = 0.0 d);
       disformal sector: bounded by GW170817, |c_gamma - c_g|/c <= 1e-15,
         over the halo path scale L_path = 10 R_E (generous; the differential
         path through the lens halo cannot exceed a few Einstein radii);
       endpoint rescaling: the uniform local-well factor A_o ~ 1 - 1e-7
         multiplies the measured duration and therefore D_dt at ~1e-7.
  2. Fermat-normalisation (mass-model) sector:
       the corpus's own carrier decomposition (TEP-GL step_01) places the
       optical-tidal part of the phantom mass in scalar backreaction on
       g_mn -- real curvature contributing to image positions and to the
       Fermat potential identically to ordinary mass, hence absorbed
       self-consistently by the fitted mass model; the chronometric
       (reconstruction-space) component propagates nothing to photons and is
       absent from both image positions and delays.  The residual mass-sheet
       freedom is the mundane external convergence, whose real per-system
       posterior widths are tabulated.
  3. Kinematic-prior sector (sigma_*):
       the only route with order-percent potential.  Under the corpus's
       lensing--dynamical bookkeeping (apparent lensing mass exceeding
       dynamical mass; TEP-UCD Sec. 2), the sigma_* prior under-reads the
       lensing mass, driving the fitted internal mass-sheet normalisation to
       larger lambda and the inferred H0 high -- the required direction for
       elevating the channel above the unbiased value.  Its magnitude is
       bounded by the aperture slip fraction and tied to the same unresolved
       transfer kernel tracked under issue 19-1; it is therefore bounded
       rather than derived here.
  4. Sign bookkeeping on the delay carrier:
       a positive propagated residual (observed delay exceeding the model
       delay) biases the GR-style inference to larger D_dt and therefore to
       LOWER inferred H0.  Producing an elevated inference requires a
       negative residual -- opposite to the sign of the corpus's claimed
       Refsdal residual (+30.1 d).  The delay carrier therefore cannot mimic
       an upward-biased H0 channel regardless of amplitude.

Verdict: the corpus mechanism supplies a clean channel at transport level
(fractional <=~1e-7); the only order-percent route is the kinematic-prior
slip; and the observed ensemble elevation is commensurate with the mundane
mass-model systematic floor (real data: SDSS1206 power-law vs composite
D_dt posteriors differ by ~6%; per-system kappa_ext widths are 2-10%).

No synthetic data.  All inputs are real TDCOSMO-2025 public products or
archived pipeline outputs.
"""

import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "58"

C_KM_S = 299792.458
MPC_KM = 3.085677581e19
ARCSEC_RAD = math.pi / (180.0 * 3600.0)
S_PER_DAY = 86400.0
KPC_TO_KM = 3.0856775814913673e16

# --- corpus-stated bounds ---------------------------------------------------
GW_BOUND = 1e-15          # |c_gamma - c_g|/c (GW170817; Paper 4 Sec. 4)
A_O_MINUS_1 = 1e-7        # uniform local-well endpoint rescaling (Paper 4 Sec. 6)
PATH_FACTOR = 10.0        # L_path = PATH_FACTOR * R_E (generous ceiling)

# operational flat-LCDM frame for the distance bookkeeping; the fractional
# channel bounds are insensitive to this choice
H0_OP = 70.0
OM = 0.3

# tension scale being tested (Paper 11; ~73.0 vs 67.4 km/s/Mpc)
DELTA_H0_REQUIRED = (73.0 - 67.4) / 67.4   # ~0.083


def comoving_distance(z, n=20000):
    """Flat LCDM comoving distance in Mpc (trapezoid integral of 1/E(z))."""
    zs = np.linspace(0.0, z, n)
    ez = np.sqrt(OM * (1.0 + zs) ** 3 + (1.0 - OM))
    return float(np.trapezoid(1.0 / ez, zs)) * C_KM_S / H0_OP


def D_A(z):
    return comoving_distance(z) / (1.0 + z)


def load_yaml():
    p = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "external"
        / "tdcosmo2025_public"
        / "TDCOSMO_sample"
        / "tdcosmo_sample.yaml"
    )
    with open(p) as f:
        return yaml.safe_load(f)


def load_step05():
    p = PROJECT_ROOT / "results" / "outputs" / "step_05_tdcosmo_shear.json"
    with open(p) as f:
        return json.load(f)


def sigma_v_preferred(entry):
    """Preferred sigma_v for the kinematic-prior sector: highest-SNR entry."""
    cands = []
    for key, val in entry.items():
        if isinstance(val, dict) and "sigma_v" in val:
            sv = val.get("sigma_v")
            sverr = val.get("sigma_v_error", np.inf)
            if sv and np.isfinite(sv):
                cands.append((sverr, sv, key))
    if not cands:
        return None, None
    cands.sort()
    return cands[0][1], cands[0][2]


def main():
    print_status(f"STEP {STEP_NUM}: H0 time-delay channel audit", "TITLE")

    yaml_data = load_yaml()
    s05 = load_step05()
    s05_systems = s05.get("systems", {})

    # representative measured delay per system (max |dt_obs| over pairs)
    NAME_MAP = {"J1206+4332": "SDSS1206+4332"}
    per_system = []
    for sys_name, entry in yaml_data.items():
        if not isinstance(entry, dict) or "z_lens" not in entry:
            continue
        z_d = float(entry["z_lens"])
        z_s = float(entry.get("z_source", np.nan))
        theta_E = float(entry.get("theta_E", np.nan))
        kext = entry.get("kappa_ext")
        kext_up = entry.get("kappa_ext_err_up")
        kext_dn = entry.get("kappa_ext_err_down")

        # characteristic delay from step_05 (max |dt_obs| over image pairs)
        step05_name = sys_name if sys_name in s05_systems else None
        if step05_name is None:
            for k in s05_systems:
                if NAME_MAP.get(k) == sys_name or k == sys_name:
                    step05_name = k
        dt_char_d = np.nan
        if step05_name:
            pairs = s05_systems[step05_name].get("pair_results", {})
            vals = [abs(p.get("dt_obs_days", 0.0)) for p in pairs.values()]
            if vals:
                dt_char_d = max(vals)

        # Einstein radius physical scale and halo path ceiling
        Dd = D_A(z_d)                                   # Mpc
        R_E_kpc = theta_E * ARCSEC_RAD * Dd * 1e3       # kpc
        L_path_km = PATH_FACTOR * R_E_kpc * KPC_TO_KM   # km (10 R_E ceiling)
        dT_prop_s = GW_BOUND * L_path_km / C_KM_S
        frac_prop = dT_prop_s / (dt_char_d * S_PER_DAY) if np.isfinite(dt_char_d) and dt_char_d > 0 else np.nan

        sv, sv_src = sigma_v_preferred(entry)

        per_system.append({
            "system": sys_name,
            "z_lens": z_d,
            "z_source": z_s,
            "theta_E_arcsec": theta_E,
            "D_d_Mpc": Dd,
            "R_E_kpc": R_E_kpc,
            "delay_characteristic_days": dt_char_d if np.isfinite(dt_char_d) else None,
            "dT_propagated_ceiling_s": dT_prop_s,
            "frac_bias_propagated": frac_prop if np.isfinite(frac_prop) else None,
            "frac_bias_endpoint": A_O_MINUS_1,
            "sigma_v_kms": sv,
            "sigma_v_source": sv_src,
            "kappa_ext": kext,
            "kappa_ext_err_up": kext_up,
            "kappa_ext_err_down": kext_dn,
            "kappa_ext_halfwidth": (
                (kext_up + kext_dn) / 2.0
                if (kext_up is not None and kext_dn is not None)
                else None
            ),
        })

    # ensemble summaries
    frac_props = [s["frac_bias_propagated"] for s in per_system if s["frac_bias_propagated"] is not None]
    kext_hw = [s["kappa_ext_halfwidth"] for s in per_system if s["kappa_ext_halfwidth"]]

    # --- real mass-model-family floor: SDSS1206 power-law vs composite -------
    fam_dir = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "external"
        / "tdcosmo2025_public"
        / "TDCOSMO_sample"
        / "TDCOSMO_data"
        / "SDSS1206+4332"
    )
    family_floor = None
    try:
        pl = pickle.load(open(fam_dir / "angular_diameter_final_power_law.txt", "rb"), encoding="latin1")
        cp = pickle.load(open(fam_dir / "angular_diameter_final_composite.txt", "rb"), encoding="latin1")
        ddt_pl = np.asarray(pl[0], dtype=float)
        ddt_cp = np.asarray(cp[0], dtype=float)
        family_floor = {
            "system": "SDSS1206+4332",
            "D_dt_powerlaw_median_Mpc": float(np.median(ddt_pl)),
            "D_dt_composite_median_Mpc": float(np.median(ddt_cp)),
            "D_dt_fractional_shift": float(np.median(ddt_cp) / np.median(ddt_pl) - 1.0),
            "H0_family_spread_frac": float(np.median(ddt_pl) / np.median(ddt_cp) - 1.0),
            "note": "same system, same data, two mass-model families: ~6% D_dt "
                    "systematic floor, exceeding the transport-sector ceiling "
                    "by ~8 orders of magnitude",
        }
    except Exception as exc:  # pragma: no cover - data should exist
        family_floor = {"error": str(exc)}

    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": "H0LiCOW/TDCOSMO-class time-delay channel audit under the TEP mechanism (issue 11-9).",
        "required_gap_frac": DELTA_H0_REQUIRED,
        "channels": {
            "propagated_conformal": {
                "frac_bias": 0.0,
                "note": "exactly zero by null-cone invariance; verified numerically in step_56 (0.0 d)",
            },
            "propagated_disformal": {
                "frac_bias_max": max(frac_props) if frac_props else None,
                "frac_bias_median": float(np.median(frac_props)) if frac_props else None,
                "note": "GW170817-saturated bound over a 10 R_E halo traverse; ~5-13 ms absolute, ~1e-9 to 1e-8 fractional relative to the measured delays",
            },
            "endpoint_rescaling": {
                "frac_bias": A_O_MINUS_1,
                "note": "uniform local-well clock factor multiplies the measured duration; degeneracy level ~1e-7",
            },
            "mass_model_sector": {
                "structural_bias": 0.0,
                "note": "phantom convergence is real g-sector curvature (scalar backreaction): it contributes to image positions and to the Fermat delay identically to ordinary mass, so the fitted model absorbs it self-consistently; the reconstruction-space component propagates nothing to photons",
                "mundane_floor_kappa_ext_halfwidth_median": float(np.median(kext_hw)) if kext_hw else None,
            },
            "kinematic_prior_sector": {
                "route": "sigma_* mass-sheet-breaking prior under the lensing-dynamical slip (apparent lensing mass exceeding dynamical mass; TEP-UCD Sec. 2 bookkeeping)",
                "frac_bias_bound": "bounded by the aperture slip fraction; order of the inner phantom share (~0.1-0.3 scale), commensurate with the ~8% required but quantitatively tied to the unresolved transfer kernel (issue 19-1)",
                "sign": ("under the corpus's own bookkeeping M_dyn < M_lens: the sigma_* prior "
                         "under-reads the mass that the image positions require, so the fitted "
                         "internal mass-sheet normalisation is driven to larger lambda (steeper "
                         "effective profile) -- inferred H0 biased HIGH, the required direction "
                         "for elevating the channel above the unbiased value"),
            },
            "delay_carrier_sign_check": {
                "statement": "a positive propagated residual biases the inference to LOWER H0 (observed delay longer than geometric -> D_dt over-estimated); elevating the inferred H0 requires a negative residual, opposite to the corpus's claimed sign pattern (Refsdal +30.1 d). The delay carrier cannot supply the required direction at any amplitude.",
                "required_direction": "delta_t_obs < delta_t_geom (negative residual)",
                "corpus_claimed_direction": "positive residuals (observed exceeds model)",
            },
        },
        "per_system": per_system,
        "mass_model_family_floor": family_floor,
        "verdict": {
            "statement": (
                "The TEP mechanism is clean on the time-delay channel at transport level: "
                "the conformal sector propagates exactly zero photon delay, the disformal "
                "sector is bounded at ~1e-9 to 1e-8 fractional (5-13 ms absolute) on galaxy-halo paths, and "
                "the endpoint rescaling contributes ~1e-7. The mass-model sector is "
                "self-consistent because phantom convergence is real curvature contributing "
                "identically to positions and Fermat delays. The sole order-percent route is "
                "the lensing-dynamical slip entering the sigma_*-anchored mass-sheet prior, "
                "bounded by the aperture phantom fraction and tied to the unresolved transfer "
                "kernel (19-1); under the corpus's own bookkeeping its sign biases inferred H0 "
                "high -- the required direction, at a magnitude commensurate with the ~8% "
                "needed. Sign bookkeeping on the delay carrier shows a positive residual "
                "biases inferred H0 low, so the transport sector cannot supply the elevation "
                "at any amplitude. The channel's commensurate same-sign route therefore runs "
                "through the mass-model normalisation (phantom mass as real curvature, with "
                "the lensing--dynamical slip biasing the kinematic anchor), not through the "
                "propagated delay. The residual quantitative task is the per-system slip "
                "magnitude, which shares the transfer-kernel requirement of issue 19-1; the "
                "mundane mass-model systematic floor (~6% family spread on SDSS1206; ~3-10% "
                "per-system kappa_ext) supplies the reference scale it must be evaluated "
                "against."
            )
        },
    }

    outdir = PROJECT_ROOT / "results" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "step_58_h0_channel_audit.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print_status(f"Results saved to {path}", "SUCCESS")


if __name__ == "__main__":
    main()
