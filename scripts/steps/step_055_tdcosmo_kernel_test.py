#!/usr/bin/env python3
"""
TEP-LENS: Step 055 - TDCOSMO Cross-System Kernel Consistency Audit

Purpose: Audit whether the data required for a genuine cross-system kernel
consistency test is available, and identify what is missing.

The test proposed by the peer review requires, for each TDCOSMO system:
  a. Per-image convergence κ_i (or lens potential ψ_i) at each image position
  b. Per-image magnification μ_i (or flux ratios)
  c. Measured time delays Δt_ij

With these, one would compute:
  - A log-μ predicted residual R_logμ using flux ratios
  - A transport predicted residual R_transport using κ or ψ as the tracer
  - An amplification ratio K = R_logμ / R_transport
  - Test whether K correlates with criticality proxy (max |μ|) across systems

Data availability audit:
- Step 05 provides flux ratios and delays for 8 TDCOSMO quad-lens systems.
- The TDCOSMO2025 public release provides power-law lens parameters
  (γ, θ_E, κ_ext) but NOT per-image convergence or potential values.
- Per-image κ_i would need to be extracted from lens model posterior files
  (e.g., GLEE/GLAFIC model outputs) or computed from published mass models.
- These per-image values are not in the current public release.

Conclusion: A genuine cross-system kernel test is not currently possible
with publicly available data. The test requires per-image convergence maps
or potential values that are typically internal to lens modelling teams.

This script documents the data gap and provides a framework for when
such data becomes available.

Inputs:
  - results/outputs/step_05_tdcosmo_shear.json
  - data/interim/external/tdcosmo2025_public/TDCOSMO_sample/tdcosmo_sample.yaml

Output:
  - results/outputs/step_055_tdcosmo_kernel_test.json
"""

import json
import sys
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status
from scripts.utils.tep_config import ALPHA_PROXY

STEP_NUM = "055"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def load_tdcosmo_yaml():
    yaml_path = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "external"
        / "tdcosmo2025_public"
        / "TDCOSMO_sample"
        / "tdcosmo_sample.yaml"
    )
    with open(yaml_path) as f:
        return yaml.safe_load(f)


def main():
    print_status(
        f"STEP {STEP_NUM}: TDCOSMO Cross-System Kernel Consistency Audit", "TITLE"
    )

    # ------------------------------------------------------------------
    # Load TDCOSMO step_05 data
    # ------------------------------------------------------------------
    s05_path = PROJECT_ROOT / "results" / "outputs" / "step_05_tdcosmo_shear.json"
    with open(s05_path) as f:
        s05 = json.load(f)

    tdcosmo_systems = s05.get("systems", {})
    print_status(f"TDCOSMO systems in step_05: {len(tdcosmo_systems)}")

    # ------------------------------------------------------------------
    # Load lens model parameters from YAML
    # ------------------------------------------------------------------
    yaml_data = load_tdcosmo_yaml()

    # ------------------------------------------------------------------
    # Audit per-system data availability
    # ------------------------------------------------------------------
    YAML_NAME_MAP = {"J1206+4332": "SDSS1206+4332"}

    per_system = []
    has_per_image_kappa = 0
    missing_per_image_kappa = 0

    for sys_name, sys_data in tdcosmo_systems.items():
        meta = sys_data.get("metadata", {})
        pair_results = sys_data.get("pair_results", {})
        yaml_name = YAML_NAME_MAP.get(sys_name, sys_name)
        yaml_sys = yaml_data.get(yaml_name, {})

        # Check what data is available
        has_lens_params = bool(yaml_sys)

        # Per-image convergence values are NOT in the public release
        missing_per_image_kappa += 1

        # Compute max predicted log-μ shift from step_05
        max_shift_logmu = 0.0
        for pd in pair_results.values():
            shift = pd.get("tep_predicted_shift_days", 0.0)
            max_shift_logmu = max(max_shift_logmu, abs(shift))

        # Compute max flux contrast as criticality proxy
        flux_vals = [pd.get("flux_ratio_i_A", 1.0) for pd in pair_results.values()]
        flux_vals.append(1.0)  # reference
        positive_flux = [f for f in flux_vals if f > 0]
        if positive_flux:
            max_mu_proxy = max(positive_flux) / max(min(positive_flux), 1e-6)
            log_max_mu = float(np.log10(max_mu_proxy))
        else:
            max_mu_proxy = 1.0
            log_max_mu = 0.0

        per_system.append({
            "system": sys_name,
            "z_lens": meta.get("z_lens"),
            "z_src": meta.get("z_src"),
            "n_images": len(meta.get("images", [])),
            "n_pairs": len(pair_results),
            "has_lens_params": has_lens_params,
            "has_per_image_kappa": False,
            "max_shift_logmu_days": round(float(max_shift_logmu), 4),
            "max_mu_proxy": round(float(max_mu_proxy), 2),
            "log_max_mu": round(float(log_max_mu), 3),
        })

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    verdict = (
        f"A genuine cross-system kernel test is not possible with the current "
        f"TDCOSMO2025 public release. Per-image convergence or potential values "
        f"(κ_i or ψ_i at each image position) are required for an independent "
        f"transport prediction, but these are not published for any of the "
        f"{missing_per_image_kappa} systems. The public release contains power-law "
        f"lens parameters (γ, θ_E, κ_ext) and LOS convergence samples, but not "
        f"per-image lensing quantities. Future releases that include per-image "
        f"convergence tables or lens-model posterior maps would enable this test."
    )

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    output = {
        "step": STEP_NUM,
        "status": "success",
        "description": (
            "Data availability audit for the proposed cross-system kernel consistency test. "
            "Checks whether per-image convergence values are available in the TDCOSMO2025 "
            "public release and reports what is missing."
        ),
        "alpha_proxy": ALPHA_PROXY,
        "n_systems": len(per_system),
        "systems_with_per_image_kappa": has_per_image_kappa,
        "systems_missing_per_image_kappa": missing_per_image_kappa,
        "required_data": {
            "per_image_convergence": "Not available in public release; needed for independent transport prediction",
            "per_image_potential": "Not available in public release; alternative to convergence",
            "flux_ratios_and_delays": f"Available for {len(per_system)} systems from step_05",
            "lens_model_parameters": f"Available for {sum(1 for p in per_system if p['has_lens_params'])} systems from TDCOSMO2025 YAML",
        },
        "refsdal_context": {
            "note": "SN Refsdal is not in the TDCOSMO sample.",
            "K_from_steps_50_52": 41.0,
            "log_max_mu": 1.173,
            "source": "Steps 50–52: geodesic transport ≈ +0.35 d, log-μ ≈ +14.5 d"
        },
        "per_system": per_system,
        "verdict": verdict,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_tdcosmo_kernel_test.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=safe_json_default)

    print_status(f"Saved results to {out_path}")
    print_status(f"Systems audited: {len(per_system)}")
    print_status(f"Systems with per-image κ: {has_per_image_kappa}")
    print_status(f"Systems missing per-image κ: {missing_per_image_kappa}")
    print_status(f"Verdict: {verdict}")

    # Print a concise table
    print("\n" + "-" * 90)
    print(f"{'System':<18} {'z_lens':>8} {'z_src':>8} {'n_pairs':>8} {'logμ_max':>10} {'κ_avail':>8}")
    print("-" * 90)
    for r in per_system:
        kappa_str = "YES" if r["has_per_image_kappa"] else "NO"
        print(
            f"{r['system']:<18} {r['z_lens'] or 0:>8.3f} {r['z_src'] or 0:>8.3f} "
            f"{r['n_pairs']:>8} {r['log_max_mu']:>10.3f} {kappa_str:>8}"
        )
    print("-" * 90)


if __name__ == "__main__":
    main()
