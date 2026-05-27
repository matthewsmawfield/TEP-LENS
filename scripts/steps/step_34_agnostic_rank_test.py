#!/usr/bin/env python3
"""
TEP-LENS: Step 34 - Agnostic Rank Test

Purpose: Test whether delays correlate with potential depth (traced by inferred
convergence kappa) without assuming the log-magnification ansatz.
This separates the physical claim ("delays scale with potential depth")
from the parametric claim ("the scaling follows 1 + alpha*log10(mu)").

Algorithm:
1. Use Kelly et al. (2023) measured delays for S1-S4-SX.
2. Use the inferred kappa distribution from Step 32 as the potential-depth tracer.
3. Compute Spearman rank correlation between delays and kappa_median.
4. Compute permutation p-value under random delay orderings.
5. Compare to the same test using mu_proxy instead of kappa.

With only n=5, the test is underpowered. The value lies in showing that even
an ansatz-free rank test points in the same direction (SX shallow, SX late).

Outputs:
- results/outputs/step_34_agnostic_rank_test.json
"""

import json
import sys
from pathlib import Path
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "34"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def permutation_pvalue(delays, tracer, n_perm=20000, seed=42):
    """
    One-sided permutation test for Spearman correlation.
    H0: delays and tracer are independent.
    H1: positive correlation (later delay = deeper potential = higher tracer).
    """
    rng = np.random.default_rng(seed)
    rho_obs, _ = stats.spearmanr(delays, tracer)

    count = 0
    for _ in range(n_perm):
        perm_delays = rng.permutation(delays)
        rho_perm, _ = stats.spearmanr(perm_delays, tracer)
        if rho_perm >= rho_obs:
            count += 1

    return float(rho_obs), float(count / n_perm)


def main():
    print_status(f"STEP {STEP_NUM}: Agnostic Rank Test", "TITLE")

    catalog_path = PROJECT_ROOT / "data" / "raw" / "sn_lensing" / "lensed_sn_catalog.json"
    with open(catalog_path) as f:
        catalog = json.load(f)

    refsdal = catalog["sn_refsdal"]
    delays_data = refsdal["time_delays_days"]
    fluxes = refsdal["magnification_proxies"]["flux_total"]

    # Measured delays relative to S1
    images = ["S1", "S2", "S3", "S4", "SX"]
    delays = np.array([
        0.0,
        delays_data["dt_S2_S1"]["value"],
        delays_data["dt_S3_S1"]["value"],
        delays_data["dt_S4_S1"]["value"],
        delays_data["dt_SX_S1"]["value"],
    ])

    # Load step_32 kappa inference
    step32_path = PROJECT_ROOT / "results" / "outputs" / "step_32_kappa_proxy_validation.json"
    kappa_medians = None
    kappa_p16 = None
    kappa_p84 = None
    if step32_path.exists():
        with open(step32_path) as f:
            s32 = json.load(f)
        kappa_summary = s32["sn_refsdal_sensitivity"]["kappa_summary_per_image"]
        kappa_medians = np.array([kappa_summary[img]["median"] for img in images])
        kappa_p16 = np.array([kappa_summary[img]["p16"] for img in images])
        kappa_p84 = np.array([kappa_summary[img]["p84"] for img in images])
    else:
        print_status("Step 32 output not found; using hardcoded kappa medians", "WARN")
        # Fallback values from step_32 median output
        kappa_medians = np.array([0.093, 0.029, 0.010, 0.181, 0.010])
        kappa_p16 = np.array([0.000, 0.000, 0.000, 0.001, 0.010])
        kappa_p84 = np.array([0.267, 0.206, 0.151, 0.348, 0.057])

    # Magnification proxies for comparison
    mu_rel = np.array([fluxes[img]["value"] for img in images])
    mu_norm = mu_rel / np.mean(mu_rel)

    # ------------------------------------------------------------------
    # Rank test 1: delays vs inferred kappa (ansatz-free potential depth)
    # ------------------------------------------------------------------
    rho_kappa, p_kappa_perm = permutation_pvalue(delays, kappa_medians, n_perm=20000)

    print_status(f"Delay vs inferred kappa:")
    print_status(f"  Spearman rho = {rho_kappa:.3f}")
    print_status(f"  Permutation p (one-sided, n=20000) = {p_kappa_perm:.4f}")

    # ------------------------------------------------------------------
    # Rank test 2: delays vs mu_proxy (parametric proxy)
    # ------------------------------------------------------------------
    rho_mu, p_mu_perm = permutation_pvalue(delays, mu_norm, n_perm=20000)

    print_status(f"Delay vs flux-proxy mu:")
    print_status(f"  Spearman rho = {rho_mu:.3f}")
    print_status(f"  Permutation p (one-sided, n=20000) = {p_mu_perm:.4f}")

    # ------------------------------------------------------------------
    # Rank test 3: inner cross only (S1-S4), excluding SX
    # ------------------------------------------------------------------
    inner_idx = [0, 1, 2, 3]  # S1-S4
    rho_kappa_inner, p_kappa_inner = permutation_pvalue(
        delays[inner_idx], kappa_medians[inner_idx], n_perm=10000
    )
    rho_mu_inner, p_mu_inner = permutation_pvalue(
        delays[inner_idx], mu_norm[inner_idx], n_perm=10000
    )

    print_status(f"Inner cross only (S1-S4) delay vs kappa:")
    print_status(f"  Spearman rho = {rho_kappa_inner:.3f}, p = {p_kappa_inner:.4f}")
    print_status(f"Inner cross only (S1-S4) delay vs mu:")
    print_status(f"  Spearman rho = {rho_mu_inner:.3f}, p = {p_mu_inner:.4f}")

    # ------------------------------------------------------------------
    # Bootstrap over kappa uncertainty: draw kappa from [p16, p84] envelope
    # ------------------------------------------------------------------
    rng = np.random.default_rng(42)
    n_boot = 10000
    boot_rhos = []
    for _ in range(n_boot):
        k_boot = rng.uniform(kappa_p16, kappa_p84)
        r_b, _ = stats.spearmanr(delays, k_boot)
        boot_rhos.append(r_b)
    boot_rhos = np.array(boot_rhos)
    rho_kappa_lo = float(np.percentile(boot_rhos, 2.5))
    rho_kappa_hi = float(np.percentile(boot_rhos, 97.5))
    p_kappa_median = float(np.median([permutation_pvalue(delays, k_boot, n_perm=1000)[1]
                                      for k_boot in [rng.uniform(kappa_p16, kappa_p84)
                                                     for _ in range(1000)]]))

    print_status(f"Kappa-uncertainty bootstrap (n={n_boot}):")
    print_status(f"  Spearman rho 95% CI: [{rho_kappa_lo:.3f}, {rho_kappa_hi:.3f}]")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    output = {
        "step": STEP_NUM,
        "status": "success",
        "all_five_images": {
            "images": images,
            "delays_days": delays.tolist(),
            "kappa_median": kappa_medians.tolist(),
            "mu_norm": mu_norm.tolist(),
            "delay_vs_kappa": {
                "spearman_rho": rho_kappa,
                "permutation_p_onesided": p_kappa_perm,
            },
            "delay_vs_mu": {
                "spearman_rho": rho_mu,
                "permutation_p_onesided": p_mu_perm,
            },
            "kappa_uncertainty_bootstrap": {
                "rho_95ci": [rho_kappa_lo, rho_kappa_hi],
            },
            "interpretation": (
                "With n=5, the ansatz-free rank test is underpowered. "
                "The delay-kappa correlation (rho={:.3f}, p={:.3f}) points in the same "
                "direction as the proxy model (SX shallow, SX late), but is not significant. "
                "The physical claim (delays scale with potential depth) and the parametric "
                "claim (log-magnification ansatz) are therefore not separately constrained by "
                "present data. Only the sign of the S4-SX contrast is probative."
                .format(rho_kappa, p_kappa_perm)
            ),
        },
        "inner_cross_only": {
            "images": [images[i] for i in inner_idx],
            "delay_vs_kappa": {
                "spearman_rho": rho_kappa_inner,
                "permutation_p_onesided": p_kappa_inner,
            },
            "delay_vs_mu": {
                "spearman_rho": rho_mu_inner,
                "permutation_p_onesided": p_mu_inner,
            },
            "interpretation": (
                "The inner cross alone shows no rank correlation with either kappa or mu. "
                "This confirms the proxy model is non-predictive without SX."
            ),
        },
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_agnostic_rank_test.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
