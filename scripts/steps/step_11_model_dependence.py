#!/usr/bin/env python3
"""
TEP-LENS: Step 11 - Model Dependence Robustness

Purpose:
1) Quantify sensitivity of headline evidence to inter-model dependence.
2) Run leave-one-out (LOO) stress tests.
3) Add exact/permutation-style tests for sign and weighted-mean residual.

This step does not introduce new astrophysical assumptions; it stress-tests
inference built from step_07 outputs.
"""

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "11"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def method_family(method_str: str) -> str:
    s = method_str.lower()
    if "lenstool" in s:
        return "LENSTOOL"
    if "glee" in s:
        return "GLEE"
    if "glafic" in s:
        return "GLAFIC"
    if "ltm" in s:
        return "LTM"
    if "wslap" in s:
        return "WSLAP+"
    return "OTHER"


def weighted_mean_and_sigma(values: np.ndarray, sigmas: np.ndarray):
    w = 1.0 / (sigmas**2)
    mu = float(np.sum(w * values) / np.sum(w))
    sigma_mu = float(1.0 / np.sqrt(np.sum(w)))
    return mu, sigma_mu


def exact_sign_flip_pvalue(values: np.ndarray, sigmas: np.ndarray):
    """
    Exact two-sided style tail under symmetric sign-flip null for weighted mean.
    Uses all 2^N sign assignments on |values|.

    DEPRECATED in v0.2: This test assumes independent, symmetric residuals.
    It is invalid under correlated lens-model systematics. Kept for traceability
    but excluded from headline results.
    """
    abs_vals = np.abs(values)
    obs_mu, _ = weighted_mean_and_sigma(values, sigmas)

    n = len(values)
    tails = 0
    total = 1 << n

    for mask in range(total):
        signs = np.ones(n)
        for i in range(n):
            if (mask >> i) & 1:
                signs[i] = -1.0
        trial = signs * abs_vals
        mu_trial, _ = weighted_mean_and_sigma(trial, sigmas)
        if mu_trial >= obs_mu:
            tails += 1

    p_one_sided = tails / total
    return float(p_one_sided)


def exact_nonzero_wilcoxon_greater(values: np.ndarray) -> float:
    """One-sided exact Wilcoxon on the manuscript's non-zero residual set."""
    nonzero = np.asarray(values, dtype=float)
    nonzero = nonzero[nonzero != 0]
    if nonzero.size == 0:
        return float("nan")

    # Exact Wilcoxon with Pratt zeros is undefined in SciPy; dropping exact
    # zeros matches the headline non-zero-residual test definition.
    return float(
        stats.wilcoxon(
            nonzero,
            alternative="greater",
            zero_method="wilcox",
            mode="exact",
        ).pvalue
    )


def wilcoxon_positive_rank_statistic(values: np.ndarray) -> float:
    """Wilcoxon W+ statistic on non-zero residuals, with average tied ranks."""
    nonzero = np.asarray(values, dtype=float)
    nonzero = nonzero[nonzero != 0]
    if nonzero.size == 0:
        return 0.0
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    return float(np.sum(ranks[nonzero > 0]))


def wilcoxon_block_bootstrap(deltas, method_families, n_bootstrap=20000, seed=20260607):
    """
    Bootstrap Wilcoxon p-values respecting method-family clustering.
    Models within the same family are treated as correlated blocks.
    Families are resampled with replacement; within each drawn family,
    individual models are also resampled with replacement.
    """
    rng = np.random.default_rng(seed)
    families = np.array(method_families)
    unique_fams = np.unique(families)

    p_values = []
    for _ in range(n_bootstrap):
        boot_fams = rng.choice(unique_fams, size=len(unique_fams), replace=True)
        boot_indices = []
        for fam in boot_fams:
            fam_idx = np.where(families == fam)[0]
            boot_idx = rng.choice(fam_idx, size=len(fam_idx), replace=True)
            boot_indices.extend(boot_idx)

        boot_deltas = deltas[boot_indices]
        p_value = exact_nonzero_wilcoxon_greater(boot_deltas)
        p_values.append(float(p_value) if np.isfinite(p_value) else 1.0)

    p_vals_arr = np.array(p_values)
    return {
        "p_median": float(np.median(p_vals_arr)),
        "p_16": float(np.percentile(p_vals_arr, 16)),
        "p_84": float(np.percentile(p_vals_arr, 84)),
        "fraction_le_0_05": float(np.mean(p_vals_arr <= 0.05)),
        "n_bootstrap": n_bootstrap,
    }


def wilcoxon_permutation_test(deltas, method_families):
    """
    Exact family-sign-flip test for the Wilcoxon statistic under method-family dependence.

    Null: signs are exchangeable at method-family level. We enumerate every
    possible family sign assignment (preserving intra-family correlation) and
    recompute the Wilcoxon statistic. This is exact under the sharp null and
    requires no superpopulation assumption, unlike the bootstrap.
    """
    families = np.array(method_families)
    unique_fams = np.unique(families)

    # Observed Wilcoxon W+ statistic (one-sided), computed directly to avoid
    # small-sample normal-approximation warnings.
    obs_w = wilcoxon_positive_rank_statistic(deltas)

    n_ge = 0
    n_total = 0
    for family_signs in product((-1.0, 1.0), repeat=len(unique_fams)):
        signs = np.ones(len(deltas))
        for fam, sign in zip(unique_fams, family_signs):
            signs[families == fam] = sign
        perm_deltas = signs * np.abs(deltas)
        perm_w = wilcoxon_positive_rank_statistic(perm_deltas)
        if perm_w >= obs_w:
            n_ge += 1
        n_total += 1

    p_one_sided = n_ge / n_total
    return {
        "p_value_one_sided": float(p_one_sided),
        "n_permutations": int(n_total),
        "n_extreme_or_more": int(n_ge),
        "n_method_families": int(len(unique_fams)),
        "method_families": [str(fam) for fam in unique_fams],
        "observed_statistic": obs_w,
        "description": (
            "Exact family-sign-flip test: enumerates all method-family sign assignments. "
            "Exact under sharp null; no superpopulation assumption."
        ),
    }


def beta_binom_tail_p(n_pos: int, n_total: int, rho: float):
    """
    Correlated sign-test null via beta-binomial with p=0.5 and ICC=rho.
    rho=0 recovers ordinary binomial.
    """
    if rho <= 1e-12:
        return float(stats.binomtest(n_pos, n_total, 0.5, alternative="greater").pvalue)

    # For beta-binomial: rho = 1 / (alpha + beta + 1) when alpha=beta for p=0.5.
    # Solve alpha+beta = 1/rho - 1, and alpha=beta.
    a_plus_b = (1.0 / rho) - 1.0
    alpha = 0.5 * a_plus_b
    beta = 0.5 * a_plus_b

    xs = np.arange(n_pos, n_total + 1)
    return float(np.sum(stats.betabinom.pmf(xs, n_total, alpha, beta)))


def main():
    print_status(f"STEP {STEP_NUM}: Model Dependence Robustness", "TITLE")

    step07_path = PROJECT_ROOT / "results" / "outputs" / "step_07_observed_vs_predicted.json"
    with open(step07_path) as f:
        s07 = json.load(f)

    models = s07["per_model_results"]
    deltas = np.array([m["delta_obs_minus_pred_days"] for m in models], dtype=float)
    sigmas = np.array([m["sigma_total_days"] for m in models], dtype=float)
    names = [m["name"] for m in models]

    n_total = len(models)
    n_pos = int(np.sum(deltas > 0))

    # Baseline metrics
    mu_all, sigma_all = weighted_mean_and_sigma(deltas, sigmas)
    z_all = mu_all / sigma_all
    p_binom_all = float(stats.binomtest(n_pos, n_total, 0.5, alternative="greater").pvalue)
    p_flip_exact = exact_sign_flip_pvalue(deltas, sigmas)
    p_wilcoxon = exact_nonzero_wilcoxon_greater(deltas)

    print_status(f"Baseline weighted mean residual: {mu_all:+.2f} ± {sigma_all:.2f} d (z={z_all:.2f})")
    print_status(f"Baseline sign test: {n_pos}/{n_total} positive, p={p_binom_all:.4f}")
    print_status(f"Wilcoxon signed-rank p(one-sided): {p_wilcoxon:.4f}")
    print_status(f"[DEPRECATED] Exact sign-flip p={p_flip_exact:.4f} — assumes independent residuals; invalid under correlation")

    # Blind-only baseline
    blind_mask = np.array([m.get("blind", True) for m in models], dtype=bool)
    blind_deltas = deltas[blind_mask]
    n_blind_total = int(blind_mask.sum())
    n_blind_pos = int(np.sum(blind_deltas > 0))
    n_nonzero_all = int(np.sum(deltas != 0))
    n_nonzero_blind = int(np.sum(blind_deltas != 0))
    p_binom_blind = float(stats.binomtest(n_blind_pos, n_blind_total, 0.5, alternative="greater").pvalue)
    p_wilcoxon_blind = exact_nonzero_wilcoxon_greater(blind_deltas)

    print_status(f"Delay-blind sign test: {n_blind_pos}/{n_blind_total} positive, p={p_binom_blind:.4f}")
    print_status(f"Delay-blind Wilcoxon signed-rank p(one-sided): {p_wilcoxon_blind:.4f}")

    # Strictly blind original-values baseline
    strictly_blind_original_deltas = np.array([
        m["delta_obs_minus_pred_original_days"]
        for m in models if m.get("blind", True)
    ], dtype=float)
    strictly_blind_original_sigmas = np.array([
        m["sigma_total_days"]
        for m in models if m.get("blind", True)
    ], dtype=float)
    n_strictly_blind_pos = int(np.sum(strictly_blind_original_deltas > 0))
    n_strictly_blind_total = len(strictly_blind_original_deltas)
    p_binom_strictly_blind = float(stats.binomtest(
        n_strictly_blind_pos, n_strictly_blind_total, 0.5, alternative="greater"
    ).pvalue)
    p_wilcoxon_strictly_blind = exact_nonzero_wilcoxon_greater(strictly_blind_original_deltas)
    print_status(f"Strictly blind sign test (original values): {n_strictly_blind_pos}/{n_strictly_blind_total} positive, p={p_binom_strictly_blind:.4f}")
    print_status(f"Strictly blind Wilcoxon (original values): p={p_wilcoxon_strictly_blind:.4f}")

    # Method-family effective sample size (Kish-style weighting)
    fams = [method_family(m.get("method", "")) for m in models]
    fam_counts = {}
    for fam in fams:
        fam_counts[fam] = fam_counts.get(fam, 0) + 1

    # Each model weighted by inverse family size
    fam_weights = np.array([1.0 / fam_counts[f] for f in fams], dtype=float)
    neff_family = float((np.sum(fam_weights) ** 2) / np.sum(fam_weights**2))

    print_status(f"Method-family effective sample size (Kish proxy): N_eff={neff_family:.2f} from N={n_total}")

    # Correlation sensitivity curve for sign-test p-values (all models)
    rho_grid = np.linspace(0.0, 0.8, 17)
    corr_curve = []
    for rho in rho_grid:
        p_corr = beta_binom_tail_p(n_pos, n_total, float(rho))
        neff_rho = n_total / (1.0 + (n_total - 1.0) * rho)
        corr_curve.append(
            {
                "rho": float(rho),
                "n_eff": float(neff_rho),
                "p_beta_binom_tail": float(p_corr),
                "z_equiv": float(stats.norm.isf(min(max(p_corr, 1e-12), 1 - 1e-12))),
            }
        )

    # Correlation sensitivity curve for delay-blind sign-test p-values
    corr_curve_blind = []
    for rho in rho_grid:
        p_corr = beta_binom_tail_p(n_blind_pos, n_blind_total, float(rho))
        neff_rho = n_blind_total / (1.0 + (n_blind_total - 1.0) * rho)
        corr_curve_blind.append(
            {
                "rho": float(rho),
                "n_eff": float(neff_rho),
                "p_beta_binom_tail": float(p_corr),
                "z_equiv": float(stats.norm.isf(min(max(p_corr, 1e-12), 1 - 1e-12))),
            }
        )

    # Correlation sensitivity curve for strictly blind sign-test p-values
    corr_curve_strictly_blind = []
    for rho in rho_grid:
        p_corr = beta_binom_tail_p(n_strictly_blind_pos, n_strictly_blind_total, float(rho))
        neff_rho = n_strictly_blind_total / (1.0 + (n_strictly_blind_total - 1.0) * rho)
        corr_curve_strictly_blind.append(
            {
                "rho": float(rho),
                "n_eff": float(neff_rho),
                "p_beta_binom_tail": float(p_corr),
                "z_equiv": float(stats.norm.isf(min(max(p_corr, 1e-12), 1 - 1e-12))),
            }
        )

    # Break-even ICC: correlation at which p-value reaches 0.05
    def find_break_even_rho(n_pos_sub, n_total_sub, target_p=0.05, n_grid=10001):
        rhos = np.linspace(0.0, 0.99, n_grid)
        ps = np.array([beta_binom_tail_p(n_pos_sub, n_total_sub, float(r)) for r in rhos])
        # ps is monotonically increasing with rho (higher correlation -> higher p)
        # Find first index where ps >= target_p
        idx = np.searchsorted(ps, target_p, side="left")
        if idx == 0:
            return 0.0  # p already >= target_p at rho=0
        if idx >= len(rhos):
            return None  # p never reaches target_p
        # Linear interpolation between rhos[idx-1] and rhos[idx]
        p_lo, p_hi = ps[idx - 1], ps[idx]
        if p_hi == p_lo:
            return float(rhos[idx - 1])
        frac = (target_p - p_lo) / (p_hi - p_lo)
        return float(rhos[idx - 1] + frac * (rhos[idx] - rhos[idx - 1]))

    rho_break_even_all = find_break_even_rho(n_pos, n_total)
    rho_break_even_delay_blind = find_break_even_rho(n_blind_pos, n_blind_total)
    rho_break_even_strictly_blind = find_break_even_rho(n_strictly_blind_pos, n_strictly_blind_total)

    n_nonzero_strictly_blind = int(np.sum(strictly_blind_original_deltas != 0))

    # Wilcoxon break-even (approximate): n_eff = n / (1 + (n-1)*rho) under ICC,
    # and p ≈ 1/2^{n_eff} for all-positive ranks.
    def wilcoxon_break_even_rho(n_nonzero, target_p=0.05):
        n_eff_target = -np.log2(target_p)
        if n_eff_target >= n_nonzero:
            return 0.0
        rho = (n_nonzero / n_eff_target - 1.0) / (n_nonzero - 1.0)
        return float(min(rho, 0.99))

    rho_break_even_wilcoxon_all = wilcoxon_break_even_rho(n_nonzero_all)
    rho_break_even_wilcoxon_delay_blind = wilcoxon_break_even_rho(n_nonzero_blind)
    rho_break_even_wilcoxon_strictly_blind = wilcoxon_break_even_rho(n_nonzero_strictly_blind)

    print_status(f"Break-even ICC (all {n_total} binomial {n_pos}/{n_total}): rho = {rho_break_even_all:.3f}")
    print_status(f"Break-even ICC (delay-blind {n_blind_total} binomial {n_blind_pos}/{n_blind_total}): rho = {rho_break_even_delay_blind:.3f}")
    print_status(f"Break-even ICC (strictly blind {n_strictly_blind_total} binomial {n_strictly_blind_pos}/{n_strictly_blind_total}): rho = {rho_break_even_strictly_blind:.3f}")

    # Family-sign-flip test (exact under sharp null, no superpopulation assumption)
    blind_fams = [f for i, f in enumerate(fams) if blind_mask[i]]
    perm_all = wilcoxon_permutation_test(deltas, fams)
    perm_delay_blind = wilcoxon_permutation_test(blind_deltas, blind_fams)
    perm_strictly_blind = wilcoxon_permutation_test(
        strictly_blind_original_deltas, blind_fams
    )
    print_status(
        f"Exact family-sign-flip Wilcoxon (all {n_total}): p = {perm_all['p_value_one_sided']:.4f} "
        f"({perm_all['n_extreme_or_more']}/{perm_all['n_permutations']})"
    )
    print_status(
        f"Exact family-sign-flip Wilcoxon (delay-blind {perm_delay_blind['n_method_families']} families): p = {perm_delay_blind['p_value_one_sided']:.4f} "
        f"({perm_delay_blind['n_extreme_or_more']}/{perm_delay_blind['n_permutations']})"
    )
    print_status(
        f"Exact family-sign-flip Wilcoxon (strictly blind {perm_strictly_blind['n_method_families']} families): p = {perm_strictly_blind['p_value_one_sided']:.4f} "
        f"({perm_strictly_blind['n_extreme_or_more']}/{perm_strictly_blind['n_permutations']})"
    )

    # Block-bootstrap Wilcoxon (realistic method-family dependence)
    boot_all = wilcoxon_block_bootstrap(deltas, fams, n_bootstrap=20000)
    boot_delay_blind = wilcoxon_block_bootstrap(blind_deltas, blind_fams, n_bootstrap=20000)
    boot_strictly_blind = wilcoxon_block_bootstrap(
        strictly_blind_original_deltas, blind_fams, n_bootstrap=20000
    )

    print_status(f"Block-bootstrap Wilcoxon (all {n_total}, family-aware): "
                 f"p_median={boot_all['p_median']:.4f}, "
                 f"p_16-84=[{boot_all['p_16']:.4f}, {boot_all['p_84']:.4f}], "
                 f"P(p<=0.05)={boot_all['fraction_le_0_05']:.2%}")
    n_blind_fams = len(np.unique(blind_fams))
    print_status(f"Block-bootstrap Wilcoxon (delay-blind {n_blind_fams} families, family-aware): "
                 f"p_median={boot_delay_blind['p_median']:.4f}, "
                 f"p_16-84=[{boot_delay_blind['p_16']:.4f}, {boot_delay_blind['p_84']:.4f}], "
                 f"P(p<=0.05)={boot_delay_blind['fraction_le_0_05']:.2%}")
    print_status(f"Block-bootstrap Wilcoxon (strictly blind {n_blind_fams} families, family-aware): "
                 f"p_median={boot_strictly_blind['p_median']:.4f}, "
                 f"p_16-84=[{boot_strictly_blind['p_16']:.4f}, {boot_strictly_blind['p_84']:.4f}], "
                 f"P(p<=0.05)={boot_strictly_blind['fraction_le_0_05']:.2%}")

    # DEPRECATED: Wilcoxon break-even via n_eff scaling is mathematically invalid.
    # The block-bootstrap above is the primary dependence-aware test.
    print_status(f"[CAVEAT] Wilcoxon break-even ICC (all {n_total}): rho ≈ {rho_break_even_wilcoxon_all:.3f} — "
                 f"approximate, not a closed-form correction; use block-bootstrap instead.")
    print_status(f"[CAVEAT] Wilcoxon break-even ICC (delay-blind): rho ≈ {rho_break_even_wilcoxon_delay_blind:.3f} — "
                 f"approximate, not a closed-form correction; use block-bootstrap instead.")
    print_status(f"[CAVEAT] Wilcoxon break-even ICC (strictly blind): rho ≈ {rho_break_even_wilcoxon_strictly_blind:.3f} — "
                 f"approximate, not a closed-form correction; use block-bootstrap instead.")
    print_status("Interpretation: if inter-model correlation exceeds these values, p > 0.05.")

    # Leave-one-out stress tests
    loo = []
    R_tep = s07["tep_prediction"]["R_tep_prediction_days"]
    for i in range(n_total):
        mask = np.ones(n_total, dtype=bool)
        mask[i] = False

        d_sub = deltas[mask]
        s_sub = sigmas[mask]
        n_sub = len(d_sub)
        n_pos_sub = int(np.sum(d_sub > 0))

        mu_sub, sigma_sub = weighted_mean_and_sigma(d_sub, s_sub)
        z_sub = mu_sub / sigma_sub
        p_sign_sub = float(stats.binomtest(n_pos_sub, n_sub, 0.5, alternative="greater").pvalue)

        chi2_gr = float(np.sum((d_sub / s_sub) ** 2))
        chi2_tep = float(np.sum(((d_sub - R_tep) / s_sub) ** 2))
        dchi2 = chi2_gr - chi2_tep

        loo.append(
            {
                "left_out": names[i],
                "weighted_mean_days": float(mu_sub),
                "weighted_sigma_days": float(sigma_sub),
                "z_from_gr": float(z_sub),
                "n_positive": int(n_pos_sub),
                "n_total": int(n_sub),
                "p_sign": float(p_sign_sub),
                "delta_chi2_gr_minus_tep": float(dchi2),
            }
        )

    p_vals = np.array([x["p_sign"] for x in loo])
    z_vals = np.array([x["z_from_gr"] for x in loo])
    dchi2_vals = np.array([x["delta_chi2_gr_minus_tep"] for x in loo])

    loo_summary = {
        "p_sign_min": float(np.min(p_vals)),
        "p_sign_max": float(np.max(p_vals)),
        "z_min": float(np.min(z_vals)),
        "z_max": float(np.max(z_vals)),
        "delta_chi2_min": float(np.min(dchi2_vals)),
        "delta_chi2_max": float(np.max(dchi2_vals)),
        "n_cases_tep_better": int(np.sum(dchi2_vals > 0)),
        "n_cases_total": int(len(dchi2_vals)),
    }

    print_status("LOO robustness summary:")
    print_status(f"  Sign-test p range: {loo_summary['p_sign_min']:.4f} to {loo_summary['p_sign_max']:.4f}")
    print_status(f"  Weighted-mean z range: {loo_summary['z_min']:.2f} to {loo_summary['z_max']:.2f}")
    print_status(
        f"  Δχ² (GR-TEP) range: {loo_summary['delta_chi2_min']:+.2f} to {loo_summary['delta_chi2_max']:+.2f}"
    )

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZE

        set_pub_style()

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIG_SIZE)

        # Correlation sensitivity
        ax1.plot(
            [c["rho"] for c in corr_curve],
            [c["p_beta_binom_tail"] for c in corr_curve],
            marker="o",
            lw=1.8,
        )
        ax1.axhline(p_binom_all, color="gray", ls="--", lw=1.2, label="Independent binomial p")
        ax1.set_xlabel("Assumed inter-model ICC (rho)")
        ax1.set_ylabel("One-sided sign-test p-value")
        ax1.set_title("Dependence sensitivity (beta-binomial)")
        ax1.set_ylim(0, 1)
        ax1.grid(alpha=0.3, ls=":")
        ax1.legend()

        # LOO sign p-values
        order = np.argsort(p_vals)
        ax2.barh(np.arange(n_total), p_vals[order], color="#4c78a8")
        ax2.axvline(0.05, color="black", ls="--", lw=1.0)
        ax2.set_yticks(np.arange(n_total))
        ax2.set_yticklabels([loo[i]["left_out"] for i in order], )
        ax2.set_xlabel("LOO one-sided sign-test p-value")
        ax2.set_title("Leave-one-out stress test")
        ax2.grid(alpha=0.3, ls=":", axis="x")

        fig.tight_layout()
        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_model_dependence.png"
        fig.savefig(out_fig)
        plt.close(fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")
        out_fig = None

    # Save
    out = {
        "step": STEP_NUM,
        "status": "success",
        "baseline": {
            "weighted_mean_days": float(mu_all),
            "weighted_sigma_days": float(sigma_all),
            "z_from_gr": float(z_all),
            "n_positive": int(n_pos),
            "n_total": int(n_total),
            "p_binomial_independent": float(p_binom_all),
            "p_wilcoxon_signed_rank": float(p_wilcoxon),
        },
        "delay_blind_baseline": {
            "n_positive": int(n_blind_pos),
            "n_total": int(n_blind_total),
            "p_binomial_independent": float(p_binom_blind),
            "p_wilcoxon_signed_rank": float(p_wilcoxon_blind),
            "n_nonzero": int(n_nonzero_blind),
            "note": "Delay-blind tier: revised values fixed before measured delay known.",
        },
        "strictly_blind_baseline": {
            "n_positive": int(n_strictly_blind_pos),
            "n_total": int(n_strictly_blind_total),
            "p_binomial_independent": float(p_binom_strictly_blind),
            "p_wilcoxon_signed_rank": float(p_wilcoxon_strictly_blind),
            "n_nonzero": int(n_nonzero_strictly_blind),
            "note": "Strictly blind tier: original pre-reappearance published predictions.",
        },
        "effective_sample_size": {
            "n_total": int(n_total),
            "method_family_counts": fam_counts,
            "n_eff_family_kish": float(neff_family),
        },
        "dependence_sensitivity_curve": corr_curve,
        "dependence_sensitivity_curve_delay_blind": corr_curve_blind,
        "dependence_sensitivity_curve_strictly_blind": corr_curve_strictly_blind,
        "dependence_aware_test_hierarchy": {
            "bridge_statement": (
                "Because the Wilcoxon statistic lacks a closed-form variance under "
                "exchangeable intra-class correlation, an exact family-sign-flip test "
                "(enumerating all method-family sign assignments) provides the most "
                "rigorous dependence-aware bound. The block-bootstrap is retained as "
                "a sensitivity exploration but is demoted from operational primary."
            ),
            "tier_1a_primary_independence": "wilcoxon_signed_rank_delay_blind",
            "tier_1b_primary_correlation_aware": "permutation_test_wilcoxon_delay_blind",
            "tier_1c_correlation_aware_floor": "permutation_test_wilcoxon_strictly_blind",
            "tier_2_corroborating_sign_test": "beta_binomial_delay_blind",
            "sensitivity_exploration": "block_bootstrap_wilcoxon",
        },
        "permutation_test_wilcoxon": {
            "description": (
                "Exact family-sign-flip test for the Wilcoxon statistic. Enumerates "
                "all method-family sign assignments; exact under sharp null; no "
                "superpopulation assumption. More rigorous than bootstrap for tiny samples."
            ),
            "all_models": perm_all,
            "delay_blind": perm_delay_blind,
            "strictly_blind": perm_strictly_blind,
        },
        "block_bootstrap_wilcoxon": {
            "description": "Family-aware block-bootstrap respecting method-family clusters. Sensitivity exploration (not operational primary).",
            "all_models": boot_all,
            "delay_blind": boot_delay_blind,
            "strictly_blind": boot_strictly_blind,
        },
        "break_even_icc": {
            "description": "Inter-model correlation (ICC) at which one-sided p-value reaches 0.05. Two distinct methods are reported: (1) exact beta-binomial for sign-test p-values, and (2) an approximate n_eff scaling heuristic for Wilcoxon p-values that lacks formal justification.",
            "beta_binomial_sign_test": {
                "description": "Exact beta-binomial ICC break-even for binomial sign-test p-values. These are rigorous.",
                "binomial_all": {
                    "n_positive": int(n_pos),
                    "n_total": int(n_total),
                    "rho_break_even": rho_break_even_all,
                },
                "binomial_delay_blind": {
                    "n_positive": int(n_blind_pos),
                    "n_total": int(n_blind_total),
                    "rho_break_even": rho_break_even_delay_blind,
                },
                "binomial_strictly_blind": {
                    "n_positive": int(n_strictly_blind_pos),
                    "n_total": int(n_strictly_blind_total),
                    "rho_break_even": rho_break_even_strictly_blind,
                },
            },
            "wilcoxon_approx_heuristic": {
                "description": f"Approximate ICC break-even for Wilcoxon p-values via n_eff = n / (1 + (n-1)*rho) and p ≈ 1/2^n_eff. This is a mathematically unfounded heuristic, NOT a beta-binomial result. The exact family-sign-flip test (delay-blind p={perm_delay_blind['p_value_one_sided']:.4f}) is the primary correlation-aware rank bound; strictly blind floor is p={perm_strictly_blind['p_value_one_sided']:.4f}.",
                "wilcoxon_all_approx": {
                    "n_nonzero": int(n_nonzero_all),
                    "rho_break_even": rho_break_even_wilcoxon_all,
                    "caveat": "Mathematically unfounded approximation. Use exact permutation_test_wilcoxon instead.",
                },
                "wilcoxon_delay_blind_approx": {
                    "n_nonzero": int(n_nonzero_blind),
                    "rho_break_even": rho_break_even_wilcoxon_delay_blind,
                    "caveat": "Mathematically unfounded approximation. Use exact permutation_test_wilcoxon instead.",
                },
                "wilcoxon_strictly_blind_approx": {
                    "n_nonzero": int(n_nonzero_strictly_blind),
                    "rho_break_even": rho_break_even_wilcoxon_strictly_blind,
                    "caveat": "Mathematically unfounded approximation. Use exact permutation_test_wilcoxon instead.",
                },
            },
        },
        "deprecated_tests": {
            "exact_sign_flip_pvalue": {
                "value": float(p_flip_exact),
                "caveat": "Assumes independent, symmetric residuals. Invalid under correlated lens-model systematics. Not reported in headline results.",
            },
        },
        "leave_one_out": loo,
        "leave_one_out_summary": loo_summary,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_model_dependence.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
