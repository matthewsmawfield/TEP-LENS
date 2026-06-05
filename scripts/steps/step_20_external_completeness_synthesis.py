#!/usr/bin/env python3
"""
TEP-LENS: Step 20 - External Completeness Synthesis

Consumes step_16 (correlated significance) and step_19 (TDCOSMO2025 standalone
ingestion) to produce an explicit external-data completeness assessment.
Aligns with the single-test benchmark paradigm, eschewing invalid meta-analytic
combinations of correlated tests.

Output:
- results/outputs/step_20_external_completeness_synthesis.json
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "20"


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def main():
    print_status(f"STEP {STEP_NUM}: External Completeness Synthesis", "TITLE")

    s16_path = PROJECT_ROOT / "results" / "outputs" / "step_16_tier_significance.json"
    s19_path = PROJECT_ROOT / "results" / "outputs" / "step_19_tdcosmo2025_ingestion.json"

    if not s16_path.exists() or not s19_path.exists():
        missing = []
        if not s16_path.exists():
            missing.append(str(s16_path))
        if not s19_path.exists():
            missing.append(str(s19_path))
        raise FileNotFoundError(f"Required inputs missing: {missing}")

    with open(s16_path) as f:
        s16 = json.load(f)
    with open(s19_path) as f:
        s19 = json.load(f)

    local = s19.get("local_summary", {})
    lens_cov = local.get("lens_coverage", {})
    n_lenses_total = len(lens_cov)
    n_lenses_cov = int(sum(1 for v in lens_cov.values() if v))
    lens_coverage_fraction = float(n_lenses_cov / n_lenses_total) if n_lenses_total > 0 else 0.0

    n_files = int(local.get("n_files", 0))
    n_post = int(local.get("posterior_like_count", 0))
    posterior_density = float(n_post / n_files) if n_files > 0 else 0.0

    interpretation = {
        "external_readiness": "complete" if lens_coverage_fraction >= 0.99 else "partial",
        "lens_coverage_fraction": lens_coverage_fraction,
        "posterior_like_count": n_post,
        "posterior_density": posterior_density,
        "note": (
            "Step 19 provides full lens-name coverage for the target TDCOSMO set. "
            "However, combining external-informed tests with core sign tests is "
            "statistically invalid due to underlying data correlation. The recommended "
            "reporting track is the single most robust non-parametric test: "
            "the exact family-sign-flip test respecting method-family clusters, "
            "supported by the independence-primary Wilcoxon and external consistency."
        ),
    }

    out = {
        "step": STEP_NUM,
        "status": "success",
        "description": "External completeness synthesis from step_16 + step_19",
        "inputs": {
            "step_16": str(s16_path.relative_to(PROJECT_ROOT)),
            "step_19": str(s19_path.relative_to(PROJECT_ROOT)),
        },
        "external_data_completeness": {
            "n_lenses_total": n_lenses_total,
            "n_lenses_covered": n_lenses_cov,
            "lens_coverage_fraction": lens_coverage_fraction,
            "n_files_local": n_files,
            "n_posterior_like": n_post,
            "posterior_density": posterior_density,
        },
        "recommended_reporting_track": {
            "name": "single_robust_test_benchmark",
            "result": s16.get("test_hierarchy", {}).get("tier_1b_primary_correlation_aware", {}),
        },
        "interpretation": interpretation,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / "step_20_external_completeness_synthesis.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Lens coverage: {n_lenses_cov}/{n_lenses_total} ({lens_coverage_fraction:.2f})")
    print_status(
        f"Recommended Reporting Track: {out['recommended_reporting_track']['name']} "
        f"(z={out['recommended_reporting_track']['result'].get('z_score', 0):.3f})"
    )
    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
