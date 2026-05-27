#!/usr/bin/env python3
"""
TEP-LENS: Step 18 - External Dataset Registry and Ingestion Readiness

Builds a machine-readable registry of external datasets relevant to TEP-LENS,
including local staging/readiness checks and immediate integration paths.

This step is offline-safe and does not perform any network I/O.

Output:
- results/outputs/step_18_external_dataset_registry.json
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status, safe_json_default

STEP_NUM = "18"




def summarize_path(path: Path, max_examples: int = 8):
    if not path.exists():
        return {
            "exists": False,
            "is_dir": False,
            "n_files": 0,
            "total_size_bytes": 0,
            "example_files": [],
        }

    if path.is_file():
        return {
            "exists": True,
            "is_dir": False,
            "n_files": 1,
            "total_size_bytes": int(path.stat().st_size),
            "example_files": [str(path.name)],
        }

    files = [p for p in path.rglob("*") if p.is_file()]
    files_sorted = sorted(files)
    return {
        "exists": True,
        "is_dir": True,
        "n_files": len(files_sorted),
        "total_size_bytes": int(sum(p.stat().st_size for p in files_sorted)),
        "example_files": [str(p.relative_to(path)) for p in files_sorted[:max_examples]],
    }


def main():
    print_status(f"STEP {STEP_NUM}: External Dataset Registry", "TITLE")

    h0licow_cache = PROJECT_ROOT / "data" / "interim" / "external" / "h0licow_distance_chains"
    tdcosmo2025_stage = PROJECT_ROOT / "data" / "interim" / "external" / "tdcosmo2025_public"
    snh0pe_raw = PROJECT_ROOT / "data" / "raw" / "snh0pe"

    datasets = [
        {
            "name": "H0LiCOW public distance chains",
            "priority": 1,
            "measurement_type": ["strong_lensing_time_delay", "distance_chains"],
            "independence_role": "External lensing-chain prior stress test (already used in step_14/15).",
            "public_source": "https://github.com/shsuyu/H0LiCOW-public/tree/master/h0licow_distance_chains",
            "local_path": str(h0licow_cache.relative_to(PROJECT_ROOT)),
            "local_status": summarize_path(h0licow_cache),
            "integration_steps": ["step_14_external_chain_ingestion.py", "step_15_external_informed_inflation.py"],
        },
        {
            "name": "TDCOSMO 2025 public likelihood package",
            "priority": 2,
            "measurement_type": ["strong_lensing_time_delay", "hierarchical_likelihoods", "H0_posteriors"],
            "independence_role": "Independent multi-lens cosmography update for external validation.",
            "public_source": "https://github.com/TDCOSMO/TDCOSMO2025_public",
            "paper": "https://www.aanda.org/articles/aa/full_html/2025/12/aa55801-25/aa55801-25.html",
            "local_path": str(tdcosmo2025_stage.relative_to(PROJECT_ROOT)),
            "local_status": summarize_path(tdcosmo2025_stage),
            "integration_steps": ["planned: step_19_tdcosmo2025_ingestion.py"],
        },
        {
            "name": "SN H0pe updates",
            "priority": 3,
            "measurement_type": ["SN_Refsdal_time_delays", "H0_inference_cross_check"],
            "independence_role": "Observed-delay and prediction updates supporting primary SN-Refsdal evidence.",
            "public_source": "manuscript-linked SN H0pe references (curated local files)",
            "local_path": str(snh0pe_raw.relative_to(PROJECT_ROOT)),
            "local_status": summarize_path(snh0pe_raw),
            "integration_steps": ["step_01_fetch_snh0pe_data.py", "step_07_observed_vs_predicted.py"],
        },
    ]

    ready = [d["name"] for d in datasets if d["local_status"]["exists"]]
    missing = [d["name"] for d in datasets if not d["local_status"]["exists"]]

    output = {
        "step": STEP_NUM,
        "status": "success",
        "description": "External dataset shortlist and local ingestion readiness for TEP-LENS",
        "datasets": datasets,
        "summary": {
            "n_datasets": len(datasets),
            "n_ready_local": len(ready),
            "ready_now": ready,
            "missing_local": missing,
            "next_action": (
                "If TDCOSMO2025 is staged under data/interim/external/tdcosmo2025_public, "
                "implement step_19_tdcosmo2025_ingestion.py to compute lens-level coverage and posterior summaries."
            ),
        },
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_external_dataset_registry.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Ready local datasets: {len(ready)}/{len(datasets)}")
    for name in missing:
        print_status(f"Missing local staging: {name}", "WARNING")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
