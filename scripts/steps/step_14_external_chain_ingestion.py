#!/usr/bin/env python3
"""
TEP-LENS: Step 14 - External Chain Ingestion (H0LiCOW/TDCOSMO ecosystem)

This step ingests publicly available strong-lens distance chains from the
H0LiCOW public repository and summarizes immediate coverage for the TEP-LENS
quasar cross-check set.

Public source:
- https://github.com/shsuyu/H0LiCOW-public/tree/master/h0licow_distance_chains
"""

import json
import sys
import warnings
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status

STEP_NUM = "14"
API_URL = "https://api.github.com/repos/shsuyu/H0LiCOW-public/contents/h0licow_distance_chains"

# Keep downloads bounded for reproducibility/runtime.
MAX_FILE_MB = 30.0
MAX_ROWS_PARSE = 250000


def safe_json_default(obj):
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def fetch_json(url: str):
    req = Request(url, headers={"User-Agent": "TEP-LENS-step14"})
    with urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


def download_file(url: str, target: Path):
    req = Request(url, headers={"User-Agent": "TEP-LENS-step14"})
    with urlopen(req, timeout=60) as r:
        data = r.read()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)


def parse_numeric_summary(path: Path):
    # Parse robustly while suppressing known blank-line warnings from numpy.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        if path.suffix.lower() == ".csv":
            arr = np.genfromtxt(path, delimiter=",", skip_header=1, max_rows=MAX_ROWS_PARSE)
        else:
            # For whitespace/tab-delimited files, let genfromtxt infer columns.
            arr = np.genfromtxt(path, comments="#", max_rows=MAX_ROWS_PARSE)

    # Handle degenerate cases from genfromtxt.
    if arr is None or (isinstance(arr, float) and np.isnan(arr)):
        raise ValueError("no numeric data parsed")

    arr = np.asarray(arr)
    if arr.size == 0:
        raise ValueError("empty numeric array")

    if arr.ndim == 1:
        arr = arr[:, None]

    # Some .dat files are actually CSV with a header (e.g., Dt,weight).
    if np.isnan(arr).all() or np.isnan(arr[:, 0]).all():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            arr_retry = np.genfromtxt(path, delimiter=",", skip_header=1, max_rows=MAX_ROWS_PARSE)
        arr_retry = np.asarray(arr_retry)
        if arr_retry.size > 0:
            if arr_retry.ndim == 1:
                arr_retry = arr_retry[:, None]
            arr = arr_retry

    col0 = arr[:, 0]
    return {
        "n_rows_used": int(arr.shape[0]),
        "n_cols": int(arr.shape[1]),
        "col0_median": float(np.median(col0)),
        "col0_p16": float(np.percentile(col0, 16)),
        "col0_p84": float(np.percentile(col0, 84)),
        "col0_std": float(np.std(col0)),
    }


def build_des0408_equivalent_bridge(project_root: Path):
    """
    DES0408 is not present as a standalone chain file in the H0LiCOW-public
    distance-chain folder. Build an external-equivalent product from the
    published DES0408 delay/flux table already curated in step_05.
    """
    p = project_root / "results" / "outputs" / "step_05_tdcosmo_shear.json"
    if not p.exists():
        return None

    with open(p) as f:
        s05 = json.load(f)

    des = s05.get("systems", {}).get("DES0408-5354")
    if not des:
        return None

    pairs = des.get("pair_results", {})
    abs_shifts = [abs(v.get("tep_predicted_shift_days", 0.0)) for v in pairs.values()]

    return {
        "source": "step_05_tdcosmo_shear.json (published DES0408 literature values)",
        "type": "external_equivalent_literature_product",
        "metadata": des.get("metadata", {}),
        "n_pairs": len(pairs),
        "median_abs_shift_days": float(np.median(abs_shifts)) if abs_shifts else None,
        "pair_keys": list(pairs.keys()),
    }


def main():
    print_status(f"STEP {STEP_NUM}: External Chain Ingestion", "TITLE")

    out_dir = PROJECT_ROOT / "data" / "interim" / "external" / "h0licow_distance_chains"
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        entries = fetch_json(API_URL)
    except Exception as e:
        print_status(f"Failed to query GitHub API: {e}", "ERROR")
        entries = []

    files = [e for e in entries if isinstance(e, dict) and e.get("type") == "file"]

    lens_aliases = {
        "HE0435": ["HE0435"],
        "WFI2033": ["wfi2033", "WFI2033"],
        "DES0408": ["DES0408", "0408"],
    }

    manifest = []
    ingested = []

    for e in files:
        name = e.get("name", "")
        size_mb = float(e.get("size", 0)) / 1e6
        download_url = e.get("download_url")

        rec = {
            "name": name,
            "size_mb": size_mb,
            "download_url": download_url,
            "downloaded": False,
            "skipped_reason": None,
            "summary": None,
        }

        if not download_url:
            rec["skipped_reason"] = "missing_download_url"
            manifest.append(rec)
            continue

        if size_mb > MAX_FILE_MB:
            rec["skipped_reason"] = f"file_too_large_gt_{MAX_FILE_MB:.0f}MB"
            manifest.append(rec)
            continue

        target = out_dir / name
        try:
            download_file(download_url, target)
            rec["summary"] = parse_numeric_summary(target)
            rec["downloaded"] = True
            ingested.append(name)
        except Exception as e:
            rec["skipped_reason"] = f"download_or_parse_error: {e}"

        manifest.append(rec)

    # Coverage wrt this manuscript's quasar systems.
    coverage = {}
    ingested_lower = [n.lower() for n in ingested]
    for lens, aliases in lens_aliases.items():
        coverage[lens] = any(any(a.lower() in n for a in aliases) for n in ingested_lower)

    n_cov = int(sum(1 for v in coverage.values() if v))

    des0408_bridge = build_des0408_equivalent_bridge(PROJECT_ROOT)
    coverage_with_bridge = dict(coverage)
    if des0408_bridge is not None:
        coverage_with_bridge["DES0408"] = True
    n_cov_with_bridge = int(sum(1 for v in coverage_with_bridge.values() if v))

    print_status(f"API files discovered: {len(files)}")
    print_status(f"Files ingested (<= {MAX_FILE_MB:.0f} MB): {len(ingested)}")
    print_status(f"Direct chain coverage for HE0435/WFI2033/DES0408: {n_cov}/3")
    for lens, ok in coverage.items():
        print_status(f"  {lens}: {'YES' if ok else 'NO'}")
    if des0408_bridge is not None:
        print_status(f"Coverage with DES0408 equivalent literature bridge: {n_cov_with_bridge}/3")

    out_fig = None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scripts.utils.plot_style import set_pub_style, COLORS, FIG_SIZES, save_fig

        labels = ["HE0435", "WFI2033", "DES0408"]
        direct_vals = [1 if coverage[k] else 0 for k in labels]
        bridge_vals = [1 if coverage_with_bridge[k] else 0 for k in labels]
        x = np.arange(len(labels))
        w = 0.36

        set_pub_style()
        fig, ax = plt.subplots(figsize=FIG_SIZES["web_standard"])
        ax.bar(x - w / 2, direct_vals, width=w, color=COLORS["gr"], label="Direct chain ingestion")
        ax.bar(x + w / 2, bridge_vals, width=w, color=COLORS["tep"], label="With literature bridge")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, 1.2)
        ax.set_ylabel("Coverage flag")
        ax.set_title("Step 14 external data coverage for TEP-LENS quasar cross-check set")
        ax.legend()

        out_fig = PROJECT_ROOT / "results" / "figures" / f"step_{STEP_NUM}_external_coverage.png"
        save_fig(fig, out_fig)
        print_status(f"Figure saved to {out_fig}")
    except Exception as e:
        print_status(f"Plotting failed: {e}", "ERROR")

    out = {
        "step": STEP_NUM,
        "status": "success",
        "source": {
            "api": API_URL,
            "max_file_mb": MAX_FILE_MB,
            "max_rows_parse": MAX_ROWS_PARSE,
        },
        "files_discovered": len(files),
        "files_ingested": len(ingested),
        "ingested_names": ingested,
        "coverage_tdcosmo_triplet": coverage,
        "coverage_count": n_cov,
        "des0408_equivalent_bridge": des0408_bridge,
        "coverage_with_bridge": coverage_with_bridge,
        "coverage_with_bridge_count": n_cov_with_bridge,
        "manifest": manifest,
        "figure": str(out_fig) if out_fig else None,
    }

    out_path = PROJECT_ROOT / "results" / "outputs" / f"step_{STEP_NUM}_external_chain_ingestion.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {out_path}")
    print_status(f"Step {STEP_NUM} complete.")


if __name__ == "__main__":
    main()
