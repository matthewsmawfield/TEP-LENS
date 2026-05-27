#!/usr/bin/env python3
"""
TEP-LENS: Step 19 - TDCOSMO 2025 Public Ingestion (Standalone)

Standalone ingestion step for TDCOSMO 2025 public materials.

Behavior:
1) Local-first: if files already exist under data/interim/external/tdcosmo2025_public,
   use them directly.
2) If local cache is missing/empty, discover files via GitHub API and download a
   bounded subset of analysis-relevant public files.
3) Build machine-readable coverage and lightweight numeric summaries.

Output:
- results/outputs/step_19_tdcosmo2025_ingestion.json
"""

import json
import sys
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import quote

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import print_status, safe_json_default

STEP_NUM = "19"

REPO = "TDCOSMO/TDCOSMO2025_public"
API_ROOT = f"https://api.github.com/repos/{REPO}/contents"
TARGET_DIR = PROJECT_ROOT / "data" / "interim" / "external" / "tdcosmo2025_public"
OUTPUT_JSON = PROJECT_ROOT / "results" / "outputs" / "step_19_tdcosmo2025_ingestion.json"

MAX_FILE_MB = 20.0
MAX_DOWNLOAD_FILES = 60
MAX_PARSE_FILES = 20
MAX_ROWS_PARSE = 200000

ALLOWED_EXT = {".json", ".csv", ".dat", ".txt", ".md", ".yml", ".yaml", ".npz"}
RELEVANCE_TOKENS = {
    "posterior", "chain", "sample", "samples", "h0", "kappa", "lens", "likelihood", "tdcosmo"
}
LENS_ALIASES = {
    "HE0435": ["he0435", "0435"],
    "WFI2033": ["wfi2033", "2033"],
    "RXJ1131": ["rxj1131", "1131"],
    "PG1115": ["pg1115", "1115"],
    "B1608": ["b1608", "1608"],
    "DES0408": ["des0408", "0408"],
    "J1206": ["j1206", "1206"],
}




def fetch_json(url: str):
    req = Request(url, headers={"User-Agent": "TEP-LENS-step19"})
    with urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode("utf-8"))


def download_file(url: str, target: Path):
    req = Request(url, headers={"User-Agent": "TEP-LENS-step19"})
    with urlopen(req, timeout=120) as r:
        data = r.read()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)


def list_repo_files(path: str = ""):
    """Recursively list files from GitHub contents API."""
    url = API_ROOT if not path else f"{API_ROOT}/{quote(path)}"
    data = fetch_json(url)

    files = []
    if isinstance(data, dict) and data.get("type") == "file":
        return [data]
    if not isinstance(data, list):
        return files

    for entry in data:
        etype = entry.get("type")
        if etype == "file":
            files.append(entry)
        elif etype == "dir":
            sub = list_repo_files(entry.get("path", ""))
            files.extend(sub)
    return files


def local_files(root: Path):
    if not root.exists():
        return []
    return sorted([p for p in root.rglob("*") if p.is_file()])


def relevance_score(path_str: str):
    s = path_str.lower()
    score = 0
    for tok in RELEVANCE_TOKENS:
        if tok in s:
            score += 1
    ext = Path(path_str).suffix.lower()
    if ext in {".json", ".csv", ".dat", ".npz"}:
        score += 1
    return score


def parse_numeric_summary(path: Path):
    ext = path.suffix.lower()
    if ext not in {".csv", ".dat", ".txt"}:
        return None

    try:
        if ext == ".csv":
            arr = np.genfromtxt(path, delimiter=",", skip_header=1, max_rows=MAX_ROWS_PARSE)
        else:
            arr = np.genfromtxt(path, comments="#", max_rows=MAX_ROWS_PARSE)

        arr = np.asarray(arr)
        if arr.size == 0:
            return None
        if arr.ndim == 1:
            arr = arr[:, None]

        # Retry CSV interpretation for dat/txt that are comma-delimited.
        if np.isnan(arr).all():
            arr_retry = np.genfromtxt(path, delimiter=",", skip_header=1, max_rows=MAX_ROWS_PARSE)
            arr_retry = np.asarray(arr_retry)
            if arr_retry.size > 0:
                if arr_retry.ndim == 1:
                    arr_retry = arr_retry[:, None]
                arr = arr_retry

        col0 = arr[:, 0]
        finite_col0 = col0[np.isfinite(col0)]
        if finite_col0.size == 0:
            return None
        return {
            "n_rows_used": int(arr.shape[0]),
            "n_cols": int(arr.shape[1]),
            "col0_median": float(np.nanmedian(finite_col0)),
            "col0_p16": float(np.nanpercentile(finite_col0, 16)),
            "col0_p84": float(np.nanpercentile(finite_col0, 84)),
            "col0_std": float(np.nanstd(finite_col0)),
        }
    except Exception:
        return None


def compute_lens_coverage(paths_rel):
    coverage = {}
    names = [p.lower() for p in paths_rel]
    for lens, aliases in LENS_ALIASES.items():
        coverage[lens] = any(any(a in n for a in aliases) for n in names)
    return coverage


def summarize_local_dataset(root: Path):
    files = local_files(root)
    rel = [str(p.relative_to(root)) for p in files]
    ext_counts = {}
    total_size = 0
    posterior_like = []

    for p, r in zip(files, rel):
        ext = p.suffix.lower() or "<noext>"
        ext_counts[ext] = ext_counts.get(ext, 0) + 1
        total_size += p.stat().st_size
        if any(tok in r.lower() for tok in ["posterior", "chain", "sample", "likelihood"]):
            posterior_like.append(r)

    lens_cov = compute_lens_coverage(rel)
    parse_targets = [p for p in files if p.suffix.lower() in {".csv", ".dat", ".txt"}][:MAX_PARSE_FILES]
    numeric_summaries = []
    for p in parse_targets:
        summary = parse_numeric_summary(p)
        if summary is not None:
            numeric_summaries.append({
                "file": str(p.relative_to(root)),
                "summary": summary,
            })

    return {
        "exists": root.exists(),
        "n_files": len(files),
        "total_size_bytes": int(total_size),
        "extension_counts": ext_counts,
        "posterior_like_files": posterior_like[:100],
        "posterior_like_count": len(posterior_like),
        "lens_coverage": lens_cov,
        "lens_coverage_count": int(sum(1 for v in lens_cov.values() if v)),
        "example_files": rel[:40],
        "numeric_file_summaries": numeric_summaries,
    }


def main():
    print_status(f"STEP {STEP_NUM}: TDCOSMO2025 Standalone Ingestion", "TITLE")

    TARGET_DIR.mkdir(parents=True, exist_ok=True)
    remote_manifest = []
    downloaded = []
    source_mode = "local_cache"
    remote_error = None

    existing = local_files(TARGET_DIR)
    if len(existing) == 0:
        source_mode = "remote_fetch"
        print_status("No local TDCOSMO2025 cache found. Attempting GitHub API discovery...", "PROCESS")
        try:
            repo_files = list_repo_files("")
            # Build simplified manifest
            for f in repo_files:
                remote_manifest.append({
                    "path": f.get("path"),
                    "name": f.get("name"),
                    "size_bytes": int(f.get("size", 0)),
                    "download_url": f.get("download_url"),
                })

            # Select bounded relevant subset for standalone ingestion
            candidates = []
            for f in remote_manifest:
                path = f.get("path") or ""
                ext = Path(path).suffix.lower()
                size_mb = float(f.get("size_bytes", 0)) / 1e6
                if ext not in ALLOWED_EXT:
                    continue
                if size_mb > MAX_FILE_MB:
                    continue
                score = relevance_score(path)
                if score <= 0 and ext not in {".md", ".yml", ".yaml"}:
                    continue
                candidates.append((score, f))

            # Highest relevance first; deterministic by path tie-breaker
            candidates.sort(key=lambda x: (-x[0], x[1].get("path", "")))
            selected = [c[1] for c in candidates[:MAX_DOWNLOAD_FILES]]

            for rec in selected:
                dl = rec.get("download_url")
                rpath = rec.get("path")
                if not dl or not rpath:
                    continue
                target = TARGET_DIR / rpath
                try:
                    download_file(dl, target)
                    downloaded.append(rpath)
                except Exception:
                    continue

            print_status(f"Remote files discovered: {len(remote_manifest)}", "INFO")
            print_status(f"Downloaded files: {len(downloaded)}", "INFO")
        except Exception as e:
            remote_error = str(e)
            print_status(f"Remote discovery failed: {e}", "WARNING")

    local_summary = summarize_local_dataset(TARGET_DIR)

    output = {
        "step": STEP_NUM,
        "status": "success" if local_summary["n_files"] > 0 else "partial",
        "description": "Standalone TDCOSMO2025 public ingestion (local-first, API fallback)",
        "source": {
            "repository": REPO,
            "api_root": API_ROOT,
            "source_mode": source_mode,
            "remote_error": remote_error,
            "max_file_mb": MAX_FILE_MB,
            "max_download_files": MAX_DOWNLOAD_FILES,
        },
        "downloaded_files": downloaded,
        "downloaded_count": len(downloaded),
        "remote_manifest_count": len(remote_manifest),
        "remote_manifest_sample": remote_manifest[:120],
        "local_summary": local_summary,
        "next_action": (
            "If lens/posterior coverage is incomplete, stage full TDCOSMO2025_public repo under "
            "data/interim/external/tdcosmo2025_public and rerun step_19."
        ),
    }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2, default=safe_json_default)

    print_status(f"Results saved to {OUTPUT_JSON}", "INFO")
    print_status(
        f"Local files: {local_summary['n_files']} | Posterior-like: {local_summary['posterior_like_count']} | "
        f"Lens coverage: {local_summary['lens_coverage_count']}/{len(LENS_ALIASES)}",
        "INFO",
    )
    if local_summary["n_files"] == 0:
        print_status("No local files available after standalone attempt.", "WARNING")
    print_status(f"Step {STEP_NUM} complete.", "INFO")


if __name__ == "__main__":
    main()
