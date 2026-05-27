# TEP-LENS Data

Input data for the strong-lensing blind-prediction residual analysis (Paper 19).

## Structure

- `raw/` — SN Refsdal, SN H0pe, TDCOSMO/H0LiCOW catalogs, lens models
- `interim/` — Pipeline intermediate JSON/CSV
- `cosmograil/` — CosmoGRAIL temporal-shear inputs (when used)

## Provenance

See `DATA_PROVENANCE.md` for source papers, URLs, and ingestion steps.

## Reproduction

```bash
python scripts/steps/run_all_steps.py
```

Individual fetch steps (e.g. `step_01_fetch_snh0pe_data.py`) download public data where network access is available.
