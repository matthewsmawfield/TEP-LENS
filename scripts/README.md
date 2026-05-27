# TEP-LENS Scripts

## Full pipeline

```bash
python scripts/steps/run_all_steps.py
```

Runs steps `01`–`20` (SN Refsdal closure through external completeness synthesis).

## Step modules

Individual steps live in `scripts/steps/`. See `run_all_steps.py` for the canonical list.

## Site build

```bash
cd site && npm ci && npm run build:markdown
```

Generates `manuscripts/19manuscript-tep-lens.md`.
