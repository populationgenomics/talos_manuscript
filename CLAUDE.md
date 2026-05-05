# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository purpose

This repository accompanies the Talos manuscript. It contains the evaluation script that benchmarks Talos against Exomiser on cohort truth sets, plus the notebooks and source data used to generate the manuscript figures. There is no library code here — every script/notebook is a leaf that consumes data and produces a tabular result or a PDF panel.

## Layout

- `Analyses/evaluation/talos_evaluation.py` — CLI that scores Talos (and Exomiser) against per-cohort gold-standard TSVs. Cohort → input-path mapping is hard-coded in `COHORT_CONFIG` near the top of the file. Inputs are pulled from `gs://cpg-*` buckets via `cloudpathlib.AnyPath`, so running it requires GCS auth and Talos installed (see below).
- `Analyses/Jupyter_notebooks/` — `Figure2.ipynb`, `Figure3_panelA.ipynb`, `Figure3_panels_BC.ipynb`. Each reads CSVs from `../../Data/{Fig2,Fig3}/` and writes PDFs to `../../Figures/{Fig2,Fig3}/`. The relative paths assume the notebook is launched from this directory.
- `Analyses/Jupyter_notebooks/pre-submission/` and `Data/pre-submission/`, `Figures/pre-submission/` — superseded artifacts kept for provenance. Don't edit; current figure work lives in the non-`pre-submission` siblings.
- `Data/` — input CSVs for the notebooks. The Fig3 CSV filenames embed a date (e.g. `Talos_solves_VCGS-prospective_261128.csv`) — when the data refreshes, both the file in `Data/Fig3/` and the `pd.read_csv(...)` path in the notebook need to change together.
- `Figures/Fig3/Fig3.afdesign` — Affinity Designer master that composes the per-panel PDFs into the final `Fig3.pdf`. Panel D originates from `Fig3_panel_D.pptx` (PowerPoint), exported to `Fig3_panel_D.pdf`. Re-running notebooks regenerates the per-panel PDFs but does **not** rebuild `Fig3.pdf`; that step is manual in Affinity Designer.

## Running things

### Notebooks (figure generation)

```bash
cd Analyses/Jupyter_notebooks
source venv/bin/activate           # venv is checked-out-but-gitignored
pip install -r requirements.txt    # only needed for a fresh venv
jupyter lab Figure2.ipynb          # or Figure3_panelA.ipynb / Figure3_panels_BC.ipynb
```

Run notebooks from `Analyses/Jupyter_notebooks/` — they use `../../Data/...` and `../../Figures/...` paths and will silently write to or fail in the wrong place if launched elsewhere.

### Evaluation script

The script imports `from talos import models` / `from talos.utils import read_json_from_path`. Talos is not on PyPI — clone https://github.com/populationgenomics/talos and `pip install .` in that repo's root into the same venv. It also reads `gs://...` paths, so `gcloud auth application-default login` (or equivalent) is required.

Canonical invocations are listed in the docstring at the top of `talos_evaluation.py`. The general shape:

```bash
python talos_evaluation.py core --cohort acute-care \
    --summary_tsv outputs/<date>_acute-care_all_talos_evaluation.tsv \
    > outputs/<date>_acute-care_all_talos_evaluation.txt
```

Positional `core` vs `all` selects `VARIANT_TYPES_BASIC` vs `VARIANT_TYPES_ALL`. `--cohort` keys into `COHORT_CONFIG` (`acute-care`, `acute-care-singletons`, `RGP`, `RGP-singletons`). `--process_full_trios_only` restricts to trios; the manuscript reports both the trio-only and all-family numbers, so most cohorts are run twice.

## Conventions worth knowing

- **Variant counts use unique `CHR-POS-REF-ALT` strings, not `ReportVariant` objects.** Talos splits a variant into one `ReportVariant` per annotated gene, so naive counting inflates totals. See the docstring at the top of `talos_evaluation.py` and `Family.talos_candidate_count` for the canonical pattern — preserve it when adding new metrics.
- **`pre-submission/` directories are frozen.** When updating figures or data, edit the top-level `Fig2/`, `Fig3/` siblings; do not touch `pre-submission/`.
- The `talos` import in `talos_evaluation.py` is intentionally narrow (models + one util). Per the comment near the import, the script is expected to keep working even as upstream Talos diverges, so avoid pulling in additional Talos surface area.
