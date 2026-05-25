# talos_manuscript

Source data, evaluation script, and figure-generation notebooks for the Talos manuscript.

There's no library code here — every script/notebook is a leaf that consumes data and produces a tabular result or a PDF panel.

## Layout

- `Analyses/evaluation/talos_evaluation.py` — CLI that scores Talos (and Exomiser) against per-cohort gold-standard TSVs. Cohort → input-path mapping is hard-coded in `COHORT_CONFIG`. Reads `gs://cpg-*` buckets via `cloudpathlib.AnyPath`, so running it needs GCS auth and Talos installed.
- `Analyses/Jupyter_notebooks/` — figure notebooks. Each reads CSVs from `../../Data/{Fig2,Fig3}/` and writes PDFs to `../../Figures/{Fig2,Fig3}/`. **Always launch notebooks from this directory** — the relative paths assume it.
  - `Figure2.ipynb` — Figure 2.
  - `Figure3.ipynb` — **Figure 3 panels A, B, C in one figure**. Outputs `Figures/Fig3/Fig3_panels_ABC.pdf`. Panel D is composed manually (see below).
- `Data/` — input CSVs. Fig3 filenames embed a date (e.g. `Talos_solves_VCGS-prospective_261128.csv`); when data refreshes, both the file and the `pd.read_csv(...)` path in the notebook need to change together.
- `Figures/Fig3/Fig3.afdesign` — Affinity Designer master that composes per-panel PDFs into the final `Fig3.pdf`. Re-running notebooks regenerates per-panel PDFs but does **not** rebuild `Fig3.pdf` — that step is manual.
- `*/pre-submission/` — frozen artifacts kept for provenance. Don't edit; current figure work lives in the non-`pre-submission` siblings.

## Setup

```bash
cd Analyses/Jupyter_notebooks
python -m venv venv             # if not yet created
source venv/bin/activate
pip install -r requirements.txt
```

For the evaluation script you also need Talos itself, which is not on PyPI:

```bash
git clone https://github.com/populationgenomics/talos
pip install ./talos              # into the same venv
gcloud auth application-default login   # for gs:// access
```

## Regenerating figures

### Figure 2

```bash
cd Analyses/Jupyter_notebooks
jupyter lab Figure2.ipynb
```

### Figure 3 (panels A + B + C)

The whole figure (panels A, B, C) is built by a **single matplotlib `Figure`** in cell `a04cd0f7` of `Figure3.ipynb`. Top-level layout: a 2×2 `GridSpec` with panel A spanning the top row and panels B/C side-by-side in the bottom row. Panels B and C are post-shifted vertically via `set_position` so they sit higher than the gridspec would naturally place them, leaving room beneath for panel D.

```bash
cd Analyses/Jupyter_notebooks
jupyter lab Figure3.ipynb
# or render headless:
python -c "
import json, os
nb = json.load(open('Figure3.ipynb'))
import matplotlib; matplotlib.use('Agg')
exec(''.join(nb['cells'][0]['source']).replace('plt.show()', ''))
"
```

Output: `Figures/Fig3/Fig3_panels_ABC.pdf`.

To tune the figure's layout (sizes, margins, panel shifts, label positions, legend nudges, etc.) see [HANDOFF.md](HANDOFF.md) — every constant in cell `a04cd0f7` is listed there with what it controls and how it interacts with neighbouring knobs.

### Figure 3 panel D

Panel D originates from `Figures/Fig3/Fig3_panel_D.pptx` (PowerPoint), exported to `Fig3_panel_D.pdf`. The matplotlib figure draws a `d` placeholder label at the position where panel D will sit; the actual content is dropped in beneath that label in Affinity Designer.

### Final composed Fig3.pdf

Open `Figures/Fig3/Fig3.afdesign` in Affinity Designer. It links the per-panel PDFs (`Fig3_panels_ABC.pdf` + `Fig3_panel_D.pdf`); re-export the master after regenerating either source. **This step is manual** — re-running notebooks does not rebuild `Fig3.pdf`.

## Running the evaluation script

`Analyses/evaluation/talos_evaluation.py` benchmarks Talos against Exomiser on cohort truth sets. Canonical invocations are in the docstring at the top of the file. The general shape:

```bash
python talos_evaluation.py core --cohort acute-care \
    --summary_tsv outputs/<date>_acute-care_all_talos_evaluation.tsv \
    > outputs/<date>_acute-care_all_talos_evaluation.txt
```

- Positional `core` vs `all` selects `VARIANT_TYPES_BASIC` vs `VARIANT_TYPES_ALL`.
- `--cohort` keys into `COHORT_CONFIG` (`acute-care`, `acute-care-singletons`, `RGP`, `RGP-singletons`).
- `--process_full_trios_only` restricts to trios. The manuscript reports both trio-only and all-family numbers, so most cohorts are run twice.

## Conventions worth knowing

- **Variant counts use unique `CHR-POS-REF-ALT` strings, not `ReportVariant` objects.** Talos splits a variant into one `ReportVariant` per annotated gene, so naive counting inflates totals. See the docstring at the top of `talos_evaluation.py` and `Family.talos_candidate_count` for the canonical pattern.
- **`pre-submission/` directories are frozen.** When updating figures or data, edit the top-level `Fig2/`, `Fig3/` siblings.
- The `talos` import in `talos_evaluation.py` is intentionally narrow (models + one util). The script is expected to keep working even as upstream Talos diverges, so avoid pulling in additional Talos surface area.
