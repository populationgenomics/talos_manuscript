# Fig 3 panel A+B+C — handoff

## Where we are

**Branch:** `R3_tweaks`.
**Working file:** [`Analyses/Jupyter_notebooks/Figure3_panelA.ipynb`](Analyses/Jupyter_notebooks/Figure3_panelA.ipynb), cell `a04cd0f7`.
**Output PDF:** [`Figures/Fig3/Fig3_panels_ABC.pdf`](Figures/Fig3/Fig3_panels_ABC.pdf) (uncommitted).

The single cell now builds **panels A, B, C in one figure**. Panel D is composed manually in Affinity Designer beneath a `d` placeholder label that the cell already draws.

## Render & inspect

```bash
cd /Users/cas/dev/talos_manuscript-clean
/Users/cas/dev/talos_manuscript-clean/Analyses/Jupyter_notebooks/venv/bin/python -c "
import json, os
nb = json.load(open('Analyses/Jupyter_notebooks/Figure3_panelA.ipynb'))
src = ''.join(nb['cells'][0]['source'])
os.chdir('Analyses/Jupyter_notebooks')
import matplotlib; matplotlib.use('Agg')
exec(src.replace('plt.show()', ''))
"
pdftoppm -png -r 180 Figures/Fig3/Fig3_panels_ABC.pdf /tmp/fig3_check
# then Read /tmp/fig3_check-1.png
```

To probe actual rendered geometry (useful when aligning labels or estimating tick widths), exec the cell source then introspect via `ax.get_position()`, `lbl.get_window_extent(renderer)`, etc. Example:
```python
print(MARGIN_LEFT * fig.bbox.width, 'px = MARGIN_LEFT')
for lbl in ax_b.get_yticklabels():
    print(lbl.get_text(), lbl.get_window_extent(renderer))
```

## Layout architecture

Top-level: `gs_top = GridSpec(2, 2)` with row 0 = panel A spanning both cols, row 1 = panel B (left) + panel C (right).

- **Panel A** = `gs_a = GridSpecFromSubplotSpec(6, 6, subplot_spec=gs_top[0, :])`. Cohort cols are sized from data via `max_cols_for_cohort()` so cells are uniform. The vertical legend on the right is drawn by **measure-and-stack** (each mini-legend is rendered, measured, then re-anchored top-down in figure coords). Don't break the measure-and-stack — it's robust to wrap-width changes.
- **Panel B** lives inside `gs_top[1, 0]` but is itself wrapped in a 2×2 `gs_b`:
  - `gs_b[0, 0]` = empty left padding (lets y-tick labels render INSIDE the figure margin and align with panel A's row labels)
  - `gs_b[0, 1]` = the actual step plot (`ax_b`)
  - `gs_b[1, *]` = empty bottom padding (room for rotated date tick labels and the x-axis title)
- **Panel C** = `gs_c = GridSpecFromSubplotSpec(2, N_YEARS+1, subplot_spec=gs_top[1, 1])` — row 0 waffles, row 1 year labels, legend col on the right.
- **Panels A, B, C are then post-shifted** vertically via `set_position` (see *Vertical shifts* below) so B and C ride higher than the gridspec would naturally place them.

### pywaffle SubplotSpec hack (panel A)

pywaffle's `plots` dict only accepts `int`/`str`/`tuple` keys, then calls `add_subplot(*loc)`. We pass each `SubplotSpec` wrapped in a single-element tuple so it unpacks cleanly. Comment lives next to the dict construction.

### pywaffle resets the anchor (panel C)

`Waffle.make_waffle(ax=…)` overrides any previously-set `ax.set_anchor('S')`. The fix in the loop is **`set_anchor('S')` AFTER `make_waffle`** — that's why bottom-alignment works now. If you switch to a different waffle library, double-check this.

## Knobs you can tune

All at the top of cell `a04cd0f7`. mm units used wherever a knob has a physical meaning so you can think in the units you'll see on screen.

### Outer canvas

- `FIG_WIDTH_MM`, `FIG_HEIGHT_MM` — total figure size (180×155 currently).
- `MARGIN_LEFT`, `MARGIN_RIGHT`, `MARGIN_TOP`, `MARGIN_BOTTOM` — outer margins as fig-fraction.
  - `MARGIN_TOP=0.95` reserves ~7.7 mm at the top for panel labels (a/b/c) to float above their panels.
  - `MARGIN_RIGHT=0.97` keeps panel A's vertical legend from clipping. Bumping toward `0.99` will reclaim white space on the right but risks clipping the legend.
  - `MARGIN_BOTTOM=0.09` leaves room for panel B's rotated date labels and x-axis title, plus the manual panel D below.

### Top-level row/column proportions

- `PANEL_A_H`, `PANEL_BC_H` — relative heights of the panel A row and the B+C row.
- `PANEL_B_W`, `PANEL_C_W` — relative widths of panel B and panel C in the BC row.
- `TOP_HSPACE`, `TOP_WSPACE` — gaps (matplotlib's "fraction of avg cell" convention).

### Panel B internal padding

- `B_LEFT_PAD` (currently `0.11`) — width of `gs_b[0, 0]` relative to the plot col. Sized so the right-anchored y-tick labels (e.g. "4000") have their **left edge** land at `MARGIN_LEFT`, aligned with panel A's row labels. To re-tune, render and check `ax_b.get_yticklabels()[N].get_window_extent(renderer)` against `MARGIN_LEFT * fig.bbox.width`.
- `B_BOT_PAD` (currently `0.45`) — height of `gs_b[1, *]` relative to the plot row. Bigger value → smaller plot, more room for the rotated x-tick labels and "Date of entry into reanalysis study" axis label below.

### Panel C waffle config

- `C_COLUMNS` (currently `5`) — tiles wide per yearly bar. Bigger = shorter, wider bars.
- `gs_c` `width_ratios=[1]*N_YEARS + [1.6]` — bump the legend col ratio if labels need more horizontal room.
- `set_anchor('S')` after `make_waffle` bottom-aligns all 4 yearly waffles to a common baseline.

### Panel labels (a / b / c / d)

- `PANEL_LABEL_X = 0.005` — shared x for `a`, `b`, `d`. Move toward 0 for further-left labels; never above `MARGIN_LEFT` or they'll move inside the panel.
- `PANEL_LABEL_Y_OFFSET = 5/FIG_HEIGHT_MM` — labels float ~5 mm above each panel's top. Reduce for tighter packing.
- `PANEL_LABEL_FONTSIZE = 9` — bold sans serif.
- **`c` uses `y_shift=B_SHIFT`** (NOT `C_SHIFT`) — this is intentional, so `c` stays horizontally aligned with `b` even when panel C's content is shifted independently.
- `PANEL_D_LABEL_Y` — absolute y for the manual-panel-D placeholder label. Adjust to wherever Affinity Designer is going to drop panel D in the final composition.

### Per-panel legend nudges (mm)

- `A_LEGEND_Y_NUDGE_MM` (currently `+4`) — lifts the entire panel A legend stack relative to the waffles. The measure-and-stack code starts from `a_leg_pos.y1 + nudge / FIG_HEIGHT_MM` and walks downward.
- `C_LEGEND_Y_NUDGE_MM` (currently `-4`) — drops panel C's legend relative to its cell. Applied AFTER the panel C shift via `leg_c.set_bbox_to_anchor(..., transform=fig.transFigure)`, so it tracks the (shifted) `ax_c_leg` position.

### Per-panel vertical shifts (mm)

After all axes/legends are built, each panel-B and panel-C axes is moved up via `set_position` by `B_SHIFT_MM` / `C_SHIFT_MM`:

- `B_SHIFT_MM` (currently `4`) — panel B plot + tick labels + axis label all move up together (the rotated tick labels and x-label are positioned relative to `ax_b`, so they follow).
- `C_SHIFT_MM` (currently `6`) — every axes in `panel_c_axes` moves up: 4 waffles, 4 year-label axes, 1 legend axes. The legend bbox is then re-anchored in figure coords so the `C_LEGEND_Y_NUDGE_MM` still applies relative to the new cell top.

If you add new panel-C axes, **append them to `panel_c_axes`** so they're included in the shift.

## Common revision recipes

### "Move panel X by N mm"

- Whole panel B: change `B_SHIFT_MM`.
- Whole panel C: change `C_SHIFT_MM`. The `c` label stays put because it uses `B_SHIFT` — change to `y_shift=C_SHIFT` if you want `c` to follow.
- Panel D (manual placeholder): change `PANEL_D_LABEL_Y`.

### "Y-tick labels (or x-tick labels) clipping the figure edge"

- For y-tick labels of panel B: increase `B_LEFT_PAD`. Verify by measuring `ax_b.get_yticklabels()[-1].get_window_extent(renderer)` against `0` (figure left edge).
- For panel A right legend: reduce `MARGIN_RIGHT`.
- Generic: reduce `LEGEND_LABEL_WRAP` so legend text wraps narrower.

### "Panel B labels overflow into panel C's height"

Bump `B_BOT_PAD` so the gs_b cell reserves more bottom row.

### "Panel A's row labels and panel B's y-tick labels don't align"

Compute the actual mm offset between them at the current `B_LEFT_PAD`, then convert to a delta:
```
delta_in_pad = (offset_px / fig.bbox.width) / cell_B_width_in_fig × (1 + current_pad)
```
…and adjust `B_LEFT_PAD` by that delta. Alternatively just bisect.

### "Bottom margin too tight / too generous"

Reduce/increase `MARGIN_BOTTOM`. Watch panel B's "Date of entry into reanalysis study" label — it's the lowest-rendered element from the cell (panel D is drawn manually).

## Gotchas

1. **The figure file is `Figure3_panelA.ipynb` but produces all of A+B+C.** Renaming to `Figure3.ipynb` and deleting the now-redundant `Figure3_panels_BC.ipynb` is on the punch list (see *Outstanding tasks* below). The sibling cell `869ea925` does `! open ../../Figures/Fig3/Fig3_panels_ABC.pdf`.
2. **CLAUDE.md still lists `Figure3_panelA.ipynb` and `Figure3_panels_BC.ipynb` as separate files.** Update when you rename.
3. **Don't add `Figure3_panels_BC.ipynb` regenerations to the new cell** — its data/figures are now produced by `Figure3_panelA.ipynb`.
4. **Pre-existing un-related changes in working tree** that should NOT be committed as part of this refactor:
   - `Figures/Fig3/Fig3_panel_D.pdf` and `.pptx` (panel D PowerPoint work, separate)
   - `Figures/Fig3/Fig3.afdesign~lock~`, `~$Fig3_panel_D.pptx` (Affinity / PowerPoint lock files)

## Outstanding tasks

1. **Rename file** `Figure3_panelA.ipynb` → `Figure3.ipynb`, update cell `869ea925` if needed, delete redundant `Figure3_panels_BC.ipynb`.
2. **Update CLAUDE.md** — the layout description currently lists the two notebooks as separate.
3. **Commit** to `R3_tweaks`. Stage only the notebook + `Fig3_panels_ABC.pdf`; leave the panel D / lock files / untracked CLAUDE.md alone.
4. **Compose final Fig3.pdf** in Affinity from `Fig3_panels_ABC.pdf` + `Fig3_panel_D.pdf`. The notebook places a `d` placeholder; the actual panel D content is added manually.

## Conversation summary so we don't go in circles

The arc was: combine A+B+C → fix B clipping (left pad + bottom pad inside `gs_b`) → bottom-align C waffles (call `set_anchor('S')` AFTER `make_waffle`) → align y-tick labels with panel A row labels (`B_LEFT_PAD = 0.11` after measuring `"4000"` width = 21 px) → fix panel A legend clipping (`MARGIN_RIGHT = 0.97`) → float panel labels 5 mm above panels (`MARGIN_TOP = 0.95` to make room) → nudge legends per panel (`A_LEGEND_Y_NUDGE_MM = +4`, `C_LEGEND_Y_NUDGE_MM = -4`) → independently shift panels B (4 mm up) and C (6 mm up) via post-render `set_position` → add `d` placeholder for the manually-composed panel D.
