---
name: plot-experiment
description: >-
  Use when the user wants to plot compiled biophysics experiment data (FXM / SMR / Coulter —
  buoyant mass, density, cell volume) into the standard ridge / box / timecourse / scatter figure
  grid. Triggers on "plot my experiment", "make the ifxm/coulter figures", "generate the plots for
  <exp>", or a directory containing a `*_compiled/experiment_data.xlsx` (iFXM) or a
  `*_coulter_sample_annotation/metadata.csv` (Coulter). Inspects the experiment's annotation schema (conditions,
  drugs, time points), generates a short plotting driver from the bundled toolkit, runs it, and
  shows the figures for fine-tuning.
---

# Plot a compiled biophysics experiment

Turn one experiment's **compiled + annotated** data into a figure grid. The toolkit inspects
whatever hand-added annotation columns exist, classifies each into a **role**, and uses those roles
to group / compare / order / color the data — so plots adapt to each experiment's schema instead of
a fixed condition/drug/time set. You generate a short driver on top of the bundled
`biophys_plot_toolkit.py`; the user fine-tunes it afterward. This skill plots — it does **not**
re-run the heavy analysis pipeline.

Bundled files (reference via `${CLAUDE_PLUGIN_ROOT}/skills/plot-experiment/`):
- `biophys_plot_toolkit.py` — the library: loaders → `infer_roles` → low-level `draw_*` →
  combinators (`plot_grouped`/`compare_groups`/`timecourse_by`/`scatter_by`/`facet`/`cross_groups`)
  → `build_plan`/`render_plan`/`autoplot` → pptx.
- `reference_driver.py` — the driver template you adapt.
- `references/data_schema.md` — the exact xlsx/csv/metadata schema + role table. **Read it first.**

## Procedure

### 1. Locate & inspect the input
The loaders read the **raw** biophys_helpers outputs directly (no reorg step). Find whichever the
experiment has (it may have only one):
- **iFXM** — a `*_compiled/` dir (from `compile_experiment.py`) holding `experiment_data.xlsx`.
- **Coulter** — a `*_coulter_sample_annotation/` dir (from `annotate_coulter_samples.py`) holding
  `metadata.csv` plus a single-cell data CSV (its filename is the original input CSV's name).

Read the `metadata` sheet / `metadata.csv` to discover the **actual** annotations — do not assume
the FL5 reference set. Use a quick Python/pandas read:
```python
import pandas as pd
pd.read_csv(r"<...>_coulter_sample_annotation/metadata.csv")          # Coulter
pd.read_excel(r"<...>_compiled/experiment_data.xlsx", sheet_name="metadata")  # iFXM
```

### 2. Infer roles and present the plan (ALWAYS get approval before generating)
The framework does not assume a fixed schema. It reads **whatever hand-added annotation columns
exist** and classifies each into a **role** via `tk.infer_roles(records)`:

| role | meaning | drives |
|------|---------|--------|
| `boolean` | yes/no, true/false, 0/1 (e.g. `is_activated`) | per-value plots **and** a cross-value comparison |
| `categorical` | low-ish-cardinality strings (e.g. `media`) | per-value plots + comparison + facet |
| `time` | name like `time_h`/`time_min`/`t_hours`; unit parsed → hours | **sequential ordering** + timecourse x-axis |
| `ordered` | numeric gradient (e.g. `dose_uM`, `passage`) — **flagged `[CONFIRM]`** | ordered grouping/comparison (after you confirm) |
| `continuous` | high-cardinality numeric | color/scatter axis only |
| `label` | `sample_name` / free-text identity | labels only |
| `structural` | `sheet_name`, `hdf5_key`, `has_*`, `*_gate_*`, … | ignored |

Do this: load the records, run `infer_roles`, build a plan with `tk.build_plan(...)`, and **show
the user `tk.render_plan(plan)`** — the role of every column, any `[CONFIRM]` gradient guesses, and
the list of proposed plots. **Always present this and get approval/overrides before generating the
driver** (this is a hard requirement). The user resolves `[CONFIRM]` columns and can re-map
anything via `overrides={col: "ordered"|"categorical"|...}`.

Behavior to convey: **every** boolean/categorical/approved-ordered column becomes its own grouping
axis (no cardinality cap — everything is plotted; reorganize on a later pass). Multiple columns are
handled **independently** by default; cross-products are available on request via `cross_groups`. A
time column orders samples sequentially and drives timecourses. Gate columns still apply; missing
gate → no cutoff; samples with no paired block are skipped; uncalibrated samples get an empty
`vol_cal`.

**Ask the user for `baseline_density`** (g/mL) — not stored in any data file, required for absolute
density. (FL5 reference used `1.008`.) Skip if no iFXM/density.

### 3. Generate the driver
1. Create the analysis output dir. Into it, **copy**
   `${CLAUDE_PLUGIN_ROOT}/skills/plot-experiment/biophys_plot_toolkit.py` (keeps each analysis
   self-contained, git-committable, reproducible independent of the plugin install).
2. Adapt `reference_driver.py` into that dir: set `EXP_NAME`, `COMPILED_DIR` (iFXM) and/or
   `COULTER_DIR` (Coulter, `None` if absent), `FIG_DIR`, `PPTX_OUT`, `BASELINE_DENSITY`, and set
   `ROLE_OVERRIDES` to the choices the user made in step 2 (resolving every `[CONFIRM]`). The
   template's default path is `infer_roles → build_plan → render_plan(print) → autoplot`; keep it
   for the standard grid, or drive the **explicit combinators** for full control:
   - `plot_grouped(group_col=…)` — per-group detail (ridge + box).
   - `compare_groups(group_col=…)` — cross-group comparison; default `agg="per_sample"`, both box
     and ridge.
   - `timecourse_by(time_col=…, series_col=…)` — unit-aware, sequentially ordered.
   - `scatter_by(prop_x, prop_y, group_col=…)` — per-cell scatter + marginals; **pass
     `load_ifxm_paired` records** (row-aligned).
   - `cross_groups(cols=(a, b))` — cross-product comparison, on request.
   - `facet(facet_col=…)` — compact one-figure grid.
   Keep `import biophys_plot_toolkit as tk` — do **not** inline helpers. No statistical outlier
   rejection is applied; tame a heavy tail with axis limits in the driver.

### 4. Run it
Ensure the deps are available (numpy, pandas, matplotlib, openpyxl, python-pptx, Pillow).
A conda env spec is bundled at `${CLAUDE_PLUGIN_ROOT}/environment.yaml`
(`conda env create -f ...` then activate `biophys_plotting`), or reuse an existing analysis env.
Run the driver. It writes PNGs to `<exp>_fig/` (see the naming grid in `data_schema.md`:
`{datatype}_{metric}_{plottype}_{col}={value}.png`, `..._by_{col}.png`,
`{datatype}_{propY}_vs_{propX}[_{col}={value}].png`) plus `<exp>_figures.pptx`.

### 5. Show & hand off
Surface the generated figures and the driver path. Tell the user they can fine-tune the driver
directly in Claude Code (`ROLE_OVERRIDES`, which columns to group/compare/cross, palettes via
`COND_COLORS`/`DRUG_COLORS`/`BOOL_COLORS`, ridge bins/overlap, scatter pairs, figure sizes) and
re-run — the copied toolkit makes it fully editable.

## Gotchas
- **`baseline_density` has no default** — `load_ifxm` raises without it. Always set it deliberately.
- **iFXM gating**: `mass` uses `bm_gate`; `density`/`vol_cal`/`vol_uncal` share one mask on the
  *uncalibrated* volume. No statistical outlier rejection is applied — only non-finite values are
  dropped. `load_ifxm` gates mass and the volume props with separate masks (so per-property arrays
  can differ in length); `load_ifxm_paired` uses one shared mask to keep arrays row-aligned — always
  use it for `scatter_by` (and pass `paired_records=` to `autoplot`), or a scatter's x/y won't pair.
- **Roles are inferred, not fixed** — `condition`/`time_h`/`drug_name` are just the *reference*
  column names; any hand-added column works. If a numeric column is misread (gradient vs category),
  fix it with `ROLE_OVERRIDES`/`overrides=`, not by renaming data.
- **Booleans** need values in {yes/no, true/false, 0/1}; checkbox columns from the GUI already are.
- **`openpyxl` is required** to read `experiment_data.xlsx` (Coulter needs only pandas' CSV reader).
- **No fixed-role drivers** — always use the column-parameterized combinators (`plot_grouped`,
  `compare_groups`, `timecourse_by`, `scatter_by`, `facet`, `cross_groups`) or `autoplot`; there is
  no `ridge_box_by_condition`/`drug_split`/`scatter_2d`.
- **iFXM sheets are keyed by `sheet_name`, not `sample_name`** (Excel's 31-char sanitized name);
  the loaders resolve sheets via the metadata `sheet_name` column automatically.
- `images.h5` (raw `h5py` BF image stacks, alongside `experiment_data.xlsx`) is **not** used by
  these plots.
