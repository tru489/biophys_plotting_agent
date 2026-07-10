---
name: plot-experiment
description: >-
  Use when the user wants to plot compiled biophysics experiment data (FXM / SMR / Coulter —
  buoyant mass, density, cell volume) into the standard ridge / box / timecourse / scatter figure
  grid. Triggers on "plot my experiment", "make the ifxm/coulter figures", "generate the plots for
  <exp>", or a directory containing `<exp>_data/ifxm/experiment_data.h5` or
  `<exp>_data/coulter/metadata.csv`. Inspects the experiment's annotation schema (conditions,
  drugs, time points), generates a short plotting driver from the bundled toolkit, runs it, and
  shows the figures for fine-tuning.
---

# Plot a compiled biophysics experiment

Turn one experiment's **compiled + annotated** data into the standard figure grid. You generate a
short driver on top of the bundled `biophys_plot_toolkit.py`; the user fine-tunes it afterward in
normal Claude Code. This skill plots — it does **not** re-run the heavy analysis pipeline.

Bundled files (reference via `${CLAUDE_PLUGIN_ROOT}/skills/plot-experiment/`):
- `biophys_plot_toolkit.py` — the plotting library (loaders, ridge/box/timecourse/scatter, pptx).
- `reference_driver.py` — the driver template you adapt.
- `references/data_schema.md` — the exact h5/metadata schema. **Read it before mapping a schema.**

## Procedure

### 1. Locate & inspect the input
Find the experiment's data root: a `<exp>_data/` dir with `coulter/{metadata.csv,data.h5}` and/or
`ifxm/experiment_data.h5`. If the user only has a raw `*_compiled/` dir, that supplies the `ifxm/`
half (`experiment_data.h5`); offer to assemble the `<exp>_data/` layout (symlink/copy `ifxm/` and
point `coulter/` at an `annotate_coulter_samples.py` output). An experiment may have only one half.

Read `coulter/metadata.csv` and/or the iFXM `/metadata` to discover the **actual** annotations —
do not assume the FL5 reference set. Use a quick Python/pandas read:
```python
import pandas as pd
pd.read_csv(r"<exp>_data/coulter/metadata.csv")
pd.read_hdf(r"<exp>_data/ifxm/experiment_data.h5", "/metadata")
```

### 2. Map the schema (report back to the user)
**Annotation column names and values vary between experiments** — do not assume the reference
schema. Inspect the actual columns and identify which one plays each role:
- **sample name** (required), **h5 key** (falls back to sample name if absent),
- **condition** grouping column, **time** column (decimal hours), **drug** column (optional),
- iFXM **gate** columns (bm lower/upper, ifxm lower/upper).

The loaders take these as keyword arguments (`sample_col`, `key_col`, `condition_col`,
`time_col`, `drug_col`, `bm_lower_col`, `bm_upper_col`, `ifxm_lower_col`, `ifxm_upper_col`), all
defaulting to the reference names. If this experiment uses different names, note the mapping so
you can pass it in step 3. The toolkit degrades gracefully on its own: a **missing condition** →
one unnamed group; **missing time** → 0.0; **missing drug** → no drug figures; **missing/NaN
gate** → treated as no cutoff (not "drop everything"). Unknown **condition/drug values** are each
auto-assigned a distinct color (via `build_color_map`), so arbitrary annotation values plot fine.

State what you found: distinct condition values, whether a drug arm exists, the set of
timepoints, replicate count, and which property families are available (coulter `volume`; iFXM
`mass` / `density` / `vol_cal` / `vol_uncal`). Flag anything missing (no gates, no
`volume_calibrated`, samples that will be skipped).

**Ask the user for `baseline_density`** (g/mL) for this experiment — it is not stored in any data
file and is required to compute absolute density. (The FL5 reference experiments used `1.008`.)
If the user has no density plots / no iFXM data, this can be skipped.

### 3. Generate the driver
1. Create the analysis output dir (default: alongside the data, or ask). Into it, **copy**
   `${CLAUDE_PLUGIN_ROOT}/skills/plot-experiment/biophys_plot_toolkit.py`. Copying (rather than
   importing from the plugin cache) keeps each analysis self-contained, git-committable, and
   reproducible independent of the plugin install.
2. Adapt `reference_driver.py` into that dir: set `DATA_DIR`, `FIG_DIR`, `PPTX_OUT`,
   `BASELINE_DENSITY`, and trim `COULTER_PROPS` / `IFXM_PROPS` / `SCATTER_PAIRS` to what the
   metadata actually contains. **If the experiment's annotation columns differ from the defaults,
   pass the mapping to the loaders** (e.g.
   `tk.load_ifxm(DATA_DIR, baseline_density=BD, condition_col="treatment", time_col="t_hours")`).
   Keep `import biophys_plot_toolkit as tk` — do **not** inline the helpers. Use `scatter_2d` with
   records from `load_ifxm_paired` (row-aligned) for any property-vs-property scatter (e.g.
   mass-vs-density, volume-vs-mass); it renders a per-condition grid of per-sample panels, each a
   scatter with marginal histograms on both axes. Pass `trim_y="mad"` (or `trim_x`) to tame
   heavy-tailed axes like density.

### 4. Run it
Ensure the deps are available (numpy, pandas, matplotlib, tables/PyTables, python-pptx, Pillow).
A conda env spec is bundled at `${CLAUDE_PLUGIN_ROOT}/environment.yaml`
(`conda env create -f ...` then activate `biophys_plotting`), or reuse an existing analysis env.
Run the driver. It writes PNGs
to `<exp>_fig/` following `{datatype}_{metric}_{plottype}[_{cond|drug}].png` (and
`{datatype}_{propY}_vs_{propX}_scatter[_{cond}].png`) plus `<exp>_figures.pptx`.

### 5. Show & hand off
Surface the generated figures and the driver path. Tell the user they can now fine-tune the driver
directly in Claude Code (palettes via `COND_COLORS`/`DRUG_COLORS`, ridge bins/overlap, which
scatter pairs, figure sizes, etc.) and re-run it — the copied toolkit makes it fully editable.

## Gotchas
- **`baseline_density` has no default** — `load_ifxm` raises without it. Always set it deliberately.
- **iFXM gating**: `mass` uses `bm_gate`; `density`/`vol_cal`/`vol_uncal` share one mask on the
  *uncalibrated* volume. `load_ifxm` applies per-property outlier rejection; `load_ifxm_paired`
  keeps arrays row-aligned (use it for `scatter_2d`).
- **Starved samples** typically share `time_h == 0`; the toolkit labels them by full `sample_name`
  so replicates are distinguishable — expected, not a bug.
- **PyTables (`tables`) is required** to read the pandas HDF5 stores.
- **Missing drug arm** (e.g. a wt experiment): `drug_split` is a no-op — no drug figures, fine.
- Two HDF5 conventions coexist: pandas `HDFStore` for DataFrames vs raw `h5py` for `images.h5`
  (images are not used by these plots).
