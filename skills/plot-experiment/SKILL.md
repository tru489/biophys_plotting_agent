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
  combinators (`plot_grouped`/`compare_groups`/`timecourse_by`/`scatter_by`/`facet`/`cross_groups`/
  `grid_heatmap`)
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

**Grid-search heatmaps — suggest, never assume.** If two annotation columns were varied *together*
across samples (a grid search, e.g. media osmolarity × drug dose), a heatmap of the per-sample mean
over the two-parameter grid (`tk.grid_heatmap`) is often the clearest summary. `render_plan` lists
any such pairs it detects (`tk.suggest_grids`: both columns boolean/categorical/ordered, ≥3 levels
on one axis, genuinely crossed rather than 1:1, ≥75% of combos present) under **"Suggested grid
heatmaps — NOT included"**. Do not add them on your own:
- If a suggestion appears **and** you judge it meaningful given the metadata (the columns really
  are two independently varied parameters, not e.g. a sample ID crossed with a batch), ask the user
  with `AskUserQuestion` whether to add the heatmap, naming the two columns and the grid shape.
  Mention any repeated combination (usually a control re-run through the session): by default the
  grid cell shows the mean of those repeats (marked "mean of n=k"); offer `show_repeats=True` to
  also show each repeat individually in a side panel on the same color scale (to see drift).
- Only on a yes, put the pair in the driver's `GRID_PAIRS` (`[(x_col, y_col)]`, or
  `[{"cols": (x_col, y_col), "show_repeats": True}]`). With no suggestion, or a no, leave it empty.
- Pairs that the detector misses can still be added if the user asks for a grid heatmap.

Behavior to convey: **every** boolean/categorical/approved-ordered column becomes its own grouping
axis (no cardinality cap — everything is plotted; reorganize on a later pass). Multiple columns are
handled **independently** by default; cross-products are available on request via `cross_groups`. A
time column orders samples sequentially and drives timecourses. Gate columns still apply; missing
gate → no cutoff; unpaired samples without a volume reading get an empty `vol`.

**Unpaired runs are supported.** A paired sample uses its matched `pair_` block (unchanged); a
**mass-only** run (`mass` only) or **volume-only** run (`vol` only) falls back to the standalone
`mass_`/`vol_` blocks, so those samples still plot. `density` and `scatter_by` need pairing, so they
simply don't appear for unpaired samples — `build_plan` omits plots for absent properties
automatically.

**There is exactly one paired format** — straight from a `*_CELLGROUPED.hdf5`'s
`analysis/density/cells` table (current SMRFXMAnalysis output). Every paired sample carries all
three: `pair_mass_pg` (pg), `pair_volume_fl` (already-calibrated fL), and
`pair_cell_density_g_per_mL` (ABSOLUTE density, computed by the hdf5 pipeline itself from its own
`media_density_g_per_mL` — no baseline to supply). A pair_ block missing any of the three raises
loudly rather than silently falling back.

### 3. Generate the driver
1. Create the analysis output dir. Into it, **copy**
   `${CLAUDE_PLUGIN_ROOT}/skills/plot-experiment/biophys_plot_toolkit.py` (keeps each analysis
   self-contained, git-committable, reproducible independent of the plugin install).
2. Adapt `reference_driver.py` into that dir: set `EXP_NAME`, `COMPILED_DIR` (iFXM) and/or
   `COULTER_DIR` (Coulter, `None` if absent), `FIG_DIR`, `PPTX_OUT`, and set
   `ROLE_OVERRIDES` to the choices the user made in step 2 (resolving every `[CONFIRM]`) and
   `GRID_PAIRS` to any grid heatmaps the user confirmed (else leave it `[]`).
   **`EXP_NAME` must be the full analysis dir name including its date prefix** (e.g.
   `2026-09-22_fl5_wnki_conc-curves` → `2026-09-22_fl5_wnki_conc-curves_figures.pptx`) — never a
   shortened, dateless form. The template derives it from the driver's own dir
   (`Path(__file__).resolve().parent.name`); keep that when the driver lives in the dated analysis
   dir, otherwise hardcode the dated name (prefix the experiment date if the dir has none).
   `save_pptx` warns if the deck name has no date prefix — treat that warning as a bug. The
   template's default path is `infer_roles → build_plan → render_plan(print) → autoplot`; keep it
   for the standard grid, or drive the **explicit combinators** for full control:
   - `plot_grouped(group_col=…)` — per-group detail (ridge + box).
   - `compare_groups(group_col=…)` — cross-group comparison; default `agg="per_sample"`, both box
     and ridge.
   - `timecourse_by(time_col=…, series_col=…)` — unit-aware, sequentially ordered.
   - `scatter_by(prop_x, prop_y, group_col=…)` — per-cell scatter + marginals; **pass
     `load_ifxm_paired` records** (row-aligned).
   - `cross_groups(cols=(a, b))` — cross-product comparison, on request.
   - `grid_heatmap(cols=(x, y), show_repeats=False)` — grid-search heatmap of per-sample means
     (`agg="mean"`; repeated cells averaged), **only when the user confirmed it** (see step 2).
   - `facet(facet_col=…)` — compact one-figure grid.
   Keep `import biophys_plot_toolkit as tk` — do **not** inline helpers. No statistical outlier
   rejection is applied; tame a heavy tail with axis limits in the driver.

**Axis-label convention (default, automatic — general principle, not just this dataset):**
never repeat on a per-sample tick something the figure already states via its title or its bold
group separators. Every box plot's per-sample x-tick shows only the annotation columns *not*
otherwise conveyed elsewhere on that same figure, not the raw `sample_name`. Concretely, in
`plot_grouped`/`compare_groups`/`cross_groups`:
   - `plot_grouped(group_col=…)`: the title already fixes `group_col`, and (when a time column
     exists) the bold separators below the axis already show time — so each tick shows only the
     *other* remaining annotation column(s).
   - `compare_groups(group_col=…)`: the bold separators already show `group_col`'s value — so
     each tick shows time plus every other remaining column, but not `group_col` again.
   - `cross_groups(cols=(a, b))`: the bold separators already show the joined `a | b` key — so
     each tick shows only whatever's left over (typically just time), not `a` or `b` again.
   This is implemented once, generically, in `biophys_plot_toolkit.py` (`_leftover_label`, wired
   into the box-plot path of all three combinators via each function's `roles` dict — it has no
   per-experiment hardcoding). Ridge plots and pooled (`agg="pool"`) comparisons are intentionally
   left showing their existing fuller/group-only labels. Apply this same "don't restate what's
   already on the plot" principle to any new combinator or hand-rolled plotting code you add here
   — it isn't specific to time/osm/drug columns, it's a general labeling rule for every axis.

### 4. Run it
Ensure the deps are available (numpy, pandas, matplotlib, openpyxl, python-pptx, Pillow).
A conda env spec is bundled at `${CLAUDE_PLUGIN_ROOT}/environment.yaml`
(`conda env create -f ...` then activate `biophys_plotting`), or reuse an existing analysis env.
Run the driver. It writes PNGs to `<exp>_fig/` (see the naming grid in `data_schema.md`:
`{datatype}_{metric}_{plottype}_{col}={value}.png`, `..._by_{col}.png`,
`{datatype}_{propY}_vs_{propX}[_{col}={value}].png`, `{datatype}_{prop}_heatmap_{colX}-x-{colY}.png`)
plus `<exp>_figures.pptx` (with `<exp>` the dated analysis dir name).

**Point counts (automatic):** every ridge row shows its `n=` (number of points plotted) in small
gray text just outside the right edge of the axes, and every box shows its `n=` just above its
highest plotted point (with headroom added to the y-axis) — both via `draw_ridge`/`draw_boxes`, so
all combinators get them. Keep this in hand-rolled ridge/box plots (use those primitives).

**Legend placement (automatic — a legend must never sit on top of plot features):** every legend
the toolkit draws goes through `tk.place_legend(ax, ...)`. It renders the figure without the
legend, and if the legend's box would cover any drawn feature (lines, boxes/whiskers, ridge fills,
jittered points, group separators) it tries the other in-axes anchor points (upper right → upper
left → lower right → … → center), and if none is clear it moves the legend just outside the axes on
the right (`bbox_inches="tight"` keeps it in the PNG). `tk._save` re-checks every legend right
before writing, so a legend added by a driver (`ax.legend(...)`) or drawn before the data is also
fixed. **In any hand-rolled plotting code you write, create legends with
`tk.place_legend(ax, handles=..., loc=<preferred>)` after drawing the data, and write the figure
with `tk._save(fig, name, FIG_DIR)`** — never a bare `ax.legend(loc="best")` +
`fig.savefig(...)`, since `"best"` does not reliably avoid boxes, fills or scatter clouds. If you
must use `fig.savefig` directly, call `tk.place_legend(ax)` immediately before it.

### 5. Show & hand off
**Look at the figures before handing off** (open the PNGs) and confirm no legend overlaps a plot
feature — if one does (e.g. a hand-rolled legend, or the printed `WARNING: legend overlaps plotted
features at every position tried`), fix it (shrink/trim the legend, enlarge the figure, or route it
through `tk.place_legend`) and re-run rather than shipping it.

Surface the generated figures and the driver path. Tell the user they can fine-tune the driver
directly in Claude Code (`ROLE_OVERRIDES`, which columns to group/compare/cross, palettes via
`COND_COLORS`/`DRUG_COLORS`/`BOOL_COLORS`, ridge bins/overlap, scatter pairs, figure sizes, and
box-tick labeling via `compare_groups`/`cross_groups`'s `label_exclude=` if the automatic
"don't restate the title/bold-separator" rule above needs a different set excluded for some
column) and re-run — the copied toolkit makes it fully editable.

## Outlier rejection (opt-in — the loaders never trim data)

By default **no statistical outlier rejection** is applied — the data is loaded verbatim so the user
can decide per experiment / per sample. The toolkit provides an opt-in transform, `reject_outliers`,
that returns cleaned records feeding straight into `build_plan`/`autoplot`/any combinator, so one
call cleans every downstream plot.

**When the user brings up outlier rejection at all** (e.g. "add outlier rejection", "trim the density
tails", "reject outliers on volume"), do **not** guess — present a menu with `AskUserQuestion`
enumerating the full spec, then wire the answer into the driver. Ask for:

1. **Which properties** to clean — any of coulter `volume`; iFXM `mass` / `density` / `vol`
   (multi-select; can differ per property).
2. **Method** (per property allowed): `mad` (modified z-score, robust — best for density),
   `iqr` (Tukey fences, robust, matches the box whiskers), `percentile` (fixed-fraction tail clip).
3. **Scope**: `per_sample` (each sample trimmed on its own stats — default) or `pooled` (one
   global cutoff across all cells).
4. **Log-space?** (`log=True`) — for log-normal mass/volume so the high tail isn't over-trimmed.
5. **Method params** if they care — `mad thresh` (3.5), `iqr k` (1.5), `percentile lower/upper`
   (1/99).

Then add lines to the driver (and re-run):
```python
ifxm        = tk.reject_outliers(ifxm, method={"density": "mad", "mass": "iqr"}, scope="per_sample")
ifxm_paired = tk.reject_outliers(ifxm_paired, method="mad", props=["density"], paired=True)  # scatter
coulter     = tk.reject_outliers(coulter, method="iqr", props=["volume"])
```
Notes: apply to the **scatter** records (`ifxm_paired`) with `paired=True` so a cell's props stay
row-aligned; for the distribution records use the plain call. `verbose=True` prints how many cells
each sample dropped. Low-level keep-masks (`tk.outlier_mask`, `keep_mad`/`keep_iqr`/`keep_percentile`)
are exposed for bespoke logic. (k-sigma/3-std is intentionally not built in — ask if the user wants it.)

## Gotchas
- **iFXM gating**: `mass` uses `bm_gate`; `density`/`vol` share one mask on `pair_volume_fl` (the
  paired sample's volume). No statistical outlier rejection is applied by the loaders — only
  non-finite values are dropped (see the opt-in `reject_outliers` above for trimming). `load_ifxm`
  gates mass and the volume props with separate masks (so per-property arrays can differ in
  length); `load_ifxm_paired` uses one shared mask to keep arrays row-aligned — always use it for
  `scatter_by` (and pass `paired_records=` to `autoplot`), or a scatter's x/y won't pair.
- **A pair_ block must carry all three columns** — `pair_mass_pg`, `pair_volume_fl`,
  `pair_cell_density_g_per_mL`. `load_ifxm`/`load_ifxm_paired` raise a clear error naming the
  sample if any are missing, rather than silently dropping data. This is the only paired format.
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
- **Legends are auto-placed, not `loc="best"`** — see "Legend placement" in step 3. Needs
  matplotlib ≥ 3.8 (`Legend.set_loc`); the bundled `environment.yaml` resolves to that or newer.
- **Box-plot tick labels are auto-trimmed to leftover info only** (see the axis-label convention
  in step 3) — if a customized/hand-rolled box plot looks like it's duplicating its own title or
  bold separators on every tick, that's a bug to fix the same way, not a one-off label edit.
