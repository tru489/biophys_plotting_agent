---
name: plot-experiment
description: >-
  Use when the user wants to plot compiled biophysics experiment data (FXM / SMR / Coulter —
  buoyant mass, density, cell volume) into the standard ridge / box / timecourse / scatter figure
  grid. Triggers on "plot my experiment", "make the ifxm/coulter figures", "generate the plots for
  <exp>", or a directory containing a `*_compiled/experiment_data.xlsx` (iFXM) or a
  `*_coulter_sample_annotation/metadata.csv` (Coulter). Inspects the experiment's annotation schema (conditions,
  drugs, time points), generates a short plotting driver from the bundled toolkit, runs it, and
  shows the figures for fine-tuning. Also triggers on "plot this like <past experiment>" — reads
  the past analysis's `*_analysis_notes.md` to recapitulate that analysis on a new dataset.
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
  `grid_heatmap`/`props_vs_time`)
  → `build_plan`/`render_plan`/`autoplot` → pptx (plus a default second `props_vs_time` deck).
- `reference_driver.py` — the driver template you adapt.
- `references/data_schema.md` — the exact xlsx/csv/metadata schema + role table. **Read it first.**
- `references/analysis_notes_template.md` — outline for the per-analysis notes file (step 6).

## Procedure

### 0. Reusing a past analysis (only when the user points to one)
If the user asks to analyze a new dataset "like" a past one and points at a past analysis dir (or
its dataset dir), find its `*_analysis_notes.md` (next to the `*_figures.pptx`) and read it, plus
the driver it names. Use it as the starting point for steps 1–3: same properties, grouping axes,
`ROLE_OVERRIDES`, `GRID_PAIRS`, scatter pairs, outlier rejection and customizations — mapped onto
the new dataset's columns. Still run steps 1–2 on the new data: inspect its metadata, and in the
plan you present call out every difference from the past analysis (renamed/missing/new columns,
different levels, a grid that's no longer crossed) and what you propose for each. Anything the
notes flag as experiment-specific gets re-decided with the user, not copied. If there's no notes
file, fall back to the past driver alone and say so.

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
| `structural` | `sheet_name`, `hdf5_key`, `has_*`, `gate_*`, `n_<block>_(un)gated`, legacy `*_gate_*`, … | ignored |

Do this: load the records, run `infer_roles`, build a plan with `tk.build_plan(...)`, and **show
the user `tk.render_plan(plan)`** — the role of every column, any `[CONFIRM]` gradient guesses, and
the list of proposed plots. **Always present this and get approval/overrides before generating the
driver** (this is a hard requirement). The user resolves `[CONFIRM]` columns and can re-map
anything via `overrides={col: "ordered"|"categorical"|...}`.

**Cell counts — always in the plan.** `render_plan` ends with a per-sample table: cells **before
gating**, **after gating**, and **after outlier rejection** (one column per property). In the plan,
show it for the recommended rejection so the user sees what each step removes before answering
the outlier question: `print(tk.render_cell_counts(tk.reject_outliers(ifxm, method="iqr",
props=[...], scope="per_sample")))` (and the same for `coulter`). With the default joint trim, a
paired sample's per-property columns are all equal. Point out any sample where
gating or rejection removes an unusually large fraction. The same table goes in the notes (step 6)
and on the first slide(s) of the main deck (`save_pptx(..., count_tables=...)`, already in the
template).

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

**Outlier rejection — ask every run (hard requirement).** Alongside the plan approval, always ask
with `AskUserQuestion` whether to apply outlier rejection. Propose, as the first
"(Recommended)" option, **Tukey fences (`iqr`, k=1.5), `per_sample` scope, on mass, volume and
density** (every one of iFXM `mass`/`density`/`vol` and Coulter `volume` that the experiment
has). Other options: **none** (plot verbatim), and **customize** (then walk through the full menu
in "Outlier rejection" below). When reusing a past analysis (step 0), also offer the past
analysis's spec as an option if it differs. Wire the answer into the driver (see that section)
before running it; never skip the question even if the user already approved the plot plan.

Behavior to convey: **every** boolean/categorical/approved-ordered column becomes its own grouping
axis (no cardinality cap — everything is plotted; reorganize on a later pass). Multiple columns are
handled **independently** by default; cross-products are available on request via `cross_groups`. A
time column orders samples sequentially and drives timecourses. **Gating is done upstream**:
`compile_experiment.py` writes every block already filtered to the AND of all gates in the
sample's gate file (`mass_pg`, `baseline_norm`, `bl_slope`, `node_dev_mean`, `volume_fl`), so only
gated cells reach any plot (ridge, box, timecourse, heatmap, scatter, props-vs-time) and the
plugin does not gate again. A sample with `has_gate = False` was compiled ungated — mention it in
the plan. Older xlsx files (with `bm_gate_*` / `ifxm_gate_*` columns) are still gated by the
loaders, the old way. Unpaired samples without a volume reading get an empty `vol`.

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
     `load_ifxm_paired` records** (row-aligned). All panels in one figure share x/y limits that
     span every panel's points. Keep that rule in any hand-rolled scatter grid (`tk._span_limits`
     + `draw_scatter_marginal(..., xlim=, ylim=)`): no per-panel autoscaling, no clipping.
   - `cross_groups(cols=(a, b))` — cross-product comparison, on request.
   - `grid_heatmap(cols=(x, y), show_repeats=False)` — grid-search heatmap of per-sample means
     (`agg="mean"`; repeated cells averaged), **only when the user confirmed it** (see step 2).
   - `facet(facet_col=…)` — compact one-figure grid.
   - `props_vs_time(paired_records, "ifxm", TIME_FIG_DIR, samples=[…], props=(…), trend=…)` —
     **always made, as a separate default deck** (not in `build_plan`; the template runs it
     after the main deck): one figure per paired sample of per-cell mass / vol / density vs. run
     time (stacked panels) with a black trend line, to show drift within a run. Figures go to
     `<exp>_vs_time_fig/` and the deck to `<exp>_props_vs_time.pptx` beside `<exp>_figures.pptx`
     — never into the main `FIG_DIR`/deck. **Do not ask whether to make it**; keep the
     template's block for every iFXM experiment with paired samples (skip it only if the user
     says to). Defaults: every paired sample, mass/vol/density, `trend="binned_median"`. **Pass
     `load_ifxm_paired` records** — they carry the per-cell `time_h` (matched from the mass_
     block's `peak_time_h`), which `reject_outliers` keeps aligned, so the outlier
     choice from step 2 applies here too. `trend`: `"binned_median"` (`bins=30`),
     `"rolling_mean"` (moving average) or `"rolling_median"` (`window=` cells, default ~n/30), or
     `None`. x-axis auto-switches to minutes for runs under 1 h (`time_unit=`). Narrow
     `samples=`/`props=` or change `trend` only when the user asks.
   Keep `import biophys_plot_toolkit as tk` — do **not** inline helpers. Add the
   `reject_outliers` lines the user chose in step 2 (none if they declined); beyond that, tame a
   heavy tail with axis limits in the driver.

**Box plot definition (default, automatic):** box = Q1–Q3, line = median, whiskers = 5th and
95th percentile of each sample's plotted cells (i.e. after outlier rejection, so with the default
Tukey trim they are percentiles of the Tukey-trimmed data). All kept cells are also drawn as
jittered points, so points past the caps are the outer 5% tails, not outliers. State this when
describing box plots in the analysis notes.

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
`{datatype}_{propY}_vs_{propX}[_{col}={value}].png`, `{datatype}_{prop}_heatmap_{colX}-x-{colY}.png`,
plus `<exp>_figures.pptx` (with `<exp>` the dated analysis dir name), and — for iFXM with paired
samples — the per-sample `{datatype}_{props}_vs_time_{sample}.png` figures in `<exp>_vs_time_fig/`
compiled into `<exp>_props_vs_time.pptx`, and (iFXM) the plotted cells to
`processed_data/<exp>_processed.xlsx` (see "Processed-data export" below). You then write
`<exp>_analysis_notes.md` beside the pptx (step 6).

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

### 6. Write the analysis notes (always, alongside the pptx)
Every time the skill produces a pptx, write `<EXP_NAME>_analysis_notes.md` **in the same dir as
`<EXP_NAME>_figures.pptx`**, following `references/analysis_notes_template.md`. Its purpose: a
future run pointed at this analysis (step 0) can read it and recapitulate the analysis on a new
dataset. Keep it short and factual — only the data-structure details that drive the plotting:
inputs and sample counts (paired / unpaired), the per-sample cell-count table the driver prints
(`render_cell_counts`: before gating / after gating / after outlier rejection), each annotation column with its role, levels and
any override, how the columns relate (crossed grid, independent, repeated controls, missing
combos), which plots/axes/pairs were made, and every non-default choice (outlier rejection, axis
limits, excluded samples, palettes, hand-rolled figures). No prose about results. When the user
fine-tunes and you re-run, update the notes to match the final driver — they must never describe
an older state.

## Outlier rejection (asked every run — the loaders never trim data)

The loaders apply **no statistical outlier rejection** — data is loaded verbatim, and trimming
happens only through the explicit `reject_outliers` transform in the driver, which returns cleaned
records feeding straight into `build_plan`/`autoplot`/any combinator, so one call cleans every
downstream plot. Step 2 always asks the user whether to use it (recommended default: Tukey fences,
per sample, on mass/volume/density, **joint across properties**):
```python
OUTLIER_SPEC = dict(method="iqr", props=["mass", "density", "vol"], scope="per_sample")
ifxm        = tk.reject_outliers(ifxm, **OUTLIER_SPEC)          # ridge / box / heatmap / timecourse
ifxm_paired = tk.reject_outliers(ifxm_paired, **OUTLIER_SPEC)   # scatter / props-vs-time
coulter     = tk.reject_outliers(coulter, method="iqr", props=["volume"], scope="per_sample")
```
(Trim `props` to the properties that exist.)

**Joint trimming (default, `joint=True`).** For each paired sample, the Tukey fence of every
selected property is computed once from the untrimmed (gated) cells. A cell outside **any** fence
is dropped from **every** property, so a cell rejected on mass is also absent from the volume and
density plots. Every plot of that sample (ridge, box, heatmap, timecourse, scatter, props-vs-time)
and its cell count show the same cells. This is a single pass: fences are not recomputed after
the joint trim. Density gets its own fence even though it's derived from mass and volume. Because
of the joint trim, a property's plotted distribution can lose cells from its interior, not just
from beyond its own fence. Unpaired mass-only / volume-only samples (and Coulter) have one property
per cell, so they are trimmed on that property alone. Both iFXM record sets get the same spec and
so the same mask; `bounds_from` is no longer needed. Use `joint=False` only if the user explicitly
asks for each property to be trimmed independently. Scatter records stay joint either way, so the
scatter rows stay aligned.

**If the user picks "customize", or brings up outlier rejection later** (e.g. "trim the density
tails", "reject outliers on volume"), do **not** guess — present a menu with `AskUserQuestion`
enumerating the full spec, then wire the answer into the driver. Ask for:

1. **Which properties** to clean — any of coulter `volume`; iFXM `mass` / `density` / `vol`
   (multi-select; can differ per property).
2. **Method** (per property allowed): `mad` (modified z-score, robust — best for density),
   `iqr` (Tukey fences, robust — the default), `percentile` (fixed-fraction tail clip).
   (Box whiskers are independent of this: they always mark the 5th/95th percentile of the data
   left after rejection — `tk.BOX_WHIS`, overridable via `draw_boxes(whis=…)`.)
3. **Scope**: `per_sample` (each sample trimmed on its own stats — default) or `pooled` (one
   global cutoff across all cells).
4. **Joint or independent**: `joint=True` (default — a paired cell outside any property's cutoff
   is dropped from all properties' plots) or `joint=False` (each property's ridges/boxes trimmed
   on its own cutoff only).
5. **Log-space?** (`log=True`) — for log-normal mass/volume so the high tail isn't over-trimmed.
6. **Method params** if they care — `mad thresh` (3.5), `iqr k` (1.5), `percentile lower/upper`
   (1/99).

Then add lines to the driver (and re-run):
```python
spec = dict(method={"density": "mad", "mass": "iqr"}, scope="per_sample")
ifxm        = tk.reject_outliers(ifxm, **spec)
ifxm_paired = tk.reject_outliers(ifxm_paired, **spec)
coulter     = tk.reject_outliers(coulter, method="iqr", props=["volume"])
```
**Processed-data export (default, iFXM — do not ask).** The template's
`tk.export_processed(ifxm, COMPILED_DIR, PROCESSED_XLSX)` runs right after the `reject_outliers`
lines and writes `processed_data/<EXP_NAME>_processed.xlsx` (in `OUT_ROOT`, beside the decks): the
cells the plots show, in `experiment_data.xlsx`'s layout. Each sample gets one sheet, under the
same sheet name. A paired sample's sheet holds `pair_cell_id`, `pair_mass_pg`, `pair_volume_fl` and
`pair_cell_density_g_per_mL`, one row per kept cell, so a cell removed by the joint trim is gone
from the whole row. An unpaired sample's sheet holds `mass_mass_pg` and/or `vol_volume_fl`. The
workbook also has a `metadata` sheet (the source rows plus `n_<block>_after_outlier_rejection`)
and a `README` sheet (source, date, outlier spec). It must get the `load_ifxm` records **after**
rejection. If rejection was declined, it holds the gated cells. It reloads with `tk.load_ifxm`.
Coulter data is not exported. Keep the block unless the user says to skip it.

Notes: always clean both iFXM record sets with the same spec. `load_ifxm_paired` records are
always trimmed jointly, even with `joint=False`, so scatter rows stay aligned. `verbose=True`
prints how many cells each sample dropped. Low-level keep-masks (`tk.outlier_mask`, `keep_mad`/`keep_iqr`/`keep_percentile`)
are exposed for bespoke logic. (k-sigma/3-std is intentionally not built in — ask if the user wants it.)

## Gotchas
- **iFXM gating is upstream**: `compile_experiment.py` gates every block (AND of all five gates)
  before writing the xlsx; the loaders only drop non-finite values, with ONE row-aligned mask per
  paired sample, so every plot of it, whether ridge, box, heatmap or scatter, shows the same cells.
  Legacy xlsx files (`bm_gate_*` / `ifxm_gate_*` columns) are gated by the loaders instead:
  `bm_gate` on `pair_mass_pg` AND `ifxm_gate` on `pair_volume_fl`, standalone blocks each on their
  own gate. No statistical outlier rejection is applied by the
  loaders (trimming is the driver's `reject_outliers` call, asked every run). Still use
  `load_ifxm_paired` for `scatter_by` / `props_vs_time` (and pass `paired_records=` to
  `autoplot`): it skips unpaired samples and carries the per-cell `time_h`.
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
