# Compiled biophysics data schema

Reference for the `plot-experiment` skill so it can map an experiment's annotation schema
without re-exploring the raw pipeline. Everything here is the **raw output** of the
`biophys_helpers` pipeline — `compile_experiment.py` (iFXM) and `annotate_coulter_samples.py`
(Coulter) — read **directly**, with no intermediate reorganization step.

> **Format note.** These scripts used to emit pandas HDF5 (`experiment_data.h5`, `data.h5`); as of
> the "get rid of h5s" change they emit an **Excel workbook** (`experiment_data.xlsx`) and **CSVs**
> instead. The loaders here read that current format. `images.h5` (raw `h5py` brightfield image
> stacks) is still HDF5 but is **not** used by these plots.

## Expected input layout

The skill plots one experiment at a time from whichever raw output dirs exist (it may have only
one half):

```
<...>_compiled/                     # from compile_experiment.py  (iFXM)
├── experiment_data.xlsx            # 'metadata' sheet + one worksheet per sample + 'README'
├── images.h5                       # raw h5py BF image stacks (NOT needed for these plots)
└── {sheet}_pair_overflow.csv       # only if a paired block exceeded Excel's row limit

<...>_coulter_sample_annotation/    # from annotate_coulter_samples.py  (Coulter)
├── metadata.csv                    # one row per Coulter sample + annotations
└── <original_input_name>.csv       # single-cell data; COLUMNS are samples, ROWS are volumes
```

`load_ifxm`/`load_ifxm_paired` take the `*_compiled/` dir (or the `.xlsx` path directly);
`load_coulter` takes the `*_coulter_sample_annotation/` dir (or its `metadata.csv`).

## Coulter — `metadata.csv` + single-cell data CSV

`metadata.csv` columns:

| column        | meaning                                                                    |
|---------------|----------------------------------------------------------------------------|
| `sample_name` | display name; also the **column header** in the data CSV (the two join by this). `rep\d+` embedded → replicate |
| *(annotations)* | every other column is **hand-added in the GUI** — e.g. `time_h`, `condition`, `drug_name` — present only if the user added it |

The data CSV (`load_coulter` auto-locates it as the one CSV that isn't `metadata.csv`; override
with `data_file=`) has **one column per sample** (header = `sample_name`) and **one row per cell**;
each column is that sample's per-cell volume (fL) array, NaN-padded to the longest column. There is
no `h5_key` and no `/data/...` key — the join is purely by sample name.

## iFXM — `experiment_data.xlsx`

The `metadata` sheet (`pd.read_excel(xlsx, sheet_name="metadata")`, one row per sample) columns
used by the toolkit:

| column                          | meaning                                                       |
|---------------------------------|---------------------------------------------------------------|
| `sample_name`                   | display name                                                  |
| `sheet_name`                    | **worksheet name** for this sample (Excel-legal, ≤31 chars, sanitized) — how the loaders find the sample's data sheet |
| `hdf5_key`                      | key into `images.h5` only (**not** the data sheet key)        |
| `condition`, `time_h`, `drug_name` (optional) | hand-added annotation columns, same meaning as Coulter |
| `gate_<attr>_lower`, `gate_<attr>_upper` | bounds of each gate already applied by `compile_experiment.py`, `<attr>` ∈ `mass_pg`, `baseline_norm`, `bl_slope`, `node_dev_mean`, `volume_fl` (NaN = not gated) — informational, the loaders don't re-apply them |
| `has_gate`                      | False = sample had no gate file and was compiled **ungated**   |
| `n_<block>_ungated`, `n_<block>_gated` | cells in each block (`mass`, `paired`, `volume`) before / after gating; feed the cell-count tables |
| `has_mass`, `has_volume`, `has_pairing`, `has_images`, `coulter_column`, `calibration_factor` | availability flags / provenance |
| legacy: `bm_gate_lower/upper`, `ifxm_gate_lower/upper`, `has_bm_gate`, `has_ifxm_gate` | older xlsx files only — blocks were written ungated, so the loaders gate with these (mass gate on `mass`, volume gate on `vol`) |

The loaders tell the two formats apart by the presence of `has_gate` / any `gate_*` column.

Each sample's worksheet (looked up by `sheet_name`) holds up to **three side-by-side blocks**
separated by one blank spacer column, distinguished by a column-name prefix:

- `vol_*` — every FXM cell (unpaired), written only for samples with no `pair_` block: currently
  the CELLGROUPED hdf5's `analysis/volume/cells` table (`vol_cell_id`, `vol_volume_fl` —
  calibrated fL — plus pass-through columns); legacy files have `vol_volume_au` / `vol_volume_fL`.
- `mass_*` — every SMR cell (unpaired): `mass_mass_pg` (+ any pass-through columns from the mass CSV).

Every block holds only gated cells (current files); see the metadata table above.
- `pair_*` — **matched cells, row-aligned per cell**, straight from a `*_CELLGROUPED.hdf5`'s
  `analysis/density/cells` table (current SMRFXMAnalysis output; the only paired format), filtered
  to `status_code == 'ok'` and mirrored essentially verbatim: `pair_cell_id`, `pair_mass_pg` (pg),
  `pair_volume_fl` (lowercase — already-calibrated fL), and **`pair_cell_density_g_per_mL`** — the
  **ABSOLUTE** density, computed by the hdf5 pipeline itself from its own `media_density_g_per_mL`
  — plus confidence/QC columns (`pair_mass_confidence_score`, `pair_volume_confidence`,
  `pair_density_confidence`, `pair_qc_fail_mask`, `pair_qc_warn_mask`, `pair_status_code`,
  `pair_status_text`, `pair_mass_cell_index`, `pair_volume_cell_index`, `pair_matched_peak_index`,
  `pair_media_density_g_per_mL`, `pair_buoyant_density`). The loader requires `pair_mass_pg`,
  `pair_volume_fl`, and `pair_cell_density_g_per_mL` to be present whenever a `pair_` block exists
  at all — it raises a clear error naming the sample if any are missing, rather than falling back.

Blocks are independent: a **paired run** has all three; a **mass-only run** has only `mass_*`; a
**volume-only run** has only `vol_*`. Read a block with `sheet.filter(regex="^pair_").dropna(how="all")`
(use an **anchored** `^pair_`/`^mass_`/`^vol_`, since `mass_` also occurs inside `pair_mass_pg`). If
a block overflowed Excel's row limit it is written full to `{sheet_name}_{prefix}_overflow.csv` and
truncated in the sheet; the loaders prefer the overflow CSV when present.

Derived properties, and where the loader sources each:

| prop key  | paired sample                                | unpaired sample (mass-only / volume-only)                          | legacy-file gate (current files: gated upstream) | units |
|-----------|-----------------------------------------------|----------------------------------------------------------------------|-------------|-------|
| `mass`    | `pair_mass_pg`                                 | `mass_mass_pg` (full SMR distribution)                                | paired: `bm_gate` AND `ifxm_gate`; unpaired: `bm_gate` | pg    |
| `density` | `pair_cell_density_g_per_mL` (already absolute) | — (empty; density **requires pairing**)                              | `bm_gate` AND `ifxm_gate` | g/mL  |
| `vol`     | `pair_volume_fl` (always present)              | `vol_volume_fl` (current) or legacy `vol_volume_fL`, else `vol_volume_au` | paired: `bm_gate` AND `ifxm_gate`; unpaired: `ifxm_gate` | fL (or AU if uncalibrated) |

**Paired-primary, standalone-fallback**: when a sample has a `pair_` block, marginal `mass`/`vol`
come from that matched subset (unchanged behavior); the standalone `mass_`/`vol_` blocks are used
**only** when there is no `pair_` block at all. `scatter_by` (via `load_ifxm_paired`) uses the
`pair_` block only, so unpaired samples never appear in scatters. **Gating is upstream**:
`compile_experiment.py` writes each block already filtered to the AND of every gate, and in both
loaders a paired sample's props share **one** row-aligned mask (all three finite), so every plot
of that sample shows the same cells. For legacy files the loaders add the gates to that mask
(`bm_gate` on `pair_mass_pg` AND `ifxm_gate` on `pair_volume_fl`; an unpaired sample's standalone
`mass_` / `vol_` blocks each on their own gate).

**Cell counts.** Every record carries `cell_counts` — per source block, cells before gating
(`n_<block>_ungated`, or the raw block length for legacy files; `None` for Coulter) and after
gating (as loaded). `reject_outliers` keeps it and appends its spec to `outlier_rejection`.
`tk.cell_count_rows(records)` / `tk.render_cell_counts(records)` turn that into the per-sample
table (before gating / after gating / after outlier rejection per prop) used in the plan, the
notes and the deck's first slide(s) (`save_pptx(..., count_tables={"iFXM": ifxm, ...})`).

**Per-cell run time (paired records only).** `pair_` rows carry no timestamp, and
`pair_mass_cell_index` does **not** index the (filtered) `mass_` block. `load_ifxm_paired` therefore
matches each paired cell to the `mass_` row with the identical `mass_pg` value and takes its
`mass_peak_time_h` (hours since run start), stored row-aligned as `record["time_h"]` (NaN if there's
no match or the mass value is ambiguous). Used by `props_vs_time`.

### Density is always absolute — no baseline to supply

`density` is `pair_cell_density_g_per_mL` directly — already absolute, computed by the hdf5
pipeline itself from its own `media_density_g_per_mL`. There is no relative-density fallback and
nothing for the driver or the skill to ask the user for.

## Annotation columns are arbitrary — roles are inferred, not fixed

Annotation columns are **hand-added in the GUI** and vary per experiment. Beyond `sample_name`
(always present) and the iFXM structural columns, *any* column may exist. The whole metadata row is
carried into each record's `meta` bag, and `infer_roles(records)` classifies every column into a
**role** that decides how it drives plots:

| role         | detected by                                              | drives                         |
|--------------|---------------------------------------------------------|--------------------------------|
| `boolean`    | values ⊆ {yes/no, true/false, 0/1, y/n, t/f}            | per-value plots + comparison   |
| `time`       | numeric + name `time`/`t_*`/suffix `_h`/`_min`/`_sec` (unit → hours; only `_sec`/`_seconds` = seconds) | sequential ordering + timecourse x |
| `ordered`    | numeric gradient (`dose`, `conc`, `passage`, …) or few distinct numeric values — **`confirm`** | ordered grouping/comparison (after confirmation) |
| `categorical`| other strings (any cardinality — no cap)                | per-value plots + comparison + facet |
| `continuous` | high-cardinality numeric, no gradient hint              | color / scatter axis only      |
| `label`      | `sample_name` / unique free text                        | labels only                    |
| `structural` | `sheet_name`, `hdf5_key`, `has_*`, `gate_*`, `n_<block>_(un)gated` (any block, e.g. `mass`/`paired`/`volume`/`vqvae`), legacy `*_gate_*`, `coulter_column`, `calibration_factor` | ignored |

Notes:
- **Booleans** (checkbox columns are literally `yes`/`no`) get the high-contrast, colorblind-safe
  `BOOL_COLORS` pair, ordered falsey→truthy.
- **Time** values are normalized to hours for ordering (`time_min` 360 → 6 h). A bare `time`/`t`
  with no unit is assumed hours and flagged for confirmation.
- **Gradient/ordered** guesses are always flagged (`confirm=True`); the driver's `ROLE_OVERRIDES`
  (or `infer_roles(overrides={col: "ordered"|"categorical"})`) pins them.
- **No cardinality cap** — every categorical value is plotted; reorganize on a later pass. Multiple
  columns are handled independently; use `cross_groups(cols=(a, b))` for a cross-product on request.
- **Grid search** — two columns varied jointly (e.g. osmolarity × drug dose) are detected by
  `suggest_grids` and listed by `render_plan` as *suggested* heatmaps; they're only plotted
  (`grid_heatmap`) once the user confirms and the pair is added to the driver's `GRID_PAIRS`.
- Replicate is parsed from `rep\d+` in `sample_name` and exposed as `meta["rep"]` for ordering/labels.
- `condition`/`drug_name` still get the reference value fixups (`drug_treat`→`drug_treated`;
  drug lowercased, `zt1a`→`zt-1a`) via `VALUE_NORMALIZERS` (overridable per loader).

### Loader arguments

Loaders no longer take per-role column args (role assignment is downstream). They take only:

| loader | args |
|--------|------|
| `load_coulter(coulter_dir, …)` | `sample_col="sample_name"`, `data_file=None` (override data-CSV auto-locate), `normalizers=VALUE_NORMALIZERS` |
| `load_ifxm` / `load_ifxm_paired(compiled_dir, …)` | `sample_col`, `sheet_col="sheet_name"` (falls back to sample), `bm_lower_col`/`bm_upper_col`/`ifxm_lower_col`/`ifxm_upper_col` (legacy files only), `normalizers` |

Grouping/ordering columns are chosen at plot time by passing `group_col=` / `series_col=` /
`facet_col=` / `time_col=` (or letting `build_plan` pick them from the inferred roles). In a legacy
file, a missing gate → unbounded (no cutoff); a missing grouping column just means no plots for it.

## Output naming grid

`<exp>` (`EXP_NAME`) is the full analysis dir name **including its date prefix**
(e.g. `2026-09-22_fl5_wnki_conc-curves`), so decks are `2026-09-22_fl5_wnki_conc-curves_figures.pptx`.
Figures go to `<exp>_fig/`; the deck to `<exp>_figures.pptx`. The default second deck (iFXM paired
samples only) puts the `props_vs_time` figures in `<exp>_vs_time_fig/` and the deck in
`<exp>_props_vs_time.pptx`. File names:

```
{datatype}_{prop}_{kind}_{col}={slug(value)}.png   # per-group detail (kind ∈ ridge|box)
{datatype}_{prop}_{kind}_by_{col}.png              # cross-group comparison (ridge & box)
{datatype}_{prop}_timecourse[_{series_col}].png    # timecourse (unit-aware, ordered)
{datatype}_{prop}_facet_{col}.png                  # compact one-figure grid
{datatype}_{propY}_vs_{propX}[_{col}={slug(value)}].png   # per-cell scatter + marginals
{datatype}_{prop}_{kind}_by_{colA}-x-{colB}.png    # cross_groups (crossing on request)
{datatype}_{prop}_heatmap_{colX}-x-{colY}.png      # grid_heatmap (grid search, user-confirmed)
{datatype}_{prop1-prop2-...}_vs_time_{slug(sample)}.png  # props_vs_time (per sample, default; in <exp>_vs_time_fig/)
```

**Processed data (iFXM, default).** `export_processed` writes
`processed_data/<exp>_processed.xlsx`, which holds the plotted cells: gated, finite, and after the
driver's `reject_outliers`. The file uses the `experiment_data.xlsx` layout above, so
`load_ifxm(<path to the .xlsx>)` reloads it. Pass the file, not the `processed_data/` folder:
given a folder, the loader looks for `experiment_data.xlsx` inside it.
- **Sample sheets:** one per sample, with the same `sheet_name`, holding only the property columns
  under their source names. A paired sample has `pair_cell_id`, `pair_mass_pg`, `pair_volume_fl`
  and `pair_cell_density_g_per_mL`, row-aligned. A joint-rejected cell is gone from the whole row.
  An unpaired sample has `mass_mass_pg` and/or `vol_volume_*`, side by side with a blank spacer.
- **`metadata`:** the source rows plus `n_<block>_after_outlier_rejection`, which is structural.
- **`README`:** source, date and outlier spec.

Loader records carry the provenance this needs: `sheet`, `columns` (prop → source column), and
for paired samples a row-aligned `cell_id` that `reject_outliers` trims with the props. Blocks over
Excel's row limit spill to `{sheet}_{prefix}_overflow.csv` beside the xlsx.

- `datatype` ∈ {`coulter`, `ifxm`}
- `prop` ∈ {`volume`} (coulter) / {`mass`, `density`, `vol`} (ifxm)
- `col` / `value` are whatever grouping column and value were used (`is_activated=yes`, `media=rpmi`, …)
- `slug(value)` lowercases and replaces non-alphanumerics with `-`
