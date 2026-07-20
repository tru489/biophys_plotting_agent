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
| `bm_gate_lower`, `bm_gate_upper`| buoyant-mass (pg) gate for `mass` (NaN when no gate)          |
| `ifxm_gate_lower`, `ifxm_gate_upper` | volume gate on **uncalibrated** volume (NaN when no gate) |
| `has_mass`, `has_volume`, `has_pairing`, `has_bm_gate`, `has_ifxm_gate`, `has_images`, `coulter_column`, `calibration_factor` | availability flags / provenance |

Each sample's worksheet (looked up by `sheet_name`) holds up to **three side-by-side blocks**
separated by one blank spacer column, distinguished by a column-name prefix:

- `vol_*` — every FXM cell (unpaired): `vol_transit_index`, `vol_volume_au`, `vol_volume_fL`.
- `mass_*` — every SMR cell (unpaired): `mass_mass_pg` (+ any other mass columns).
- `pair_*` — **matched cells, row-aligned per cell** — this is the only block the toolkit reads:
  `pair_transit_index`, `pair_mass_pg`, `pair_volume_au`, `pair_volume_fL` (calibrated volume,
  present **only** when a Coulter calibration ran), `pair_buoyant_density`.

Read a sample's paired block with `sheet.filter(regex="^pair_").dropna(how="all")` (use an anchored
`^pair_`, since `mass_` also occurs inside `pair_mass_pg`). If a paired block overflowed Excel's row
limit it is written full to `{sheet_name}_pair_overflow.csv` and truncated in the sheet; the loaders
prefer the overflow CSV when present. Samples with **no paired block** are skipped (e.g. a no-iFXM
proliferating control).

Derived properties (as built by the loaders, all from the one `pair_` block):

| prop key    | source (pair_ block, prefix stripped)     | gate            | units  |
|-------------|-------------------------------------------|-----------------|--------|
| `mass`      | `mass_pg`                                 | `bm_gate`       | pg     |
| `density`   | `buoyant_density + baseline_density`      | `ifxm_gate`     | g/mL   |
| `vol_cal`   | `volume_fL` (empty if uncalibrated)       | `ifxm_gate`     | fL     |
| `vol_uncal` | `volume_au`                               | `ifxm_gate`     | AU     |

The three iFXM `ifxm_gate` properties share **one** mask computed on the uncalibrated `volume_au`
(cells are row-aligned within the paired block). Do not re-gate per property.

### `baseline_density` — not in any file

`density = buoyant_density + baseline_density`. `buoyant_density` (`pair_buoyant_density`) is
**RELATIVE**; the baseline (g/mL) is **experiment-specific and stored nowhere in the data**.
`load_ifxm` requires it explicitly (no default). The skill must ask the user for the correct value
per experiment. The FL5 reference experiments used `1.008`.

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
| `structural` | `sheet_name`, `hdf5_key`, `has_*`, `*_gate_*`, `coulter_column`, `calibration_factor` | ignored |

Notes:
- **Booleans** (checkbox columns are literally `yes`/`no`) get the high-contrast, colorblind-safe
  `BOOL_COLORS` pair, ordered falsey→truthy.
- **Time** values are normalized to hours for ordering (`time_min` 360 → 6 h). A bare `time`/`t`
  with no unit is assumed hours and flagged for confirmation.
- **Gradient/ordered** guesses are always flagged (`confirm=True`); the driver's `ROLE_OVERRIDES`
  (or `infer_roles(overrides={col: "ordered"|"categorical"})`) pins them.
- **No cardinality cap** — every categorical value is plotted; reorganize on a later pass. Multiple
  columns are handled independently; use `cross_groups(cols=(a, b))` for a cross-product on request.
- Replicate is parsed from `rep\d+` in `sample_name` and exposed as `meta["rep"]` for ordering/labels.
- `condition`/`drug_name` still get the reference value fixups (`drug_treat`→`drug_treated`;
  drug lowercased, `zt1a`→`zt-1a`) via `VALUE_NORMALIZERS` (overridable per loader).

### Loader arguments

Loaders no longer take per-role column args (role assignment is downstream). They take only:

| loader | args |
|--------|------|
| `load_coulter(coulter_dir, …)` | `sample_col="sample_name"`, `data_file=None` (override data-CSV auto-locate), `normalizers=VALUE_NORMALIZERS` |
| `load_ifxm` / `load_ifxm_paired(compiled_dir, baseline_density, …)` | `sample_col`, `sheet_col="sheet_name"` (falls back to sample), `bm_lower_col`/`bm_upper_col`/`ifxm_lower_col`/`ifxm_upper_col`, `normalizers` |

Grouping/ordering columns are chosen at plot time by passing `group_col=` / `series_col=` /
`facet_col=` / `time_col=` (or letting `build_plan` pick them from the inferred roles). Missing gate
→ unbounded (no cutoff); a missing grouping column just means no plots for it.

## Output naming grid

Figures go to `<exp>_fig/`; the deck to `<exp>_figures.pptx`. File names:

```
{datatype}_{prop}_{kind}_{col}={slug(value)}.png   # per-group detail (kind ∈ ridge|box)
{datatype}_{prop}_{kind}_by_{col}.png              # cross-group comparison (ridge & box)
{datatype}_{prop}_timecourse[_{series_col}].png    # timecourse (unit-aware, ordered)
{datatype}_{prop}_facet_{col}.png                  # compact one-figure grid
{datatype}_{propY}_vs_{propX}[_{col}={slug(value)}].png   # per-cell scatter + marginals
{datatype}_{prop}_{kind}_by_{colA}-x-{colB}.png    # cross_groups (crossing on request)
```

- `datatype` ∈ {`coulter`, `ifxm`}
- `prop` ∈ {`volume`} (coulter) / {`mass`, `density`, `vol_cal`, `vol_uncal`} (ifxm)
- `col` / `value` are whatever grouping column and value were used (`is_activated=yes`, `media=rpmi`, …)
- `slug(value)` lowercases and replaces non-alphanumerics with `-`
