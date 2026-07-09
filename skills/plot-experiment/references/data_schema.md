# Compiled biophysics data schema

Reference for the `plot-experiment` skill so it can map an experiment's annotation schema
without re-exploring the raw pipeline. Everything here is the **compiled + annotated** output of
the `biophys_helpers` pipeline (`compile_experiment.py`, `annotate_coulter_samples.py`,
`pair_bm_runs.py`), reorganized into a per-experiment analysis directory.

## Expected input layout

The skill plots one experiment at a time from a `<exp>_data/` directory:

```
<exp>_data/
├── coulter/
│   ├── metadata.csv          # one row per Coulter sample + annotations
│   └── data.h5               # pandas HDFStore, /data/{h5_key} per sample
└── ifxm/
    ├── experiment_data.h5    # pandas HDFStore, /metadata + /samples/{hdf5_key}/...
    └── images.h5             # raw h5py BF image stacks (NOT needed for these plots)
```

An experiment may have only `coulter/` or only `ifxm/`. If instead you're pointed at a raw
`*_compiled/` dir (from `compile_experiment.py`), that dir already holds `experiment_data.h5` +
`images.h5` — those become the `ifxm/` half; the Coulter half comes from an
`annotate_coulter_samples.py` output (`metadata.csv` + `data.h5`).

## Coulter — `coulter/metadata.csv` + `coulter/data.h5`

`metadata.csv` columns:

| column       | meaning                                                        |
|--------------|----------------------------------------------------------------|
| `sample_name`| display name (`rep\d+` embedded → replicate)                   |
| `h5_key`     | HDF5-safe key into `data.h5` at `/data/{h5_key}`               |
| `time_h`     | timepoint in **decimal hours** (`.25`, `1`, `6.33`, `24`; `0` for starved) |
| `condition`  | free-text condition (`activated`, `starved`, ...)             |
| `drug_name`  | **optional** — drug for treated samples (may be absent entirely) |

`data.h5`: `pd.read_hdf(h5, "/data/{h5_key}")`; the **first column** is the per-cell volume (fL)
array. Coulter volume is already gated upstream.

## iFXM — `ifxm/experiment_data.h5`

`/metadata` (one row per sample) columns used by the toolkit:

| column                          | meaning                                       |
|---------------------------------|-----------------------------------------------|
| `sample_name`                   | display name                                  |
| `hdf5_key`                      | key prefix into `/samples/{hdf5_key}/...`     |
| `condition`, `time_h`, `drug_name` (optional) | same annotation meaning as Coulter |
| `bm_gate_lower`, `bm_gate_upper`| buoyant-mass (pg) gate for `mass`             |
| `ifxm_gate_lower`, `ifxm_gate_upper` | volume gate applied on **uncalibrated** volume |
| `has_mass`, `has_volume`, `has_pairing`, ... | availability flags               |

Per-sample subtables (only `pairing` + `volume_calibrated` are needed here):

- `/samples/{hdf5_key}/pairing` — per-cell paired rows. Columns used:
  `matched_mass` (pg), `buoyant_density`, `volume` (uncalibrated, vol_au).
- `/samples/{hdf5_key}/volume_calibrated` — includes `volume_fL` (calibrated, fL).
  Present only when Coulter calibration ran; **samples lacking `pairing` or `volume_calibrated`
  are skipped** (e.g. a no-iFXM proliferating control).

Derived properties (as built by the loaders):

| prop key    | source                                   | gate            | units  |
|-------------|------------------------------------------|-----------------|--------|
| `mass`      | `matched_mass`                           | `bm_gate`       | pg     |
| `density`   | `buoyant_density + baseline_density`     | `ifxm_gate`     | g/mL   |
| `vol_cal`   | `volume_fL`                              | `ifxm_gate`     | fL     |
| `vol_uncal` | `volume`                                 | `ifxm_gate`     | fL/AU  |

The three iFXM `ifxm_gate` properties share **one** mask computed on the uncalibrated `volume`
(cells are row-aligned across the paired subtables). Do not re-gate per property.

### `baseline_density` — not in any file

`density = buoyant_density + baseline_density`. The baseline (g/mL) is **experiment-specific and
stored nowhere in the data**. `load_ifxm` requires it explicitly (no default). The skill must ask
the user for the correct value per experiment. The FL5 reference experiments used `1.008`.

## Annotation conventions (the part that varies per experiment)

Conditions, drugs, and timepoints are **hand-added columns**, not a fixed schema:

- `condition` is free-text. Normalized by `_norm_cond`: `drug_treat`/`drug_treated` → `drug_treated`.
- `drug_name` is optional and only present when there is a drug arm. Normalized by `_norm_drug`:
  lowercased, `zt1a` → `zt-1a`.
- `time_h` is decimal hours. Starved samples typically share `time_h == 0`, so ridge/box code
  labels starved rows by full `sample_name` to tell replicates apart.
- Replicate is parsed from `rep\d+` in `sample_name` (default `rep1`).

The sample names encode the same info (e.g. `activated_1uM-wnk463_20h15min_rep1`), but the skill
should read the **metadata columns**, not parse names, to discover what conditions / drugs /
timepoints are actually present in a given experiment.

### Column names are defaults, not guarantees

Every column name above is the *reference* name and the loader default — **not fixed**. The
loaders (`load_coulter`, `load_ifxm`, `load_ifxm_paired`) accept each role as a keyword argument
so a differently-named schema just needs the mapping passed in:

| role         | coulter default | iFXM default      | loader arg        |
|--------------|-----------------|-------------------|-------------------|
| sample name  | `sample_name`   | `sample_name`     | `sample_col` (required) |
| h5 key       | `h5_key`        | `hdf5_key`        | `key_col` (falls back to sample) |
| condition    | `condition`     | `condition`       | `condition_col`   |
| time (hours) | `time_h`        | `time_h`          | `time_col`        |
| drug         | `drug_name`     | `drug_name`       | `drug_col`        |
| bm gate      | —               | `bm_gate_lower/upper`   | `bm_lower_col` / `bm_upper_col` |
| ifxm gate    | —               | `ifxm_gate_lower/upper` | `ifxm_lower_col` / `ifxm_upper_col` |

Graceful degradation when a column is absent: **condition** → one unnamed group; **time** → 0.0;
**drug** → no drug figures; **gate** (missing or NaN) → unbounded (no cutoff), so an ungated
experiment is not silently emptied. Unknown **condition/drug values** are auto-assigned distinct
colors by `build_color_map`, seeded from `COND_COLORS`/`DRUG_COLORS`. Value normalization
(`_norm_cond`, `_norm_drug`) is overridable via the loaders' `normalize_cond`/`normalize_drug`
arguments (pass `lambda x: x` to disable).

## Output naming grid

Figures go to `<exp>_fig/`; the deck to `<exp>_figures.pptx`. File names:

```
{datatype}_{metric}_{plottype}[_{cond|drug}].png
{datatype}_{propY}_vs_{propX}_scatter[_{cond}].png
```

- `datatype` ∈ {`coulter`, `ifxm`}
- `metric` ∈ {`volume`} (coulter) / {`mass`, `density`, `vol_cal`, `vol_uncal`} (ifxm)
- `plottype` ∈ {`ridge`, `box`, `timecourse`}; drug arm adds `drug_ridge`/`drug_box`/`drug_timecourse`
- condition groups are whatever `condition` values exist (`activated`, `starved`,
  `drug_treated`, `proliferating`, ...)
```
