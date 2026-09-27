"""
biophys_plot_toolkit — reusable plotting library for compiled FXM/SMR/Coulter experiments.

Ported from a hand-written reference analysis and generalized so a short per-experiment driver can
load the raw biophys_helpers outputs and produce the standard figure grid without re-deriving the
plotting internals.

Three layers, low -> high:
  * Loaders           -> records: {sample, props:{name: ndarray}, meta:{col: value, ..., rep}}
                         `meta` carries EVERY annotation column, so any column can drive plots.
  * infer_roles       -> classifies each metadata column into a role (boolean / categorical /
                         time / ordered / continuous / label / structural) so the driver (and
                         Claude) can decide how to group, compare, order and color the data.
  * draw_*            -> low-level primitives: draw on a passed-in axis, no semantics.
  * plot_grouped / compare_groups / timecourse_by / scatter_by / facet / cross_groups /
    grid_heatmap      -> mid-level combinators parameterized by WHICH column(s) to use.
                         (suggest_grids flags grid-search column pairs for grid_heatmap.)
  * build_plan / render_plan / autoplot
                      -> high-level: infer a plot plan from the roles, show it, execute it.

The loaders read the RAW biophys_helpers outputs directly (no reorg step):
  Coulter — annotate_coulter_samples.py's '*_coulter_sample_annotation/' dir (metadata.csv +
            a single-cell data CSV; columns are samples, rows are per-cell volumes).
            One property "volume" (fL, gated upstream).
  iFXM    — compile_experiment.py's '*_compiled/experiment_data.xlsx' (a 'metadata' sheet + one
            worksheet per sample). A sample's PAIRED ('pair_') block gives, row-aligned, per matched
            cell — straight from a *_CELLGROUPED.hdf5's analysis/density/cells table (current
            SMRFXMAnalysis output; the only paired format): mass (pair_mass_pg, pg), volume
            (pair_volume_fl, already-calibrated fL), and ABSOLUTE density
            (pair_cell_density_g_per_mL, g/mL — computed by the hdf5 pipeline itself from its own
            media_density_g_per_mL; there is no baseline_density to supply). A sample with no
            paired block (mass-only / volume-only run) falls back to the standalone MASS ('mass_')
            and/or VOLUME ('vol_') blocks for the distribution props (density needs pairing, so it
            is empty there). scatter (load_ifxm_paired) uses paired samples only.

No statistical outlier rejection is applied anywhere — the loaders drop only non-finite (NaN/inf)
values. The only intentional data exclusions are the metadata-driven gates (bm_gate / ifxm_gate)
and samples skipped because they lack a paired ('pair_') block.
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from pathlib import Path

# ---------------------------------------------------------------------------
# Styling defaults (override in the driver if an experiment needs different keys)
# ---------------------------------------------------------------------------
COND_COLORS = {
    "activated":     "#0072B2",
    "starved":       "#E69F00",
    "drug_treated":  "#009E73",
    "proliferating": "#CC79A7",
}
DRUG_COLORS = {  # keys match the lowercased output of _norm_drug
    "dmso":       "#0072B2",
    "1um-wnk463": "#E69F00",
    "2um-zt-1a":  "#009E73",
}
# High-contrast, colorblind-accommodating pair for booleans: [falsey, truthy] (Okabe-Ito).
BOOL_COLORS = ["#0072B2", "#D55E00"]  # blue (false) vs vermillion (true)
FALLBACK_COLOR = "#999999"

# Colorblind-safe cycle (Okabe-Ito + a few extensions) for auto-assigning colors to unknown values.
_AUTO_PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9",
    "#F0E442", "#000000", "#8C613C", "#666666", "#B8C9D0", "#5975A4",
]

# (prop_key, axis_label) lists — defaults matching the reference experiments.
COULTER_PROPS = [("volume", "Volume (fL)")]
IFXM_PROPS = [
    ("mass",    "Buoyant mass (pg)"),
    ("density", "Density (g/mL)"),
    ("vol",     "Volume (fL)"),
]

# ---------------------------------------------------------------------------
# Value normalization
# ---------------------------------------------------------------------------
def _is_blank(x) -> bool:
    return x is None or (isinstance(x, float) and np.isnan(x)) or str(x).strip() == ""


def _norm_cond(c) -> str:
    c = "" if _is_blank(c) else str(c).strip()
    return "drug_treated" if c in ("drug_treat", "drug_treated") else c


def _norm_drug(d) -> str:
    if _is_blank(d):
        return ""
    d = str(d).strip().lower()
    return d.replace("zt1a", "zt-1a")


def _rep(sample_name: str) -> str:
    m = re.search(r"rep(\d+)", str(sample_name))
    return f"rep{m.group(1)}" if m else "rep1"


# Applied to matching metadata columns at load time (by column name). Overridable per loader.
# Keeps the reference condition/drug fixups while everything else passes through stripped.
VALUE_NORMALIZERS = {"condition": _norm_cond, "drug_name": _norm_drug}


def _clean_value(v):
    """Light normalization for a raw metadata cell: blanks -> '', strings stripped, else as-is."""
    if _is_blank(v):
        return ""
    return v.strip() if isinstance(v, str) else v


def _build_meta(row, sample, normalizers) -> dict:
    """The generic annotation bag for one record: every column, lightly normalized, plus rep."""
    meta = {}
    for col in row.index:
        val = _clean_value(row[col])
        fn = normalizers.get(col)
        meta[col] = fn(val) if fn else val
    meta["rep"] = _rep(sample)
    return meta


# ---------------------------------------------------------------------------
# Robust metadata access
# ---------------------------------------------------------------------------
def rget(r, col, default=None):
    """Value of annotation `col` for record `r` (from its meta bag), or `default`."""
    v = r.get("meta", {}).get(col, default)
    return default if (col != "rep" and _is_blank(v)) else v


def _get(row, col, default=""):
    if col and col in row.index and not _is_blank(row[col]):
        return row[col]
    return default


def _get_float(row, col, default):
    v = _get(row, col, None)
    try:
        return float(v) if v is not None else default
    except (TypeError, ValueError):
        return default


def _gate_bound(row, col, default):
    """A gate bound as float, or `default` (±inf) when the column is missing or NaN, so a
    missing/ungated experiment means 'no cutoff' rather than dropping every cell."""
    if col and col in row.index and not _is_blank(row[col]):
        try:
            return float(row[col])
        except (TypeError, ValueError):
            return default
    return default


def _require_col(meta, col, role):
    if col not in meta.columns:
        raise KeyError(
            f"could not find the {role} column '{col}'. Available columns: "
            f"{list(meta.columns)}. Pass the right name via the loader's *_col argument."
        )


def _to_floats(values):
    """Return values as a list of floats if every non-blank value is numeric, else None."""
    out = []
    for v in values:
        if _is_blank(v):
            continue
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            return None
    return out


# ---------------------------------------------------------------------------
# Color maps
# ---------------------------------------------------------------------------
def build_color_map(values, base=None):
    """Color map covering every value in `values`. Known keys in `base` keep their color; unknown
    values get the next distinct palette color (deterministic, order-stable)."""
    out = dict(base or {})
    used = set(out.values())
    i = 0
    for v in values:
        if v in out:
            continue
        while i < len(_AUTO_PALETTE) and _AUTO_PALETTE[i] in used:
            i += 1
        color = _AUTO_PALETTE[i % len(_AUTO_PALETTE)]
        out[v] = color
        used.add(color)
        i += 1
    return out


def color_map_for(records, col, roles=None, base=None):
    """Color map for the values of `col` present in `records`. Booleans get the high-contrast
    BOOL_COLORS pair; `condition`/`drug`-named columns seed from COND_COLORS / DRUG_COLORS."""
    role = _role_for(records, col, roles)
    values = group_order(records, col, roles)
    if role["role"] == "boolean" and len(values) <= 2:
        # falsey -> BOOL_COLORS[0], truthy -> BOOL_COLORS[1] (values already ordered falsey->truthy)
        return {v: BOOL_COLORS[i] for i, v in enumerate(values)}
    if base is None:
        if col == "condition":
            base = COND_COLORS
        elif "drug" in col.lower():
            base = DRUG_COLORS
    return build_color_map(values, base)


# ---------------------------------------------------------------------------
# Loaders  ->  list of records {sample, props:{name: arr}, meta:{col: value, ..., rep}}
# ---------------------------------------------------------------------------

# Columns inside each iFXM sample sheet's blocks, after the block prefix is stripped. This is the
# compile_experiment.py output contract. Each sample sheet holds up to three side-by-side blocks:
#   PAIRED ('pair_')  — matched cells, row-aligned per cell, straight from a *_CELLGROUPED.hdf5's
#                       analysis/density/cells table (current SMRFXMAnalysis output; the only
#                       paired format): mass_pg (pg), volume_fl (already-calibrated fL),
#                       cell_density_g_per_mL (ABSOLUTE density, computed by the hdf5 pipeline
#                       itself from its own media_density_g_per_mL).
#   MASS   ('mass_')  — every SMR cell (unpaired): mass_pg (+ pass-through mass_* columns).
#   VOLUME ('vol_')   — every FXM cell (unpaired): volume_au, volume_fL.
# Density is pairing-only; the standalone blocks never carry it.
_PAIR_MASS = "mass_pg"                # buoyant mass (pg)
_PAIR_VOL  = "volume_fl"              # already-calibrated volume (fL) -> single `vol` prop
_PAIR_DENS = "cell_density_g_per_mL"  # ABSOLUTE density (g/mL)
# Standalone-block value columns (after the 'mass_'/'vol_' prefix is stripped) — unpaired
# mass-only / volume-only runs only; every paired sample uses the PAIRED columns above instead.
# The VOLUME block can carry a calibrated (fL) and/or uncalibrated (AU) reading; both feed the
# single `vol` prop, preferring calibrated fL when both are present.
_MASS_STANDALONE = "mass_pg"    # MASS block buoyant mass (pg)   -> 'mass_mass_pg' on the sheet
_VOL_UNCAL = "volume_au"        # VOLUME block uncalibrated (AU) -> 'vol_volume_au'
_VOL_CAL   = "volume_fL"        # VOLUME block calibrated (fL)   -> 'vol_volume_fL' (if calibrated)


def _find_coulter_data_csv(d: Path, meta_path: Path) -> Path:
    """The single-cell data CSV keeps the input file's name, so it is the one CSV in the
    directory that is not metadata.csv. Raise clearly if it is ambiguous or missing."""
    csvs = [p for p in sorted(d.glob("*.csv")) if p.resolve() != meta_path.resolve()]
    if len(csvs) == 1:
        return csvs[0]
    if not csvs:
        raise FileNotFoundError(
            f"no single-cell data CSV found in {d} (only metadata.csv). Pass data_file=.")
    raise ValueError(
        f"multiple candidate data CSVs in {d}: {[p.name for p in csvs]}; pass data_file= to pick one.")


def load_coulter(coulter_dir, *, sample_col="sample_name", data_file=None,
                 normalizers=VALUE_NORMALIZERS) -> list:
    """Load Coulter single-cell volumes from a '*_coulter_sample_annotation/' directory produced
    by annotate_coulter_samples.py (or point straight at its metadata.csv).

    metadata.csv has `sample_name` plus whatever annotation columns were hand-added. The per-cell
    data is a separate CSV whose COLUMNS are samples (headers == sample_name values) and whose ROWS
    are single-cell volumes; it is auto-located as the non-metadata CSV in the dir (`data_file`
    overrides). Volume units follow the input CSV (fL in the standard Coulter pipeline).

    Every metadata column is carried into each record's `meta` bag (role assignment happens later
    via infer_roles), so only `sample_col` is required.
    """
    p = Path(coulter_dir)
    if p.is_file():
        meta_path, d = p, p.parent
    else:
        meta_path, d = p / "metadata.csv", p
    if not meta_path.exists():
        raise FileNotFoundError(
            f"metadata.csv not found at {meta_path}. Point load_coulter at the "
            f"'*_coulter_sample_annotation/' directory from annotate_coulter_samples.py.")
    meta = pd.read_csv(meta_path)
    _require_col(meta, sample_col, "sample-name")

    data_path = Path(data_file) if data_file else _find_coulter_data_csv(d, meta_path)
    data = pd.read_csv(data_path)

    recs = []
    for _, r in meta.iterrows():
        sample = _get(r, sample_col, "")
        if sample not in data.columns:
            raise KeyError(
                f"sample {sample!r} has no column in {data_path.name}. "
                f"Data columns: {list(data.columns)[:12]}{' ...' if data.shape[1] > 12 else ''}")
        arr = data[sample].to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]  # drop NaN padding / non-finite; no statistical outlier rejection
        recs.append({
            "sample": sample,
            "props":  {"volume": arr},
            "meta":   _build_meta(r, sample, normalizers),
        })
    return recs


def _open_ifxm_xlsx(compiled_dir) -> Path:
    """Resolve experiment_data.xlsx from a '*_compiled/' dir (or accept the .xlsx path directly)."""
    p = Path(compiled_dir)
    xlsx = p if p.suffix.lower() == ".xlsx" else p / "experiment_data.xlsx"
    if not xlsx.exists():
        raise FileNotFoundError(
            f"experiment_data.xlsx not found at {xlsx}. Point load_ifxm at the "
            f"'*_compiled/' directory produced by compile_experiment.py.")
    return xlsx


def _read_block(xls, xlsx_path: Path, sheet_name, prefix, _sheet_cache=None) -> pd.DataFrame:
    """One prefixed block ('pair'|'mass'|'vol') of a sample sheet as a DataFrame (rows dropped where
    all block columns are NaN, prefix stripped), or None if that block is absent. A block that
    overflowed Excel's row limit is written by compile_experiment.py to a sibling
    '{sheet}_{prefix}_overflow.csv'; that full copy is preferred when present. `_sheet_cache` is an
    optional {sheet_name: DataFrame} dict so the three blocks of a sheet share one read_excel."""
    overflow = Path(xlsx_path).parent / f"{sheet_name}_{prefix}_overflow.csv"
    if overflow.exists():
        blk = pd.read_csv(overflow)
    else:
        if _sheet_cache is not None and sheet_name in _sheet_cache:
            sheet = _sheet_cache[sheet_name]
        else:
            sheet = pd.read_excel(xls, sheet_name=sheet_name)
            if _sheet_cache is not None:
                _sheet_cache[sheet_name] = sheet
        blk = sheet.filter(regex=rf"^{prefix}_").dropna(how="all")
    if blk.empty:
        return None
    pre = f"{prefix}_"
    blk = blk.rename(columns=lambda c: c[len(pre):] if str(c).startswith(pre) else c)
    return blk


_EMPTY = np.array([], dtype=float)


def _gate_clean(a, mask):
    """Apply a keep-mask (if lengths match) then drop non-finite. Empty arrays pass through."""
    a = a[mask] if a.size == mask.size else a
    return a[np.isfinite(a)]


def _load_ifxm_records(compiled_dir, sample_col, sheet_col, gate_cols,
                       normalizers, paired):
    """Shared iFXM reader. `paired`=True keeps all props row-aligned under one mask from the PAIRED
    block (for scatter); unpaired samples are skipped. `paired`=False builds distribution props:
    from the PAIRED block when present (matched subset — unchanged), else falling back to the
    standalone MASS / VOLUME blocks so mass-only and volume-only runs still load (density stays
    pairing-only)."""
    bm_lo_c, bm_hi_c, ix_lo_c, ix_hi_c = gate_cols
    xlsx = _open_ifxm_xlsx(compiled_dir)
    recs = []
    with pd.ExcelFile(xlsx) as xls:
        meta = pd.read_excel(xls, sheet_name="metadata")
        _require_col(meta, sample_col, "sample-name")
        skey = sheet_col if sheet_col in meta.columns else sample_col
        cache = {}

        for _, r in meta.iterrows():
            pair = _read_block(xls, xlsx, r[skey], "pair", cache)
            has_pair = pair is not None
            if has_pair:
                missing = [c for c in (_PAIR_MASS, _PAIR_VOL, _PAIR_DENS) if c not in pair.columns]
                if missing:
                    raise KeyError(
                        f"sample {r[skey]!r} has a pair_ block but is missing "
                        f"{[f'pair_{c}' for c in missing]}. Every paired sample must carry "
                        f"pair_mass_pg, pair_volume_fl, and pair_cell_density_g_per_mL.")

            if paired and not has_pair:
                continue  # scatter needs a paired block; unpaired samples have no row-aligned pairs

            bm_lo = _gate_bound(r, bm_lo_c, -np.inf)
            bm_hi = _gate_bound(r, bm_hi_c,  np.inf)
            ix_lo = _gate_bound(r, ix_lo_c, -np.inf)
            ix_hi = _gate_bound(r, ix_hi_c,  np.inf)
            sample = _get(r, sample_col, "")
            meta_bag = _build_meta(r, sample, normalizers)

            if has_pair:
                mass = pair[_PAIR_MASS].to_numpy(dtype=float)
                vol  = pair[_PAIR_VOL].to_numpy(dtype=float)
                dens = pair[_PAIR_DENS].to_numpy(dtype=float)

                if paired:
                    # single mask over all three props keeps them length-matched & row-aligned
                    mask = (np.isfinite(mass) & np.isfinite(dens) & np.isfinite(vol)
                            & (mass >= bm_lo) & (mass <= bm_hi)
                            & (vol >= ix_lo) & (vol <= ix_hi))
                    props = {
                        "mass":    mass[mask],
                        "density": dens[mask],
                        "vol":     vol[mask],
                    }
                else:
                    bm_mask = np.isfinite(mass) & (mass >= bm_lo) & (mass <= bm_hi)
                    ix_mask = np.isfinite(vol) & (vol >= ix_lo) & (vol <= ix_hi)
                    props = {
                        "mass":    _gate_clean(mass, bm_mask),
                        "density": _gate_clean(dens, ix_mask),
                        "vol":     _gate_clean(vol, ix_mask),
                    }
            else:
                # Unpaired sample (mass-only / volume-only): fall back to the standalone blocks.
                massblk = _read_block(xls, xlsx, r[skey], "mass", cache)
                volblk  = _read_block(xls, xlsx, r[skey], "vol", cache)
                if massblk is None and volblk is None:
                    continue  # sample has no tabular iFXM data at all

                mass = massblk[_MASS_STANDALONE].to_numpy(dtype=float) \
                    if massblk is not None and _MASS_STANDALONE in massblk.columns else _EMPTY
                # Single `vol` prop: prefer the calibrated (fL) reading, else the uncalibrated (AU) one.
                if volblk is not None and _VOL_CAL in volblk.columns:
                    vol = volblk[_VOL_CAL].to_numpy(dtype=float)
                elif volblk is not None and _VOL_UNCAL in volblk.columns:
                    vol = volblk[_VOL_UNCAL].to_numpy(dtype=float)
                else:
                    vol = _EMPTY

                bm_mask = np.isfinite(mass) & (mass >= bm_lo) & (mass <= bm_hi)
                ix_mask = np.isfinite(vol) & (vol >= ix_lo) & (vol <= ix_hi)
                props = {
                    "mass":    _gate_clean(mass, bm_mask),
                    "density": _EMPTY,               # density requires pairing
                    "vol":     _gate_clean(vol, ix_mask),
                }
            recs.append({"sample": sample, "props": props, "meta": meta_bag})
    return recs


def load_ifxm(compiled_dir, *, sample_col="sample_name",
              sheet_col="sheet_name", bm_lower_col="bm_gate_lower", bm_upper_col="bm_gate_upper",
              ifxm_lower_col="ifxm_gate_lower", ifxm_upper_col="ifxm_gate_upper",
              normalizers=VALUE_NORMALIZERS) -> list:
    """Load iFXM distribution data from a '*_compiled/' dir's experiment_data.xlsx.

    For each sample (worksheet named by `sheet_col`): if it has a PAIRED ('pair_') block, mass /
    density / vol come from that matched subset (mass bm-gated; density/vol ifxm-gated on volume) —
    pair_mass_pg (pg), pair_volume_fl (already-calibrated fL), and pair_cell_density_g_per_mL
    (ABSOLUTE density, computed by the hdf5 pipeline itself). A pair_ block missing any of those
    three columns raises — every paired sample must carry all three. If a sample has NO paired block
    (a mass-only or volume-only run), the standalone MASS ('mass_') and/or VOLUME ('vol_') blocks
    are used instead — `mass` from the full SMR distribution and/or `vol` from the full FXM
    distribution (calibrated fL preferred, else uncalibrated AU), with `density` empty (density
    requires pairing). Samples with no tabular iFXM data are skipped. Every metadata column is
    carried into each record's `meta`.
    """
    return _load_ifxm_records(
        compiled_dir, sample_col, sheet_col,
        (bm_lower_col, bm_upper_col, ifxm_lower_col, ifxm_upper_col), normalizers, paired=False)


def load_ifxm_paired(compiled_dir, *, sample_col="sample_name",
                     sheet_col="sheet_name", bm_lower_col="bm_gate_lower",
                     bm_upper_col="bm_gate_upper", ifxm_lower_col="ifxm_gate_lower",
                     ifxm_upper_col="ifxm_gate_upper", normalizers=VALUE_NORMALIZERS) -> list:
    """Like load_ifxm, but keeps per-cell arrays row-ALIGNED across properties (one shared mask from
    the PAIRED block: pair_mass_pg, pair_volume_fl, pair_cell_density_g_per_mL), so a cell's mass /
    density / volume stay paired. Use for scatter_by. Only samples with a paired block appear
    (unpaired mass-only / volume-only samples are skipped — there is nothing to correlate)."""
    return _load_ifxm_records(
        compiled_dir, sample_col, sheet_col,
        (bm_lower_col, bm_upper_col, ifxm_lower_col, ifxm_upper_col), normalizers, paired=True)


# ---------------------------------------------------------------------------
# Outlier rejection — OPT-IN. The loaders NEVER call these; the data is loaded verbatim (only
# non-finite values + metadata gates removed). Call reject_outliers yourself, per experiment or
# per sample, to trim distributions for visualization. Each method derives keep-bounds (lo, hi)
# from a reference array, so 'per_sample' vs 'pooled' is the same code with a different reference.
# ---------------------------------------------------------------------------
_OUTLIER_DEFAULTS = {
    "mad":        {"thresh": 3.5},          # modified z-score: |0.6745*(x-median)/MAD| <= thresh
    "iqr":        {"k": 1.5},               # Tukey fences: [Q1 - k*IQR, Q3 + k*IQR]
    "percentile": {"lower": 1.0, "upper": 99.0},   # drop below `lower` / above `upper` percentile
}


def _outlier_bounds(ref: np.ndarray, method: str, params: dict):
    """Keep-bounds (lo, hi) in the working space, from finite reference values `ref`."""
    p = dict(_OUTLIER_DEFAULTS[method], **(params or {}))
    if ref.size == 0:
        return -np.inf, np.inf
    if method == "mad":
        med = np.median(ref)
        mad = np.median(np.abs(ref - med))
        if mad == 0:
            return -np.inf, np.inf
        half = p["thresh"] * mad / 0.6745
        return med - half, med + half
    if method == "iqr":
        q1, q3 = np.percentile(ref, [25, 75])
        iqr = q3 - q1
        return q1 - p["k"] * iqr, q3 + p["k"] * iqr
    if method == "percentile":
        lo = np.percentile(ref, p["lower"]) if p["lower"] > 0 else -np.inf
        hi = np.percentile(ref, p["upper"]) if p["upper"] < 100 else np.inf
        return lo, hi
    raise ValueError(f"unknown outlier method {method!r} (use 'mad', 'iqr', or 'percentile')")


def _to_outlier_space(a, log):
    """Values in the working space + a validity mask. log-space keeps strictly-positive finites."""
    a = np.asarray(a, float)
    finite = np.isfinite(a)
    if not log:
        return a, finite
    valid = finite & (a > 0)
    out = np.full(a.shape, np.nan)
    out[valid] = np.log(a[valid])
    return out, valid


def outlier_mask(a, method="iqr", *, log=False, ref=None, **params) -> np.ndarray:
    """Boolean KEEP-mask for 1-D `a` (True = keep). Bounds are computed from `ref` (default: `a`
    itself → per-sample). method: 'mad' | 'iqr' | 'percentile'. params (override the defaults):
    mad `thresh=3.5`; iqr `k=1.5`; percentile `lower=1, upper=99`. `log=True` computes bounds in
    log-space (non-positive values are dropped). Non-finite values are always rejected."""
    vals, valid = _to_outlier_space(a, log)
    rvals, rvalid = _to_outlier_space(a if ref is None else ref, log)
    lo, hi = _outlier_bounds(rvals[rvalid], method, params)
    return valid & (vals >= lo) & (vals <= hi)


def keep_mad(a, thresh=3.5, *, log=False, ref=None):
    return outlier_mask(a, "mad", log=log, ref=ref, thresh=thresh)


def keep_iqr(a, k=1.5, *, log=False, ref=None):
    return outlier_mask(a, "iqr", log=log, ref=ref, k=k)


def keep_percentile(a, lower=1.0, upper=99.0, *, log=False, ref=None):
    return outlier_mask(a, "percentile", log=log, ref=ref, lower=lower, upper=upper)


def reject_outliers(records, method="iqr", *, props=None, paired=False, scope="per_sample",
                    log=False, params=None, verbose=False) -> list:
    """Return NEW records with outliers removed from the selected props. OPT-IN — nothing is
    trimmed unless you call this; the returned list feeds straight into build_plan / autoplot /
    any combinator, so one call cleans every downstream plot.

    method : a method name for all selected props, or a dict {prop: method} (e.g.
             {"density": "mad", "mass": "iqr"}). See outlier_mask for methods/params.
    props  : which props to clean (default: every prop present). e.g. ["density"], ["mass","vol"].
    scope  : 'per_sample' (bounds from each sample's own values; default) or 'pooled' (bounds from
             every cell of that prop across all records — one global cutoff).
    paired : True for row-aligned records (load_ifxm_paired): a keep-mask is built from each
             selected prop, AND-combined, and applied to EVERY prop so a cell's props stay aligned.
    log    : compute bounds in log-space (for log-normal mass/volume).
    params : method params applied globally (merged over the method defaults). For different params
             per prop, call reject_outliers once per props subset.
    verbose: print how many cells/rows each sample dropped (rejection is never silent).

    Records are copied (never mutated); `meta` is preserved.
    """
    if not records:
        return records
    present = {k for r in records for k in r["props"]}
    sel = [p for p in (props if props is not None else present) if p in present]
    kw = dict(params or {})

    def meth(p):
        return method[p] if isinstance(method, dict) else method

    pooled_ref = {}
    if scope == "pooled":
        for p in sel:
            arrs = [np.asarray(r["props"].get(p, []), float) for r in records]
            arrs = [a for a in arrs if a.size]
            pooled_ref[p] = np.concatenate(arrs) if arrs else np.array([])

    out = []
    for r in records:
        newprops = dict(r["props"])
        if paired:
            length = next((np.asarray(r["props"][p], float).size for p in sel
                           if np.asarray(r["props"].get(p, []), float).size), None)
            if length is None:
                out.append({**r, "props": newprops})
                continue
            keep = np.ones(length, bool)
            for p in sel:
                a = np.asarray(r["props"].get(p, []), float)
                if a.size != length:
                    continue  # e.g. empty vol (unpaired mass-only sample) — don't break the joint mask
                ref = pooled_ref.get(p) if scope == "pooled" else None
                keep &= outlier_mask(a, meth(p), log=log, ref=ref, **kw)
            for p, a in r["props"].items():
                a = np.asarray(a, float)
                newprops[p] = a[keep] if a.size == length else a
            if verbose:
                print(f"  reject[{r['sample']}] paired: {length - int(keep.sum())}/{length} rows")
        else:
            for p in sel:
                a = np.asarray(r["props"].get(p, []), float)
                if a.size == 0:
                    continue
                ref = pooled_ref.get(p) if scope == "pooled" else None
                keep = outlier_mask(a, meth(p), log=log, ref=ref, **kw)
                newprops[p] = a[keep]
                if verbose:
                    print(f"  reject[{r['sample']}] {p} ({meth(p)}): "
                          f"{a.size - int(keep.sum())}/{a.size}")
        out.append({**r, "props": newprops})
    return out


# ---------------------------------------------------------------------------
# Role inference — classify each metadata column so plots can be chosen intelligently
# ---------------------------------------------------------------------------
DEFAULT_STRUCTURAL = {"sheet_name", "hdf5_key", "coulter_column", "calibration_factor"}
STRUCTURAL_PATTERNS = [re.compile(p) for p in (r"^has_", r"_gate_lower$", r"_gate_upper$",
                                               r"^bm_gate", r"^ifxm_gate")]
IDENTITY_COLS = {"sample_name"}

_TIME_NAME = re.compile(r"(?i)(^|_)(t|time|elapsed)(_|$)")
# suffix -> unit, longest first so '_hours' wins over '_h'
_TIME_SUFFIX = [("_hours", "h"), ("_hour", "h"), ("_hrs", "h"), ("_hr", "h"), ("_h", "h"),
                ("_minutes", "min"), ("_mins", "min"), ("_min", "min"),
                ("_seconds", "s"), ("_sec", "s")]
_UNIT_PER_HOUR = {"h": 1.0, "min": 60.0, "s": 3600.0}
_GRADIENT_HINTS = re.compile(
    r"(?i)(dose|conc|concentration|passage|day|cycle|generation|gen|temp|ph|dilution)")

_BOOL_SETS = [{"yes", "no"}, {"true", "false"}, {"0", "1"}, {"y", "n"}, {"t", "f"}]
_TRUTHY = {"yes", "true", "1", "y", "t"}
_ORDERED_MAX_CARD = 15   # numeric non-time non-hint below this is proposed as an ordered gradient


def _is_structural(col: str) -> bool:
    return col in DEFAULT_STRUCTURAL or any(p.search(col) for p in STRUCTURAL_PATTERNS)


def _time_unit(col: str):
    """(unit, confirm) for a time column name, or (None, None) if the name is not time-like.
    confirm=True when the unit had to be assumed (generic 'time'/'t' with no unit suffix)."""
    c = col.lower()
    for suf, u in _TIME_SUFFIX:
        if c.endswith(suf):
            return u, False
    if _TIME_NAME.search(col):
        return "h", True   # generic time name, assume hours but flag for confirmation
    return None, None


def _role(col, role, values, order, **kw):
    info = {"col": col, "role": role, "values": list(values), "order": list(order),
            "unit": None, "to_hours": None, "cardinality": len(set(values)),
            "use_for": set(), "reason": "", "confirm": False}
    info.update(kw)
    return info


def classify_column(col, values) -> dict:
    """Classify a single column from its NAME and its non-blank VALUES into a RoleInfo dict.
    Pure and dtype-agnostic (numeric-ness is inferred by float-coercion), so it works from either a
    metadata DataFrame column or the values in loaded records."""
    vals = [v for v in values if not _is_blank(v)]
    uniq = list(dict.fromkeys(vals))                 # unique, order-preserving
    low = {str(v).strip().lower() for v in uniq}
    nums = _to_floats(uniq)

    if _is_structural(col):
        return _role(col, "structural", uniq, uniq, reason="pipeline/structural column — ignored")
    if col in IDENTITY_COLS:
        return _role(col, "label", uniq, uniq, use_for={"label"}, reason="sample identity")

    # boolean (checkbox columns are literally 'yes'/'no')
    if low and any(low <= s for s in _BOOL_SETS) and len(low) <= 2:
        order = sorted(uniq, key=lambda v: str(v).strip().lower() in _TRUTHY)  # falsey -> truthy
        return _role(col, "boolean", uniq, order,
                     use_for={"group", "compare", "facet", "color", "series"},
                     reason=f"boolean ({'/'.join(map(str, order))})")

    # time
    unit, confirm = _time_unit(col)
    if nums is not None and unit is not None:
        per_h = _UNIT_PER_HOUR[unit]
        order = sorted(uniq, key=lambda v: float(v))
        return _role(col, "time", uniq, order, unit=unit, to_hours=(lambda v, k=per_h: float(v) / k),
                     use_for={"order", "series", "color"}, confirm=confirm,
                     reason=f"time in {unit}" + (" (unit assumed — confirm)" if confirm else ""))

    # ordered / gradient (numeric, non-time)
    if nums is not None:
        hinted = bool(_GRADIENT_HINTS.search(col))
        if hinted or len(uniq) <= _ORDERED_MAX_CARD:
            order = sorted(uniq, key=lambda v: float(v))
            why = "named like a gradient" if hinted else f"{len(uniq)} distinct numeric values"
            return _role(col, "ordered", uniq, order,
                         use_for={"order", "group", "compare", "color", "series"}, confirm=True,
                         reason=f"numeric gradient? ({why}) — confirm ordered vs categorical")
        return _role(col, "continuous", uniq, sorted(uniq, key=lambda v: float(v)),
                     use_for={"color"}, reason="continuous numeric — color/scatter axis only")

    # categorical (any cardinality; no cap)
    order = sorted(uniq, key=lambda v: (-vals.count(v), str(v)))   # frequency desc, then alpha
    note = "" if len(uniq) < len(vals) else " (one value per sample)"
    return _role(col, "categorical", uniq, order,
                 use_for={"group", "compare", "facet", "color", "series"},
                 reason=f"categorical, {len(uniq)} values{note}")


def infer_roles(records, *, overrides=None, skip=("rep",)) -> dict:
    """Classify every metadata column present across `records` into a RoleInfo. `overrides` is a
    {col: role_name} map that pins a column's role (e.g. force a numeric column to 'ordered' or
    'categorical' after the user confirms the plan). `skip` omits synthesized columns from the
    report (rep is always available for ordering/labels regardless)."""
    cols = []
    for r in records:
        for c in r.get("meta", {}):
            if c not in cols:
                cols.append(c)
    roles = {}
    for col in cols:
        if col in skip:
            continue
        info = classify_column(col, [rget(r, col) for r in records])
        ov = (overrides or {}).get(col)
        if ov and ov != info["role"]:
            info = classify_column(col, [rget(r, col) for r in records])  # recompute base
            info["role"] = ov
            info["confirm"] = False
            info["reason"] = f"role overridden to {ov}"
            info["use_for"] = {
                "boolean": {"group", "compare", "facet", "color", "series"},
                "categorical": {"group", "compare", "facet", "color", "series"},
                "ordered": {"order", "group", "compare", "color", "series"},
                "time": {"order", "series", "color"},
                "continuous": {"color"},
                "label": {"label"}, "structural": set(),
            }.get(ov, info["use_for"])
        roles[col] = info
    return roles


def _role_for(records, col, roles):
    """RoleInfo for `col`: from `roles` if given, else inferred on the fly from record values."""
    if roles and col in roles:
        return roles[col]
    return classify_column(col, [rget(r, col) for r in records])


# ---------------------------------------------------------------------------
# Ordering + labels
# ---------------------------------------------------------------------------
def group_order(records, col, roles=None) -> list:
    """Values of `col` present in `records`, in canonical plotting order (role-defined)."""
    role = _role_for(records, col, roles)
    present = list(dict.fromkeys(rget(r, col) for r in records if rget(r, col) is not None))
    ordered = [v for v in role["order"] if v in present]
    ordered += [v for v in present if v not in ordered]
    return ordered


def _sort_key(r, col, role):
    v = rget(r, col)
    if role["role"] in ("time", "ordered", "continuous"):
        try:
            return (0, float(v))
        except (TypeError, ValueError):
            return (1, 0.0)
    order = role["order"]
    return (order.index(v), 0.0) if v in order else (len(order), 0.0)


def sort_records(records, by, roles=None) -> list:
    """Sort records by a list of (col, ...) — each col ordered per its role (time/ordered numeric,
    else canonical categorical order). `by` may be a list of column names or (col, _) pairs."""
    cols = [b[0] if isinstance(b, (tuple, list)) else b for b in by]
    role_map = {c: _role_for(records, c, roles) for c in cols}

    def key(r):
        return tuple(_sort_key(r, c, role_map[c]) for c in cols)
    return sorted(records, key=key)


def _time_label(t) -> str:
    try:
        t = float(t)
    except (TypeError, ValueError):
        return str(t)
    return f"{int(t)}h" if t == int(t) else f"{t}h"


def value_label(col, value, roles=None, records=None) -> str:
    """Short label for a single value of `col`, role-aware (time -> '6h', boolean -> 'is_x=yes')."""
    role = roles.get(col) if roles and col in roles else (
        _role_for(records, col, None) if records is not None else {"role": "categorical", "unit": None})
    if role["role"] == "time":
        conv = role.get("to_hours")
        return _time_label(conv(value) if conv else value)
    if role["role"] == "boolean":
        return f"{col}={value}"
    return str(value)


def _detail_label(r, roles, time_col) -> str:
    """Per-sample row/box label inside a single group: time (if any) + sample name."""
    parts = []
    if time_col is not None and rget(r, time_col) is not None:
        parts.append(value_label(time_col, rget(r, time_col), roles))
    parts.append(r["sample"])
    return " ".join(p for p in parts if p)


def _compare_label(r, group_col, roles, time_col) -> str:
    """Per-sample label in a cross-group comparison: group value + time + sample name."""
    parts = [str(rget(r, group_col))]
    if time_col is not None and rget(r, time_col) is not None:
        parts.append(value_label(time_col, rget(r, time_col), roles))
    parts.append(r["sample"])
    return " ".join(p for p in parts if p)


def _leftover_label(r, roles, exclude) -> str:
    """Box-tick label holding only the annotation info NOT already conveyed elsewhere on the
    figure (the fixed title value, and/or the bold group separators) -- i.e. every remaining
    boolean/categorical/ordered/time column, in metadata-sheet order, joined as 'a | b'. Falls
    back to the sample name if every column got excluded. General principle: never repeat on an
    axis tick what the title or a bold separator already says."""
    parts = []
    for col, info in roles.items():
        if col in exclude or info["role"] in ("structural", "label"):
            continue
        v = rget(r, col)
        if v is None:
            continue
        parts.append(value_label(col, v, roles))
    return " | ".join(parts) if parts else r["sample"]


# ---------------------------------------------------------------------------
# Legend placement — a legend never sits on top of plotted data
# ---------------------------------------------------------------------------
# Inside-the-axes anchor points, tried in this order (conventional spots first). If every one of
# them lands on something, the legend goes just outside the axes on the right.
_LEGEND_INSIDE_LOCS = ("upper right", "upper left", "lower right", "lower left", "center right",
                       "center left", "upper center", "lower center", "center")
_LEGEND_OUTSIDE = ("center left", (1.0, 0.5))    # (loc, anchor in axes fractions)
_INK_TOL = 8            # per-channel distance from the background that counts as "something drawn"
_LEGEND_PAD_PX = 2      # clearance kept around the legend box


def _ink_mask(ax, leg, renderer) -> np.ndarray:
    """Boolean image (row 0 = top) of every pixel that differs from the empty background — lines,
    boxes, fills, scatter points, separators — rendered with `leg` hidden. Uses the axes facecolor
    inside the axes and the figure facecolor elsewhere."""
    fig = ax.figure
    was_visible = leg.get_visible()
    leg.set_visible(False)
    try:
        fig.draw(renderer)
    finally:
        leg.set_visible(was_visible)
    img = np.asarray(renderer.buffer_rgba(), dtype=np.int16)[..., :3]

    def differs(color):
        return (np.abs(img - np.array(to_rgb(color)) * 255).max(axis=-1) > _INK_TOL)

    ink = differs(fig.get_facecolor())
    if ax.patch.get_visible():
        h = img.shape[0]
        bb = ax.bbox
        y0, y1 = int(np.floor(h - bb.y1)), int(np.ceil(h - bb.y0))
        x0, x1 = int(np.floor(bb.x0)), int(np.ceil(bb.x1))
        ink[max(y0, 0):y1, max(x0, 0):x1] = differs(ax.get_facecolor())[max(y0, 0):y1,
                                                                        max(x0, 0):x1]
    return ink


def _legend_ink(leg, ink, renderer) -> int:
    """How many drawn-feature pixels fall inside the legend's box (plus a small margin)."""
    bb = leg.get_window_extent(renderer)
    h, w = ink.shape
    pad = _LEGEND_PAD_PX
    x0, x1 = max(int(np.floor(bb.x0)) - pad, 0), min(int(np.ceil(bb.x1)) + pad, w)
    y0, y1 = max(int(np.floor(h - bb.y1)) - pad, 0), min(int(np.ceil(h - bb.y0)) + pad, h)
    return int(ink[y0:y1, x0:x1].sum()) if x1 > x0 and y1 > y0 else 0


def _move_legend(ax, leg, loc, anchor=None) -> None:
    leg.set_loc(loc)
    if anchor:
        leg.set_bbox_to_anchor(anchor, transform=ax.transAxes)
    else:
        leg.set_bbox_to_anchor(None)


def place_legend(ax, handles=None, labels=None, **kwargs):
    """Create (or, if `ax` already has one and no handles/labels are given, re-place) `ax`'s legend
    so that it does not overlap any plotted feature (lines, boxes, fills, points, separators).

    The current/requested position is kept if it is already clear; otherwise the 9 in-axes anchor
    points are tried in turn, then just outside the axes on the right (which savefig's
    bbox_inches='tight' includes). Overlap is measured on a rendered copy of the figure, so it sees
    everything actually drawn so far — call it AFTER the data is drawn (the toolkit's own `_save`
    also re-checks every legend just before writing a PNG). Returns the Legend (None if `ax` has no
    labeled artists). `loc=` may be passed to choose the preferred first position; other kwargs go
    to `ax.legend` (defaults: frameon=False, fontsize=8). Needs matplotlib >= 3.8."""
    leg = ax.get_legend()
    if handles is not None or labels is not None or leg is None:
        if handles is None and labels is None and not ax.get_legend_handles_labels()[0]:
            return None
        kwargs.setdefault("frameon", False)
        kwargs.setdefault("fontsize", 8)
        kwargs["loc"] = "upper right" if kwargs.get("loc", "best") == "best" else kwargs["loc"]
        leg = ax.legend(handles=handles, labels=labels, **kwargs)
    if leg is None or not leg.get_visible():
        return leg

    fig = ax.figure
    renderer = RendererAgg(int(np.ceil(fig.bbox.width)), int(np.ceil(fig.bbox.height)), fig.dpi)
    ink = _ink_mask(ax, leg, renderer)
    if _legend_ink(leg, ink, renderer) == 0:
        return leg

    candidates = [(loc, None) for loc in _LEGEND_INSIDE_LOCS] + [_LEGEND_OUTSIDE]
    best, best_score = None, None
    for loc, anchor in candidates:
        _move_legend(ax, leg, loc, anchor)
        score = _legend_ink(leg, ink, renderer)
        if score == 0:
            return leg
        if best_score is None or score < best_score:
            best, best_score = (loc, anchor), score
    _move_legend(ax, leg, *best)
    print(f"  WARNING: legend overlaps plotted features at every position tried "
          f"({best_score} px at best) — shrink it or enlarge the figure")
    return leg


def _clear_legends(fig) -> None:
    """Safety net run by `_save`: re-place any legend that ended up on top of the data (e.g. the
    data was drawn after the legend, or a driver added its own `ax.legend`)."""
    for ax in fig.axes:
        if ax.get_legend() is not None:
            place_legend(ax)


# ---------------------------------------------------------------------------
# Low-level primitives (operate on a passed-in ax; no semantics)
# ---------------------------------------------------------------------------
def _run_separators(ax, keys, label_fn=None) -> None:
    """Bold labels under the axis with dark vertical separators between runs of equal `keys`."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    n = len(keys)
    # long rotated tick labels reach below the default label row; push the bold labels below them
    max_tick = max((len(t.get_text()) for t in ax.get_xticklabels()), default=0)
    y_lab = -0.30 - (0.12 if max_tick > 20 else 0.0)
    prev, start = object(), 0
    for i, k in enumerate(list(keys) + [object()]):
        if i == n or k != prev:
            if i > 0 and start < n:
                mid = (start + i - 1) / 2
                lab = label_fn(prev) if label_fn else str(prev)
                if i - start <= 2:      # narrow run: stack "a | b" so neighbours don't collide
                    lab = lab.replace(" | ", "\n")
                ax.text(mid, y_lab, lab, ha="center", va="top", transform=trans,
                        fontsize=9, fontweight="bold")
                if i < n:
                    ax.axvline(i - 0.5, color="black", lw=1.2, alpha=0.8, zorder=1)
            start, prev = i, k


def draw_ridge(ax, arrays, labels, colors, xlabel, overlap: float = 1.7) -> None:
    """Ridge plot of per-row histograms (shared bins, max-normalized), stacked top-down."""
    lo = min(v.min() for v in arrays)
    hi = max(v.max() for v in arrays)
    bins = np.linspace(lo, hi, 41)
    n = len(arrays)
    for i, (vals, c) in enumerate(zip(arrays, colors)):
        counts, _ = np.histogram(vals, bins=bins, density=True)
        h = counts / counts.max() * overlap if counts.max() > 0 else counts
        stair = np.concatenate([h, h[-1:]])
        base = n - 1 - i
        ax.fill_between(bins, base, base + stair, step="post", color=c, alpha=0.6, zorder=n - i)
        ax.step(bins, base + stair, where="post", color="black", lw=0.8, zorder=n - i)
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel(xlabel)


def draw_boxes(ax, arrays, labels, colors, ylabel, sep_keys=None, sep_label_fn=None) -> None:
    """Boxplot + jittered datapoints, one box per array. Optional run separators under the axis."""
    n = len(arrays)
    for i, (vals, c) in enumerate(zip(arrays, colors)):
        jitter = np.random.uniform(-0.18, 0.18, len(vals))
        ax.scatter(np.full(len(vals), i, float) + jitter, vals,
                   color=c, alpha=0.1, s=4, zorder=2, linewidths=0)
        ax.boxplot(vals, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor=c, alpha=0.5),
                   medianprops=dict(color="black", linewidth=1.5),
                   whiskerprops=dict(color=c), capprops=dict(color=c))
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
    if sep_keys is not None:
        _run_separators(ax, sep_keys, sep_label_fn)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylabel(ylabel)


def draw_ecdf(ax, arrays, labels, colors, xlabel) -> None:
    """Overlaid empirical CDFs, one line per array."""
    for vals, c, lab in zip(arrays, colors, labels):
        x = np.sort(np.asarray(vals, float))
        if x.size == 0:
            continue
        y = np.arange(1, x.size + 1) / x.size
        ax.plot(x, y, color=c, lw=1.6, label=str(lab))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cumulative fraction")
    place_legend(ax)


def draw_timecourse(ax, series: dict, colors: dict, xlabel, ylabel) -> None:
    """Per-replicate means as points + per-x average as a line, one series per key.
    `series`: {key: [(x, y), ...]}."""
    for sval, pts in sorted(series.items(), key=lambda kv: str(kv[0])):
        c = colors.get(sval, FALLBACK_COLOR)
        t = np.array([p[0] for p in pts], float)
        y = np.array([p[1] for p in pts], float)
        ax.scatter(t, y, color=c, s=30, alpha=0.7, zorder=3)
        uniq = sorted(set(t))
        avg = [y[t == u].mean() for u in uniq]
        ax.plot(uniq, avg, color=c, lw=1.5, marker="o", ms=6, zorder=2, label=str(sval))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    place_legend(ax)


def draw_scatter_marginal(fig, subplot_spec, x, y, color, xlabel, ylabel, title) -> None:
    """A per-cell scatter with marginal histograms (x on top, y on right), inside subplot_spec."""
    inner = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=subplot_spec, width_ratios=[3, 1], height_ratios=[1, 3],
        hspace=0.03, wspace=0.03)
    ax_sc = fig.add_subplot(inner[1, 0])
    ax_xh = fig.add_subplot(inner[0, 0], sharex=ax_sc)
    ax_yh = fig.add_subplot(inner[1, 1], sharey=ax_sc)
    fig.add_subplot(inner[0, 1]).axis("off")

    ax_sc.scatter(x, y, color=color, s=2, alpha=0.3, linewidths=0)
    ax_xh.hist(x, bins=30, color=color, alpha=0.7)
    ax_yh.hist(y, bins=30, orientation="horizontal", color=color, alpha=0.7)

    plt.setp(ax_xh.get_xticklabels(), visible=False)
    plt.setp(ax_yh.get_yticklabels(), visible=False)
    ax_xh.tick_params(bottom=False, labelsize=6)
    ax_yh.tick_params(left=False, labelsize=6)
    ax_xh.set_title(title, fontsize=8)
    ax_sc.set_xlabel(xlabel, fontsize=7)
    ax_sc.set_ylabel(ylabel, fontsize=7)
    ax_sc.tick_params(labelsize=6)


# Backward-compatible private aliases (older code referenced the underscore names).
_ridge, _boxes = draw_ridge, draw_boxes


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def _slug(value) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "-", str(value).strip().lower()).strip("-") or "na"


def _save(fig, name: str, out_dir) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _clear_legends(fig)
    fig.savefig(out_dir / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_dir.name}/{name}")


def _with_prop(records: list, prop: str) -> list:
    return [r for r in records if len(r["props"].get(prop, [])) > 0]


def _find_time_col(records, roles):
    """The first column whose role is 'time', if any (used to order samples within groups)."""
    roles = roles or infer_roles(records)
    for col, info in roles.items():
        if info["role"] == "time":
            return col
    return None


# ---------------------------------------------------------------------------
# Mid-level combinators (parameterized by WHICH column to group / compare / order by)
# ---------------------------------------------------------------------------
def plot_grouped(records, prop, ylabel, datatype, fig_dir, *, group_col, roles=None,
                 kinds=("ridge", "box"), time_col=None, colors=None) -> None:
    """DETAIL within each group: one figure per value of `group_col`; each ridge row / box is a
    sample in that group, ordered by time (if a time column exists) then rep.
    Files: {datatype}_{prop}_{kind}_{group_col}={slug(value)}.png"""
    recs = _with_prop(records, prop)
    if not recs:
        return
    roles = roles or infer_roles(records)
    time_col = time_col if time_col is not None else _find_time_col(records, roles)
    colors = colors or color_map_for(records, group_col, roles)
    for v in group_order(recs, group_col, roles):
        sub = [r for r in recs if rget(r, group_col) == v]
        order_by = [c for c in (time_col, "rep") if c]
        sub = sort_records(sub, order_by, roles) if order_by else sub
        if not sub:
            continue
        arrays = [r["props"][prop] for r in sub]
        col = colors.get(v, FALLBACK_COLOR)
        labels = [_detail_label(r, roles, time_col) for r in sub]
        title = f"{datatype} {prop} — {value_label(group_col, v, roles)}"
        tag = f"{group_col}={_slug(v)}"
        if "ridge" in kinds:
            fig, ax = plt.subplots(figsize=(8, max(3, len(sub) * 0.55)))
            draw_ridge(ax, arrays, labels, [col] * len(sub), ylabel)
            ax.set_title(title)
            _save(fig, f"{datatype}_{prop}_ridge_{tag}.png", fig_dir)
        if "box" in kinds:
            times = [rget(r, time_col) for r in sub] if time_col else None
            # time is already conveyed by the bold time separators below the axis, so the
            # per-box tick only needs whatever annotation columns are still unaccounted for.
            box_exclude = {group_col, time_col} if time_col else {group_col}
            box_labels = [_leftover_label(r, roles, box_exclude) for r in sub]
            fig, ax = plt.subplots(figsize=(max(5, len(sub) * 0.9), 5))
            draw_boxes(ax, arrays, box_labels, [col] * len(sub), ylabel,
                       sep_keys=times, sep_label_fn=_time_label if times else None)
            ax.set_title(title)
            _save(fig, f"{datatype}_{prop}_box_{tag}.png", fig_dir)


def compare_groups(records, prop, ylabel, datatype, fig_dir, *, group_col, roles=None,
                   kinds=("box", "ridge"), agg="per_sample", colors=None, time_col=None,
                   label_exclude=None) -> None:
    """COMPARISON across the values of `group_col`. agg='per_sample' (default) draws one box/ridge
    row per sample, colored by its group value and separated by group; agg='pool' pools all cells
    per group value into a single box/ridge row. Both a box and a ridge are produced by default.
    `label_exclude`: extra annotation columns to leave out of the per-box tick label beyond
    `group_col` itself (used by cross_groups, whose synthetic group_col already encodes them).
    Files: {datatype}_{prop}_{kind}_by_{group_col}.png"""
    recs = _with_prop(records, prop)
    if not recs:
        return
    roles = roles or infer_roles(records)
    time_col = time_col if time_col is not None else _find_time_col(records, roles)
    colors = colors or color_map_for(records, group_col, roles)
    values = group_order(recs, group_col, roles)
    if len(values) < 1:
        return

    if agg == "pool":
        arrays, labels, cols = [], [], []
        for v in values:
            pooled = np.concatenate([r["props"][prop] for r in recs if rget(r, group_col) == v])
            if len(pooled) == 0:
                continue
            arrays.append(pooled)
            labels.append(value_label(group_col, v, roles))
            cols.append(colors.get(v, FALLBACK_COLOR))
        sep_keys = None
        handles = None
        box_labels = labels
    else:  # per_sample
        sub = sort_records(recs, [c for c in (group_col, time_col, "rep") if c], roles)
        arrays = [r["props"][prop] for r in sub]
        labels = [_compare_label(r, group_col, roles, time_col) for r in sub]
        # group_col's value is already the bold separator below the axis, so the per-box
        # tick only needs whatever annotation columns are still unaccounted for.
        box_exclude = {group_col, *(label_exclude or ())}
        box_labels = [_leftover_label(r, roles, box_exclude) for r in sub]
        cols = [colors.get(rget(r, group_col), FALLBACK_COLOR) for r in sub]
        sep_keys = [rget(r, group_col) for r in sub]
        handles = [Patch(facecolor=colors.get(v, FALLBACK_COLOR),
                         label=value_label(group_col, v, roles)) for v in values]
    if not arrays:
        return

    title = f"{datatype} {prop} — by {group_col}"
    if "box" in kinds:
        fig, ax = plt.subplots(figsize=(max(5, len(arrays) * 0.9), 5))
        draw_boxes(ax, arrays, box_labels, cols, ylabel,
                   sep_keys=sep_keys, sep_label_fn=(lambda v: str(v)) if sep_keys else None)
        if handles:
            place_legend(ax, handles=handles)
        ax.set_title(title)
        _save(fig, f"{datatype}_{prop}_box_by_{group_col}.png", fig_dir)
    if "ridge" in kinds:
        fig, ax = plt.subplots(figsize=(9, max(3, len(arrays) * 0.55)))
        draw_ridge(ax, arrays, labels, cols, ylabel)
        if handles:
            place_legend(ax, handles=handles)
        ax.set_title(title)
        _save(fig, f"{datatype}_{prop}_ridge_by_{group_col}.png", fig_dir)


def timecourse_by(records, prop, ylabel, datatype, fig_dir, *, time_col=None, roles=None,
                  series_col=None, colors=None) -> None:
    """Timecourse of per-sample means vs a unit-aware time axis, colored by `series_col`.
    Files: {datatype}_{prop}_timecourse[_{series_col}].png"""
    recs = _with_prop(records, prop)
    if not recs:
        return
    roles = roles or infer_roles(records)
    time_col = time_col if time_col is not None else _find_time_col(records, roles)
    if time_col is None:
        return
    trole = _role_for(records, time_col, roles)
    conv = trole.get("to_hours") or (lambda v: float(v))
    colors = colors or (color_map_for(records, series_col, roles) if series_col else {})

    series = {}
    for r in recs:
        try:
            x = conv(rget(r, time_col))
        except (TypeError, ValueError):
            continue
        key = rget(r, series_col) if series_col else "all"
        series.setdefault(key, []).append((x, r["props"][prop].mean()))

    fig, ax = plt.subplots(figsize=(8, 5))
    draw_timecourse(ax, series, colors, "Time (h)", ylabel)
    suffix = f"_{series_col}" if series_col else ""
    ax.set_title(f"{datatype} {prop} — timecourse" + (f" by {series_col}" if series_col else ""))
    _save(fig, f"{datatype}_{prop}_timecourse{suffix}.png", fig_dir)


def scatter_by(records, prop_x, prop_y, xlabel, ylabel, datatype, fig_dir, *, group_col=None,
               color_col=None, roles=None, ncols=4) -> None:
    """Per-cell scatter of prop_y vs prop_x with marginal histograms — a grid of per-sample panels.
    One figure per value of `group_col` (or a single figure if group_col is None). Panels colored
    by `color_col` (defaults to group_col). Pass records from load_ifxm_paired so x and y are
    row-aligned. Files: {datatype}_{prop_y}_vs_{prop_x}[_{group_col}={slug(value)}].png"""
    roles = roles or infer_roles(records)
    color_col = color_col or group_col
    colors = color_map_for(records, color_col, roles) if color_col else {}

    def _xy(r):
        x = np.asarray(r["props"].get(prop_x, []), float)
        y = np.asarray(r["props"].get(prop_y, []), float)
        if x.size == 0 or y.size == 0:
            return None
        if x.size != y.size:
            print(f"  WARNING: {r['sample']} {prop_x}/{prop_y} lengths differ "
                  f"({x.size} vs {y.size}) — skipping (use load_ifxm_paired?)")
            return None
        keep = np.isfinite(x) & np.isfinite(y)
        return (x[keep], y[keep]) if keep.sum() else None

    usable = [(r, xy) for r in records if (xy := _xy(r)) is not None]
    if not usable:
        return
    time_col = _find_time_col(records, roles)

    if group_col is None:
        groups = [(None, usable)]
    else:
        groups = [(v, [ru for ru in usable if rget(ru[0], group_col) == v])
                  for v in group_order([r for r, _ in usable], group_col, roles)]

    xy_by_sample = {r["sample"]: xy for r, xy in usable}
    for gval, sub in groups:
        if not sub:
            continue
        order_by = [c for c in (time_col, "rep") if c]
        ordered_recs = sort_records([r for r, _ in sub], order_by, roles) if order_by \
            else [r for r, _ in sub]
        sub = [(r, xy_by_sample[r["sample"]]) for r in ordered_recs]

        n = len(sub)
        nc = min(n, ncols)
        nrows = int(np.ceil(n / nc))
        fig = plt.figure(figsize=(nc * 3.5, nrows * 3.5))
        outer = gridspec.GridSpec(nrows, nc, figure=fig, hspace=0.55, wspace=0.45)
        for idx, (r, (x, y)) in enumerate(sub):
            ri, ci = divmod(idx, nc)
            c = colors.get(rget(r, color_col), FALLBACK_COLOR) if color_col else FALLBACK_COLOR
            draw_scatter_marginal(fig, outer[ri, ci], x, y, c, xlabel, ylabel,
                                  _detail_label(r, roles, time_col))
        for idx in range(n, nrows * nc):
            ri, ci = divmod(idx, nc)
            fig.add_subplot(outer[ri, ci]).axis("off")
        gsuffix = f"_{group_col}={_slug(gval)}" if group_col is not None else ""
        gtitle = f" — {value_label(group_col, gval, roles)}" if group_col is not None else ""
        fig.suptitle(f"{datatype} {prop_y} vs {prop_x}{gtitle}", fontsize=11)
        _save(fig, f"{datatype}_{prop_y}_vs_{prop_x}{gsuffix}.png", fig_dir)


def facet(records, prop, ylabel, datatype, fig_dir, *, facet_col, roles=None, inner="box",
          ncols=4, colors=None) -> None:
    """Compact grid: one panel per value of `facet_col` in a SINGLE figure (each panel pools that
    value's cells into one inner box/ridge). File: {datatype}_{prop}_facet_{facet_col}.png"""
    recs = _with_prop(records, prop)
    if not recs:
        return
    roles = roles or infer_roles(records)
    colors = colors or color_map_for(records, facet_col, roles)
    values = group_order(recs, facet_col, roles)
    n = len(values)
    if n == 0:
        return
    nc = min(n, ncols)
    nrows = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nrows, nc, figsize=(nc * 3.2, nrows * 3.2), squeeze=False)
    for idx, v in enumerate(values):
        ax = axes[idx // nc][idx % nc]
        arr = np.concatenate([r["props"][prop] for r in recs if rget(r, facet_col) == v])
        c = colors.get(v, FALLBACK_COLOR)
        if inner == "ridge":
            draw_ridge(ax, [arr], [""], [c], ylabel)
        else:
            draw_boxes(ax, [arr], [""], [c], ylabel)
        ax.set_title(value_label(facet_col, v, roles), fontsize=9)
    for idx in range(n, nrows * nc):
        axes[idx // nc][idx % nc].axis("off")
    fig.suptitle(f"{datatype} {prop} — by {facet_col}", fontsize=11)
    _save(fig, f"{datatype}_{prop}_facet_{facet_col}.png", fig_dir)


def cross_groups(records, prop, ylabel, datatype, fig_dir, *, cols, roles=None,
                 kinds=("box", "ridge"), colors=None) -> None:
    """CROSSING on request: compare the cross-product of two (or more) columns. Adds a synthetic
    joined key (e.g. 'activated | DMEM') and runs compare_groups on it, per_sample.
    File: {datatype}_{prop}_{kind}_by_{colA}-x-{colB}.png"""
    cross_col = "-x-".join(cols)
    tagged = []
    for r in records:
        key = " | ".join(str(rget(r, c)) for c in cols)
        rr = dict(r)
        rr["meta"] = dict(r["meta"])
        rr["meta"][cross_col] = key
        tagged.append(rr)
    compare_groups(tagged, prop, ylabel, datatype, fig_dir, group_col=cross_col, roles=roles,
                   kinds=kinds, agg="per_sample", colors=colors, label_exclude=set(cols))


# ---------------------------------------------------------------------------
# Grid search: two columns varied jointly -> heatmap of per-sample means
# ---------------------------------------------------------------------------
_GRID_ROLES = ("boolean", "categorical", "ordered")


def suggest_grids(records, roles=None, *, min_coverage=0.75, exclude=()) -> list:
    """Detect pairs of annotation columns that look like a GRID SEARCH (two parameters varied
    together across samples), for grid_heatmap. Only a SUGGESTION — never plotted by default; the
    user must confirm before a pair goes into build_plan(grid_pairs=...).
    A pair (a, b) qualifies when both are boolean/categorical/ordered with >= 2 levels (>= 3 on at
    least one), the samples actually cross them (more distinct combos than either column has
    levels, so the two aren't 1:1 / confounded), and >= `min_coverage` of the a x b combos exist.
    Time columns are never proposed. The column with more levels goes on x (ties: metadata order).
    Returns [{cols:(x, y), shape:(nx, ny), filled, coverage, repeats:{(xv, yv): n}}], best first."""
    roles = roles or infer_roles(records)
    cands = [c for c, i in roles.items() if i["role"] in _GRID_ROLES and c not in exclude
             and len(group_order(records, c, roles)) >= 2]
    out = []
    for i, a in enumerate(cands):
        for b in cands[i + 1:]:
            recs = [r for r in records if rget(r, a) is not None and rget(r, b) is not None]
            na, nb = len(group_order(recs, a, roles)), len(group_order(recs, b, roles))
            if max(na, nb) < 3:
                continue
            counts = {}
            for r in recs:
                key = (rget(r, a), rget(r, b))
                counts[key] = counts.get(key, 0) + 1
            if len(counts) <= max(na, nb):
                continue
            coverage = len(counts) / (na * nb)
            if coverage < min_coverage:
                continue
            x, y, nx, ny = (a, b, na, nb) if na >= nb else (b, a, nb, na)
            repeats = {(k if x == a else k[::-1]): n for k, n in counts.items() if n > 1}
            out.append({"cols": (x, y), "shape": (nx, ny), "filled": len(counts),
                        "coverage": coverage, "repeats": repeats})
    return sorted(out, key=lambda s: (-s["coverage"], -s["filled"]))


def _grid_tick(col, value, roles) -> str:
    """Heatmap axis tick: numeric ordered values compactly ('0.0625', '0' not '0.0')."""
    if roles.get(col, {}).get("role") == "ordered":
        try:
            return f"{float(value):g}"
        except (TypeError, ValueError):
            pass
    return value_label(col, value, roles)


def _wrap_label(label, width=14) -> str:
    """Break a long sample/annotation label onto lines at '_', ' ', '-' or '|' so it fits a cell."""
    parts = re.split(r"(?<=[_\s|-])", str(label))
    lines, cur = [], ""
    for p in parts:
        if cur and len(cur) + len(p) > width:
            lines.append(cur)
            cur = ""
        cur += p
    return "\n".join(lines + [cur]).strip()


def _grid_fmt(values) -> str:
    """Enough decimals to resolve the spread of the plotted values (density ~4, mass/vol ~0-1)."""
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    spread = np.ptp(v) if v.size > 1 else (abs(v[0]) if v.size else 1.0)
    dec = int(np.clip(np.ceil(-np.log10(spread)) + 2, 0, 6)) if spread > 0 else 2
    return f"{{:.{dec}f}}"


def _annotate_cells(ax, arr, cmap, norm, fmt, notes=None) -> None:
    """Print each finite cell's value (plus an optional small note line) in a contrasting color."""
    for (i, j), v in np.ndenumerate(arr):
        if not np.isfinite(v):
            continue
        r, g, b, _ = cmap(norm(v))
        color = "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else "black"
        note = (notes or {}).get((i, j))
        ax.text(j, i, fmt.format(v) + (f"\n{note}" if note else ""), ha="center", va="center",
                fontsize=8 if note else 10, color=color)


def grid_heatmap(records, prop, label, datatype, fig_dir, *, cols, roles=None, agg="mean",
                 show_repeats=False, cmap="viridis", fmt=None) -> None:
    """GRID SEARCH over two jointly-varied columns: a heatmap with cols[0] on x and cols[1] on y.
    Each sample is reduced to one number (agg='mean' of its cells, or 'median'); each grid cell is
    that sample value, or — when several samples share the cell (typically a control condition
    re-run through the session) — the average of those samples, marked 'mean of n=k'. Missing
    combinations are grey. show_repeats=True adds a side panel on the SAME color scale showing
    each repeated sample individually, in run order (to check control drift over the session).
    Use only when the user confirmed the grid (see suggest_grids) — never added by default.
    File: {datatype}_{prop}_heatmap_{colX}-x-{colY}.png"""
    x_col, y_col = cols
    recs = [r for r in _with_prop(records, prop)
            if rget(r, x_col) is not None and rget(r, y_col) is not None]
    if not recs:
        return
    roles = roles or infer_roles(records)
    stat = {"mean": np.mean, "median": np.median}[agg]
    xs, ys = group_order(recs, x_col, roles), group_order(recs, y_col, roles)
    time_col = _find_time_col(records, roles)

    cell = {}                                       # (xv, yv) -> [records], in run order
    for r in (sort_records(recs, [time_col], roles) if time_col else recs):
        cell.setdefault((rget(r, x_col), rget(r, y_col)), []).append(r)
    grid = np.full((len(ys), len(xs)), np.nan)
    notes = {}
    for (xv, yv), rs in cell.items():
        i, j = ys.index(yv), xs.index(xv)
        grid[i, j] = np.mean([stat(r["props"][prop]) for r in rs])
        if len(rs) > 1:
            notes[(i, j)] = f"(mean of n={len(rs)})"
    repeated = [(k, rs) for k, rs in cell.items() if len(rs) > 1] if show_repeats else []

    rep_arr = None
    if repeated:
        width = max(len(rs) for _, rs in repeated)
        rep_arr = np.full((len(repeated), width), np.nan)
        for i, (_, rs) in enumerate(repeated):
            rep_arr[i, :len(rs)] = [stat(r["props"][prop]) for r in rs]
    shown = np.concatenate([grid.ravel()] + ([rep_arr.ravel()] if rep_arr is not None else []))
    finite = shown[np.isfinite(shown)]
    vmin, vmax = finite.min(), finite.max()
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    fmt = fmt or _grid_fmt(finite)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#e6e6e6")
    norm = plt.Normalize(vmin, vmax)

    grid_w = max(4.0, 1.25 * len(xs))
    h = max(3.5, 0.85 * len(ys) + 1.8)
    if rep_arr is not None:
        rep_w = max(2.5, 1.4 * rep_arr.shape[1])
        fig, (ax, axr) = plt.subplots(1, 2, figsize=(grid_w + rep_w + 2.5, h), gridspec_kw=dict(
            width_ratios=[grid_w, rep_w], wspace=0.35))
    else:
        fig, ax = plt.subplots(figsize=(grid_w + 2, h))
        axr = None

    im = ax.imshow(np.ma.masked_invalid(grid), cmap=cm, norm=norm, aspect="auto")
    _annotate_cells(ax, grid, cm, norm, fmt, notes)
    ax.set_xticks(range(len(xs)), [_grid_tick(x_col, v, roles) for v in xs])
    ax.set_yticks(range(len(ys)), [_grid_tick(y_col, v, roles) for v in ys])
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)

    if axr is not None:
        axr.imshow(np.ma.masked_invalid(rep_arr), cmap=cm, norm=norm, aspect="auto")
        # label each repeat by whatever annotation is left to tell them apart (else sample name):
        # on the x ticks when there's one repeated condition, else inside each cell
        rep_labels = [[_wrap_label(_leftover_label(r, roles, {x_col, y_col})) for r in rs]
                      for _, rs in repeated]
        ticks = [f"#{j + 1}" for j in range(rep_arr.shape[1])]
        if len(repeated) == 1:
            ticks = [f"{t}\n{lab}" for t, lab in zip(ticks, rep_labels[0])]
            _annotate_cells(axr, rep_arr, cm, norm, fmt)
        else:
            _annotate_cells(axr, rep_arr, cm, norm, fmt, {
                (i, j): lab for i, labs in enumerate(rep_labels) for j, lab in enumerate(labs)})
        axr.set_yticks(range(len(repeated)),
                       [f"{_grid_tick(x_col, xv, roles)} | {_grid_tick(y_col, yv, roles)}"
                        for (xv, yv), _ in repeated])
        axr.set_xticks(range(rep_arr.shape[1]), ticks, fontsize=7 if len(repeated) == 1 else None)
        axr.set_xlabel("repeat, in run order")
        axr.set_title("Repeated controls (individually)", fontsize=10)
        fig.colorbar(im, ax=[ax, axr], label=label, fraction=0.04, pad=0.03)
    else:
        fig.colorbar(im, ax=ax, label=label, fraction=0.05, pad=0.03)
    ax.set_title(f"{datatype} {prop} — {agg} per sample, {x_col} × {y_col}", fontsize=10)
    _save(fig, f"{datatype}_{prop}_heatmap_{x_col}-x-{y_col}.png", fig_dir)


# ---------------------------------------------------------------------------
# High-level: infer a plot plan, show it, execute it
# ---------------------------------------------------------------------------
def build_plan(records, datatype, *, roles=None, props=None, scatter_pairs=None,
               grid_pairs=None, include_ordered=True) -> dict:
    """Infer a plot plan (list of PlotSpecs) from the column roles. For every boolean / categorical
    (and, if include_ordered, approved ordered) column: a plot_grouped + a compare_groups per prop.
    If a time column exists: a timecourse per prop, split by each grouping column. Scatters are
    added for scatter_pairs. `props` is filtered to those non-empty in at least one record, so a
    mass-only / volume-only experiment proposes no no-op plots for absent properties.
    Grid heatmaps are added ONLY for the user-confirmed `grid_pairs` — each an (x_col, y_col) tuple
    or {"cols": (x, y), "show_repeats": True, ...} (extra keys go to grid_heatmap). Grid-like
    column pairs found by suggest_grids but not in grid_pairs are listed under
    'grid_suggestions' for render_plan to offer — never plotted by default. Returns
    {roles, props, plots:[{fn,prop,kwargs,rationale}], grid_suggestions}."""
    roles = roles or infer_roles(records)
    present = {k for r in records for k, v in r["props"].items() if len(v) > 0}
    props = [p for p in (props if props is not None else present) if p in present]
    time_col = _find_time_col(records, roles)

    group_cols = [c for c, i in roles.items()
                  if i["role"] in ("boolean", "categorical")
                  or (include_ordered and i["role"] == "ordered")]

    plots = []
    for prop in props:
        for gc in group_cols:
            plots.append({"fn": "plot_grouped", "prop": prop, "kwargs": {"group_col": gc},
                          "rationale": f"per-{gc} detail"})
            plots.append({"fn": "compare_groups", "prop": prop, "kwargs": {"group_col": gc},
                          "rationale": f"compare across {gc}"})
        if time_col:
            plots.append({"fn": "timecourse_by", "prop": prop,
                          "kwargs": {"time_col": time_col, "series_col": None},
                          "rationale": f"timecourse over {time_col}"})
            for gc in group_cols:
                plots.append({"fn": "timecourse_by", "prop": prop,
                              "kwargs": {"time_col": time_col, "series_col": gc},
                              "rationale": f"timecourse over {time_col}, split by {gc}"})
    for (px, py, xl, yl) in (scatter_pairs or []):
        if px not in present or py not in present:
            continue  # e.g. a mass-only experiment has no density/volume to scatter against
        gc = group_cols[0] if group_cols else None
        plots.append({"fn": "scatter_by", "prop": f"{py}_vs_{px}",
                      "kwargs": {"prop_x": px, "prop_y": py, "xlabel": xl, "ylabel": yl,
                                 "group_col": gc},
                      "rationale": "per-cell scatter" + (f" per {gc}" if gc else "")})
    confirmed = set()
    for gp in (grid_pairs or []):
        kw = dict(gp) if isinstance(gp, dict) else {"cols": tuple(gp)}
        kw["cols"] = tuple(kw["cols"])
        confirmed.add(frozenset(kw["cols"]))
        for prop in props:
            plots.append({"fn": "grid_heatmap", "prop": prop, "kwargs": kw,
                          "rationale": f"grid heatmap {kw['cols'][0]} × {kw['cols'][1]}"})
    suggestions = [s for s in suggest_grids(records, roles)
                   if frozenset(s["cols"]) not in confirmed]
    return {"roles": roles, "props": props, "plots": plots, "grid_suggestions": suggestions}


def render_plan(plan) -> str:
    """Human-readable summary of inferred roles + proposed plots, with any confirm-me flags —
    the text to show the user for approval before generating/executing the driver."""
    roles = plan["roles"]
    lines = ["Inferred metadata roles:"]
    for col, i in roles.items():
        flag = "  [CONFIRM]" if i["confirm"] else ""
        unit = f" [{i['unit']}]" if i.get("unit") else ""
        lines.append(f"  - {col}: {i['role']}{unit} — {i['reason']}{flag}")
    confirms = [c for c, i in roles.items() if i["confirm"]]
    if confirms:
        lines.append("")
        lines.append("Needs your confirmation: " + ", ".join(confirms)
                     + "  (pass overrides={col: 'ordered'|'categorical'|...} to pin)")
    lines.append("")
    lines.append(f"Proposed plots ({len(plan['plots'])}):")
    for p in plan["plots"]:
        lines.append(f"  - {p['fn']}({p['prop']}) — {p['rationale']}")
    sugg = plan.get("grid_suggestions") or []
    if sugg:
        lines.append("")
        lines.append("Suggested grid heatmaps — NOT included; ask the user, add to GRID_PAIRS "
                     "only if confirmed:")
        for s in sugg:
            (x, y), (nx, ny) = s["cols"], s["shape"]
            rep = ""
            if s["repeats"]:
                rep = "; repeated: " + ", ".join(
                    f"{_grid_tick(x, a, roles)} | {_grid_tick(y, b, roles)} ×{n}"
                    for (a, b), n in s["repeats"].items())
            lines.append(f"  - {x} × {y}: {nx}×{ny} grid, {s['filled']}/{nx * ny} combos "
                         f"filled{rep}")
    return "\n".join(lines)


_COMBINATORS = None


def autoplot(records, plan, datatype, fig_dir, prop_labels=None, paired_records=None) -> None:
    """Execute a plan's PlotSpecs. Distribution/timecourse plots run on `records`; scatter_by specs
    run on `paired_records` (row-aligned, from load_ifxm_paired) — pass it whenever the plan has
    scatter pairs, or those specs are skipped. `prop_labels`: {prop: axis_label} for y-axis labels."""
    global _COMBINATORS
    if _COMBINATORS is None:
        _COMBINATORS = {"plot_grouped": plot_grouped, "compare_groups": compare_groups,
                        "timecourse_by": timecourse_by, "scatter_by": scatter_by,
                        "facet": facet, "grid_heatmap": grid_heatmap}
    labels = prop_labels or {}
    roles = plan["roles"]
    warned = False
    for spec in plan["plots"]:
        fn = _COMBINATORS[spec["fn"]]
        prop = spec["prop"]
        kw = dict(spec["kwargs"], roles=roles)
        if spec["fn"] == "scatter_by":
            if paired_records is None:
                if not warned:
                    print("  (skipping scatter specs — pass paired_records=load_ifxm_paired(...))")
                    warned = True
                continue
            fn(paired_records, datatype=datatype, fig_dir=fig_dir, **kw)
        else:
            fn(records, prop, labels.get(prop, prop), datatype, fig_dir, **kw)


# ---------------------------------------------------------------------------
# PowerPoint export
# ---------------------------------------------------------------------------
def save_pptx(fig_dir, out_path) -> None:
    """Compile every PNG in fig_dir into a 16:9 deck, one image per slide (centered, aspect-fit).
    The deck name should start with the analysis dir's date (e.g. '2026-09-22_<exp>_figures.pptx')
    — a warning is printed if it doesn't."""
    from PIL import Image as PILImage
    from pptx import Presentation
    from pptx.util import Inches, Emu

    fig_dir, out_path = Path(fig_dir), Path(out_path)
    if not re.match(r"\d{4}-?\d{2}-?\d{2}", out_path.name):
        print(f"  WARNING: {out_path.name} has no date prefix — set EXP_NAME to the full analysis "
              f"dir name (e.g. '2026-09-22_<exp>')")
    prs = Presentation()
    slide_w, slide_h = Inches(13.33), Inches(7.5)
    prs.slide_width, prs.slide_height = slide_w, slide_h
    blank_layout = prs.slide_layouts[6]

    pngs = sorted(fig_dir.glob("*.png"))
    for png in pngs:
        img_w, img_h = PILImage.open(png).size
        if img_w / img_h > slide_w / slide_h:
            w, h = slide_w, Emu(int(slide_w * img_h / img_w))
        else:
            h, w = slide_h, Emu(int(slide_h * img_w / img_h))
        left = Emu(int((slide_w - w) // 2))
        top  = Emu(int((slide_h - h) // 2))
        slide = prs.slides.add_slide(blank_layout)
        slide.shapes.add_picture(str(png), left=left, top=top, width=w, height=h)

    prs.save(str(out_path))
    print(f"  saved {out_path.name} ({len(pngs)} slides)")
