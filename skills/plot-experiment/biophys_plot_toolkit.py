"""
biophys_plot_toolkit — reusable plotting library for compiled FXM/SMR/Coulter experiments.

Ported from the hand-written reference analysis
  C:/code/scratch/biophys_data_analysis/20260629_fl5_wt-gfpgem_ifxm_analysis/analyze_ifxm_data.py
and generalized so a short per-experiment driver can load compiled/annotated data and produce
the standard figure grid without re-deriving any of the plotting internals.

The whole toolkit is organized around a uniform *record* abstraction:
    {sample, cond, drug, time_h, rep, props: {name: np.ndarray}}
Every loader emits records of this shape; every plot function consumes them, so plotting code is
agnostic to whether the data came from Coulter or iFXM.

Data conventions (see references/data_schema.md for the full schema):
  Coulter — one property "volume" (fL), gated upstream.
  iFXM    — mass (matched_mass, pg), density (buoyant_density + baseline_density),
            vol_cal (volume_fL), vol_uncal (volume). Mass gated by bm_gate; the three iFXM
            properties share a single mask computed on the *uncalibrated* volume.

NOTE on baseline_density: the absolute density baseline is variable between experiments and is
NOT stored in any data file, so `load_ifxm` REQUIRES it as an explicit argument. There is no
default — passing nothing raises, by design, so an experiment is never silently plotted against
the wrong baseline.
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from pathlib import Path

# ---------------------------------------------------------------------------
# Styling defaults (override in the driver if an experiment needs different keys)
# ---------------------------------------------------------------------------
COND_COLORS = {
    "activated":     "#4878d0",
    "starved":       "#ee854a",
    "drug_treated":  "#6acc65",
    "proliferating": "#956cb4",
}
DRUG_COLORS = {  # keys match the lowercased output of _norm_drug
    "dmso":       "#4878d0",
    "1um-wnk463": "#ee854a",
    "2um-zt-1a":  "#6acc65",
}
FALLBACK_COLOR = "#aaaaaa"

# (prop_key, axis_label) lists — defaults matching the reference experiments.
COULTER_PROPS = [("volume", "Volume (fL)")]
IFXM_PROPS = [
    ("mass",      "Buoyant mass (pg)"),
    ("density",   "Density (g/mL)"),
    ("vol_cal",   "Calibrated volume (fL)"),
    ("vol_uncal", "Volume (fL)"),
]

# Sentinel so a forgotten baseline_density fails loudly rather than silently defaulting.
_REQUIRED = object()


# ---------------------------------------------------------------------------
# Metadata normalization
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


# ---------------------------------------------------------------------------
# Robust metadata access — annotation column NAMES vary between experiments, so
# every column is looked up by a (configurable) name with a graceful fallback.
# ---------------------------------------------------------------------------
def _get(row, col, default=""):
    """Value of `row[col]`, or `default` if the column is absent/blank."""
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
    """A gate bound as float, or `default` (±inf) when the column is missing or NaN.
    This makes a missing/ungated experiment mean 'no cutoff' instead of dropping every cell."""
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


# Colorblind-friendly cycle used to auto-assign colors to unknown condition/drug values.
_AUTO_PALETTE = [
    "#4878d0", "#ee854a", "#6acc65", "#956cb4", "#d65f5f", "#82c6e2",
    "#d5bb67", "#8c613c", "#dc7ec0", "#797979", "#b8c9d0", "#5975a4",
]


def build_color_map(values, base=None):
    """Return a color map covering every value in `values`.

    Known keys in `base` keep their color; any new/unknown value (a condition or drug name you
    invented for this experiment) gets the next distinct palette color, so arbitrary annotation
    values render in different colors instead of all-gray. Order-stable and deterministic.
    """
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


# ---------------------------------------------------------------------------
# Outlier rejection (for visualization only)
# ---------------------------------------------------------------------------
def _reject_3sigma(a: np.ndarray) -> np.ndarray:
    """Drop points outside mean +/- 3 sigma (applied to the volume properties)."""
    a = a[np.isfinite(a)]
    if a.size == 0:
        return a
    mu, sd = a.mean(), a.std()
    return a[(a >= mu - 3 * sd) & (a <= mu + 3 * sd)]


def _reject_mad(a: np.ndarray, thresh: float = 3.5) -> np.ndarray:
    """Drop points by modified z-score |0.6745*(x-median)/MAD| > thresh.
    Robust to the heavy tails in the density data (applied to density only)."""
    a = a[np.isfinite(a)]
    if a.size == 0:
        return a
    med = np.median(a)
    mad = np.median(np.abs(a - med))
    if mad == 0:
        return a
    return a[np.abs(0.6745 * (a - med) / mad) <= thresh]


# Mask-returning variants (keep two arrays row-aligned when trimming a scatter axis).
def _mask_3sigma(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    finite = np.isfinite(a)
    if finite.sum() == 0:
        return finite
    mu, sd = a[finite].mean(), a[finite].std()
    return finite & (a >= mu - 3 * sd) & (a <= mu + 3 * sd)


def _mask_mad(a: np.ndarray, thresh: float = 3.5) -> np.ndarray:
    a = np.asarray(a, float)
    finite = np.isfinite(a)
    if finite.sum() == 0:
        return finite
    med = np.median(a[finite])
    mad = np.median(np.abs(a[finite] - med))
    if mad == 0:
        return finite
    return finite & (np.abs(0.6745 * (a - med) / mad) <= thresh)


def _trim_mask(a: np.ndarray, how) -> np.ndarray:
    """Keep-mask for outlier trimming on one scatter axis. how in {None, 'mad', '3sigma'}."""
    if how is None:
        return np.isfinite(np.asarray(a, float))
    if how == "mad":
        return _mask_mad(a)
    if how == "3sigma":
        return _mask_3sigma(a)
    raise ValueError(f"unknown trim mode {how!r} (use None, 'mad', or '3sigma')")


# ---------------------------------------------------------------------------
# Loaders  ->  list of records {sample, cond, drug, time_h, rep, props:{name: arr}}
# ---------------------------------------------------------------------------
def load_coulter(data_dir, *, sample_col="sample_name", key_col="h5_key",
                 condition_col="condition", time_col="time_h", drug_col="drug_name",
                 normalize_cond=_norm_cond, normalize_drug=_norm_drug) -> list:
    """Load Coulter single-cell volumes from <data_dir>/coulter/{metadata.csv, data.h5}.

    Column names are configurable because annotation schemas vary between experiments. Only
    `sample_col` is required; if `key_col` is absent it falls back to the sample column, and a
    missing `condition_col` / `time_col` / `drug_col` degrades gracefully (condition -> unnamed
    group, time -> 0.0, drug -> none). Pass `normalize_cond`/`normalize_drug=lambda x: x` to
    disable the default value fixups.

    Each /data/{key} table's first column is the per-cell volume (fL) array.
    """
    d = Path(data_dir) / "coulter"
    h5 = d / "data.h5"
    meta = pd.read_csv(d / "metadata.csv")
    _require_col(meta, sample_col, "sample-name")
    key_col = key_col if key_col in meta.columns else sample_col

    recs = []
    for _, r in meta.iterrows():
        arr = pd.read_hdf(h5, f"/data/{r[key_col]}").iloc[:, 0].to_numpy()
        arr = _reject_3sigma(arr)  # volume plot -> 3-sigma outlier rejection
        sample = _get(r, sample_col, "")
        recs.append({
            "sample": sample,
            "cond":   normalize_cond(_get(r, condition_col, "")),
            "drug":   normalize_drug(_get(r, drug_col, "")),
            "time_h": _get_float(r, time_col, 0.0),
            "rep":    _rep(sample),
            "props":  {"volume": arr},
        })
    return recs


def load_ifxm(data_dir, baseline_density=_REQUIRED, *, sample_col="sample_name",
              key_col="hdf5_key", condition_col="condition", time_col="time_h",
              drug_col="drug_name", bm_lower_col="bm_gate_lower", bm_upper_col="bm_gate_upper",
              ifxm_lower_col="ifxm_gate_lower", ifxm_upper_col="ifxm_gate_upper",
              normalize_cond=_norm_cond, normalize_drug=_norm_drug) -> list:
    """Load paired iFXM data from <data_dir>/ifxm/experiment_data.h5.

    baseline_density (REQUIRED): fluid baseline added to measured buoyant_density to get the
    absolute density (g/mL). It varies per experiment and is not stored in the data, so it must
    be supplied explicitly.

    Annotation/gate column names are configurable (schemas vary between experiments). Only
    `sample_col` is required; `key_col` falls back to the sample column if absent, and any
    missing annotation or gate column degrades gracefully (missing gate -> unbounded, i.e. no
    cutoff, rather than dropping every cell).

    Per sample, reads /samples/{key}/pairing (matched_mass, buoyant_density, volume) and
    /samples/{key}/volume_calibrated (volume_fL). Samples missing either subtable are skipped
    (e.g. a no-iFXM proliferating control). Mass is gated by the bm gate; the three iFXM
    properties share one mask computed on the uncalibrated volume (cells are row-aligned).
    """
    if baseline_density is _REQUIRED:
        raise ValueError(
            "load_ifxm requires baseline_density (g/mL) — it varies per experiment and is not "
            "stored in the data files. Set it explicitly at the top of your driver."
        )

    h5 = Path(data_dir) / "ifxm" / "experiment_data.h5"
    meta = pd.read_hdf(h5, "/metadata")
    _require_col(meta, sample_col, "sample-name")
    key_col = key_col if key_col in meta.columns else sample_col

    recs = []
    with pd.HDFStore(h5, "r") as store:
        for _, r in meta.iterrows():
            base = f"/samples/{r[key_col]}"
            if f"{base}/pairing" not in store or f"{base}/volume_calibrated" not in store:
                continue  # e.g. proliferating_culture_noifxm has no paired iFXM data
            pr = store[f"{base}/pairing"]
            vc = store[f"{base}/volume_calibrated"]

            mass = pr["matched_mass"].to_numpy()
            dens = pr["buoyant_density"].to_numpy() + baseline_density
            vun  = pr["volume"].to_numpy()
            vcal = vc["volume_fL"].to_numpy()

            bm_lo = _gate_bound(r, bm_lower_col, -np.inf)
            bm_hi = _gate_bound(r, bm_upper_col,  np.inf)
            ix_lo = _gate_bound(r, ifxm_lower_col, -np.inf)
            ix_hi = _gate_bound(r, ifxm_upper_col,  np.inf)
            bm_mask = np.isfinite(mass) & (mass >= bm_lo) & (mass <= bm_hi)
            # one logical mask from the uncalibrated volume, shared across the paired properties
            ifxm_mask = np.isfinite(vun) & (vun >= ix_lo) & (vun <= ix_hi)

            def clean(a, m):
                a = a[m]
                return a[np.isfinite(a)]

            sample = _get(r, sample_col, "")
            recs.append({
                "sample": sample,
                "cond":   normalize_cond(_get(r, condition_col, "")),
                "drug":   normalize_drug(_get(r, drug_col, "")),
                "time_h": _get_float(r, time_col, 0.0),
                "rep":    _rep(sample),
                "props": {
                    "mass":      clean(mass, bm_mask),
                    # density plot -> MAD-based outlier rejection for visualization
                    "density":   _reject_mad(clean(dens, ifxm_mask)),
                    # volume plots -> additional 3-sigma outlier rejection
                    "vol_cal":   _reject_3sigma(clean(vcal, ifxm_mask)),
                    "vol_uncal": _reject_3sigma(clean(vun, ifxm_mask)),
                },
            })
    return recs


def load_ifxm_paired(data_dir, baseline_density=_REQUIRED, *, sample_col="sample_name",
                     key_col="hdf5_key", condition_col="condition", time_col="time_h",
                     drug_col="drug_name", bm_lower_col="bm_gate_lower",
                     bm_upper_col="bm_gate_upper", ifxm_lower_col="ifxm_gate_lower",
                     ifxm_upper_col="ifxm_gate_upper",
                     normalize_cond=_norm_cond, normalize_drug=_norm_drug) -> list:
    """Like load_ifxm, but keeps per-cell arrays row-ALIGNED across properties (no per-property
    outlier rejection, which would desync lengths). Use this for scatter_2d so a cell's mass,
    density and volume stay paired. One shared finite+gate mask is applied to every property.
    Takes the same configurable column names as load_ifxm.

    Returns records whose props hold equal-length arrays for: mass, density, vol_cal, vol_uncal.
    """
    if baseline_density is _REQUIRED:
        raise ValueError("load_ifxm_paired requires baseline_density (g/mL). Set it in your driver.")

    h5 = Path(data_dir) / "ifxm" / "experiment_data.h5"
    meta = pd.read_hdf(h5, "/metadata")
    _require_col(meta, sample_col, "sample-name")
    key_col = key_col if key_col in meta.columns else sample_col

    recs = []
    with pd.HDFStore(h5, "r") as store:
        for _, r in meta.iterrows():
            base = f"/samples/{r[key_col]}"
            if f"{base}/pairing" not in store or f"{base}/volume_calibrated" not in store:
                continue
            pr = store[f"{base}/pairing"]
            vc = store[f"{base}/volume_calibrated"]

            mass = pr["matched_mass"].to_numpy()
            dens = pr["buoyant_density"].to_numpy() + baseline_density
            vun  = pr["volume"].to_numpy()
            vcal = vc["volume_fL"].to_numpy()

            bm_lo = _gate_bound(r, bm_lower_col, -np.inf)
            bm_hi = _gate_bound(r, bm_upper_col,  np.inf)
            ix_lo = _gate_bound(r, ifxm_lower_col, -np.inf)
            ix_hi = _gate_bound(r, ifxm_upper_col,  np.inf)
            # single mask keeps every property the same length and row-aligned per cell
            mask = (np.isfinite(mass) & np.isfinite(dens) & np.isfinite(vun) & np.isfinite(vcal)
                    & (mass >= bm_lo) & (mass <= bm_hi)
                    & (vun >= ix_lo) & (vun <= ix_hi))

            sample = _get(r, sample_col, "")
            recs.append({
                "sample": sample,
                "cond":   normalize_cond(_get(r, condition_col, "")),
                "drug":   normalize_drug(_get(r, drug_col, "")),
                "time_h": _get_float(r, time_col, 0.0),
                "rep":    _rep(sample),
                "props": {
                    "mass":      mass[mask],
                    "density":   dens[mask],
                    "vol_cal":   vcal[mask],
                    "vol_uncal": vun[mask],
                },
            })
    return recs


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def _time_label(t: float) -> str:
    return f"{int(t)}h" if float(t) == int(t) else f"{t}h"


def _save(fig, name: str, out_dir) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_dir.name}/{name}")


def _with_prop(records: list, prop: str) -> list:
    return [r for r in records if len(r["props"].get(prop, [])) > 0]


# ---------------------------------------------------------------------------
# Plot primitives (operate on a passed-in ax)
# ---------------------------------------------------------------------------
def _timepoint_separators(ax, times: list, n: int) -> None:
    """Bold per-timepoint labels under the axis with dark vertical separators."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    prev, start = None, 0
    for i, t in enumerate(list(times) + [None]):
        if t != prev:
            if prev is not None:
                mid = (start + i - 1) / 2
                ax.text(mid, -0.30, _time_label(prev), ha="center", va="top",
                        transform=trans, fontsize=9, fontweight="bold")
                if i < n:
                    ax.axvline(i - 0.5, color="black", lw=1.2, alpha=0.8, zorder=1)
            start, prev = i, t


def _ridge(ax, arrays: list, labels: list, colors: list, xlabel: str,
           overlap: float = 1.7) -> None:
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
        ax.fill_between(bins, base, base + stair, step="post", color=c,
                        alpha=0.6, zorder=n - i)
        ax.step(bins, base + stair, where="post", color="black", lw=0.8, zorder=n - i)
    ax.set_yticks([n - 1 - i for i in range(n)])
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel(xlabel)


def _boxes(ax, arrays: list, labels: list, colors: list, ylabel: str,
           times: list = None) -> None:
    """Boxplot + jittered datapoints, one box per array."""
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
    if times is not None:
        _timepoint_separators(ax, times, n)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylabel(ylabel)


def _timecourse(ax, records: list, prop: str, color_map: dict, key: str,
                ylabel: str) -> None:
    """Per-replicate means as points + per-timepoint average as a line, one series/key."""
    series = {}
    for r in records:
        arr = r["props"].get(prop)
        if arr is None or len(arr) == 0:
            continue
        series.setdefault(r[key], []).append((r["time_h"], arr.mean()))

    for sval, pts in sorted(series.items()):
        c = color_map.get(sval, FALLBACK_COLOR)
        t = np.array([p[0] for p in pts])
        y = np.array([p[1] for p in pts])
        ax.scatter(t, y, color=c, s=30, alpha=0.7, zorder=3)
        uniq = sorted(set(t))
        avg = [y[t == u].mean() for u in uniq]
        ax.plot(uniq, avg, color=c, lw=1.5, marker="o", ms=6, zorder=2, label=str(sval))
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=8)


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------
def _full_label(r: dict) -> str:
    """Time + drug (if any) + rep — so every row shows both timepoint and drug."""
    parts = [_time_label(r["time_h"])]
    if r["drug"]:
        parts.append(r["drug"])
    parts.append(r["rep"])
    return " ".join(parts)


def _box_tick(r: dict) -> str:
    """Drug (if any) + rep; the timepoint is carried by the separators."""
    return f"{r['drug']} {r['rep']}" if r["drug"] else r["rep"]


def _ridge_label(r: dict) -> str:
    # starved samples share time 0, so show full names to tell them apart
    return r["sample"] if r["cond"] == "starved" else _full_label(r)


def _box_label(r: dict) -> str:
    return r["sample"] if r["cond"] == "starved" else _box_tick(r)


# ---------------------------------------------------------------------------
# Plot drivers (create the figure, delegate to a primitive, title, save)
# ---------------------------------------------------------------------------
def ridge_box_by_condition(records, prop, ylabel, datatype, fig_dir,
                           cond_colors=None) -> None:
    """One ridge figure and one box figure per condition."""
    recs = _with_prop(records, prop)
    cond_colors = cond_colors or build_color_map(sorted(set(r["cond"] for r in recs)), COND_COLORS)
    for cond in sorted(set(r["cond"] for r in recs)):
        sub = sorted([r for r in recs if r["cond"] == cond],
                     key=lambda r: (r["time_h"], r["rep"]))
        if not sub:
            continue
        arrays = [r["props"][prop] for r in sub]
        colors = [cond_colors.get(cond, FALLBACK_COLOR)] * len(sub)
        times  = [r["time_h"] for r in sub]

        fig, ax = plt.subplots(figsize=(8, max(3, len(sub) * 0.55)))
        _ridge(ax, arrays, [_ridge_label(r) for r in sub], colors, ylabel)
        ax.set_title(f"{datatype} {prop} — {cond}")
        _save(fig, f"{datatype}_{prop}_ridge_{cond}.png", fig_dir)

        fig, ax = plt.subplots(figsize=(max(5, len(sub) * 0.9), 5))
        _boxes(ax, arrays, [_box_label(r) for r in sub], colors, ylabel, times=times)
        ax.set_title(f"{datatype} {prop} — {cond}")
        _save(fig, f"{datatype}_{prop}_box_{cond}.png", fig_dir)


def timecourse(records, prop, ylabel, datatype, fig_dir, cond_colors=None) -> None:
    """One timecourse figure combining all conditions (series colored by condition)."""
    recs = _with_prop(records, prop)
    if not recs:
        return
    cond_colors = cond_colors or build_color_map(sorted(set(r["cond"] for r in recs)), COND_COLORS)
    fig, ax = plt.subplots(figsize=(8, 5))
    _timecourse(ax, recs, prop, cond_colors, "cond", ylabel)
    ax.set_title(f"{datatype} {prop} — timecourse")
    _save(fig, f"{datatype}_{prop}_timecourse.png", fig_dir)


def drug_split(records, prop, ylabel, datatype, fig_dir, drug_colors=None) -> None:
    """Ridge / box / timecourse for the drug-treated arm, split & colored by drug_name.
    No-op when there are no drug-treated records (e.g. a wt experiment without a drug arm)."""
    # Trigger on any record carrying a drug annotation, regardless of the condition label —
    # so this works whether the drug arm is called "drug_treated" or something experiment-specific.
    recs = [r for r in _with_prop(records, prop) if r["drug"]]
    if not recs:
        return
    drug_colors = drug_colors or build_color_map(
        sorted(set(r["drug"] for r in recs)), DRUG_COLORS)
    sub = sorted(recs, key=lambda r: (r["drug"], r["time_h"], r["rep"]))
    arrays = [r["props"][prop] for r in sub]
    colors = [drug_colors.get(r["drug"], FALLBACK_COLOR) for r in sub]
    labels = [_full_label(r) for r in sub]  # time + drug + rep
    handles = [Patch(facecolor=drug_colors[d], label=d)
               for d in drug_colors if any(r["drug"] == d for r in sub)]

    fig, ax = plt.subplots(figsize=(9, max(3, len(sub) * 0.55)))
    _ridge(ax, arrays, labels, colors, ylabel)
    ax.legend(handles=handles, frameon=False, fontsize=8)
    ax.set_title(f"{datatype} {prop} — drug treated")
    _save(fig, f"{datatype}_{prop}_drug_ridge.png", fig_dir)

    fig, ax = plt.subplots(figsize=(max(5, len(sub) * 0.9), 5))
    _boxes(ax, arrays, labels, colors, ylabel)
    ax.legend(handles=handles, frameon=False, fontsize=8)
    ax.set_title(f"{datatype} {prop} — drug treated")
    _save(fig, f"{datatype}_{prop}_drug_box.png", fig_dir)

    fig, ax = plt.subplots(figsize=(8, 5))
    _timecourse(ax, recs, prop, drug_colors, "drug", ylabel)
    ax.set_title(f"{datatype} {prop} — drug timecourse")
    _save(fig, f"{datatype}_{prop}_drug_timecourse.png", fig_dir)


# ---------------------------------------------------------------------------
# 2-D property-vs-property scatter with marginal histograms (per-cell)
# ---------------------------------------------------------------------------
def scatter_2d(records, prop_x, prop_y, xlabel, ylabel, datatype, fig_dir,
               cond_colors=None, drug_colors=None, trim_x=None, trim_y=None,
               ncols=4) -> None:
    """Per-cell scatter of prop_y vs prop_x with marginal histograms — one figure per condition,
    laid out as a grid of per-sample panels. Each panel is a main scatter with a histogram of
    prop_x above it and a horizontal histogram of prop_y to its right (nested GridSpec with
    3:1 / 1:3 ratios). Panels are colored by drug when present, else by condition, and titled per
    sample. Styled to match density_mass_scatter in the reference analysis.

    IMPORTANT: pass records from load_ifxm_paired (not load_ifxm) so prop_x and prop_y are
    row-aligned per cell and equal length. Records where either property is empty or the two
    lengths differ are skipped (with a warning) rather than mis-paired.

    trim_x / trim_y ('mad' | '3sigma' | None): reject outliers on that axis, applied jointly so
    the pair stays aligned. Use trim_y='mad' for density (heavy tails), mirroring the reference.

    Files: {datatype}_{prop_y}_vs_{prop_x}_{cond}.png
    """
    cond_colors = cond_colors or build_color_map(
        sorted(set(r["cond"] for r in records)), COND_COLORS)
    drug_colors = drug_colors or build_color_map(
        sorted(set(r["drug"] for r in records if r["drug"])), DRUG_COLORS)

    def _xy(r):
        x = np.asarray(r["props"].get(prop_x, []), float)
        y = np.asarray(r["props"].get(prop_y, []), float)
        if x.size == 0 or y.size == 0:
            return None
        if x.size != y.size:
            print(f"  WARNING: {r['sample']} {prop_x}/{prop_y} lengths differ "
                  f"({x.size} vs {y.size}) — skipping (did you use load_ifxm_paired?)")
            return None
        keep = _trim_mask(x, trim_x) & _trim_mask(y, trim_y)
        if keep.sum() == 0:
            return None
        return x[keep], y[keep]

    usable = [(r, xy) for r in records if (xy := _xy(r)) is not None]
    if not usable:
        return

    for cond in sorted(set(r["cond"] for r, _ in usable)):
        sub = sorted([ru for ru in usable if ru[0]["cond"] == cond],
                     key=lambda ru: (ru[0]["time_h"], ru[0]["drug"], ru[0]["rep"]))
        if not sub:
            continue

        n = len(sub)
        nc = min(n, ncols)
        nrows = int(np.ceil(n / nc))
        fig = plt.figure(figsize=(nc * 3.5, nrows * 3.5))
        outer = gridspec.GridSpec(nrows, nc, figure=fig, hspace=0.55, wspace=0.45)

        for idx, (r, (x, y)) in enumerate(sub):
            ri, ci = divmod(idx, nc)
            inner = gridspec.GridSpecFromSubplotSpec(
                2, 2, subplot_spec=outer[ri, ci],
                width_ratios=[3, 1], height_ratios=[1, 3],
                hspace=0.03, wspace=0.03,
            )
            ax_sc = fig.add_subplot(inner[1, 0])
            ax_xh = fig.add_subplot(inner[0, 0], sharex=ax_sc)
            ax_yh = fig.add_subplot(inner[1, 1], sharey=ax_sc)
            fig.add_subplot(inner[0, 1]).axis("off")

            c = drug_colors.get(r["drug"], FALLBACK_COLOR) if r["drug"] \
                else cond_colors.get(cond, FALLBACK_COLOR)

            ax_sc.scatter(x, y, color=c, s=2, alpha=0.3, linewidths=0)
            ax_xh.hist(x, bins=30, color=c, alpha=0.7)
            ax_yh.hist(y, bins=30, orientation="horizontal", color=c, alpha=0.7)

            plt.setp(ax_xh.get_xticklabels(), visible=False)
            plt.setp(ax_yh.get_yticklabels(), visible=False)
            ax_xh.tick_params(bottom=False)
            ax_yh.tick_params(left=False)
            ax_xh.set_title(_full_label(r), fontsize=8)
            ax_sc.set_xlabel(xlabel, fontsize=7)
            ax_sc.set_ylabel(ylabel, fontsize=7)
            ax_sc.tick_params(labelsize=6)
            ax_xh.tick_params(labelsize=6)
            ax_yh.tick_params(labelsize=6)

        for idx in range(n, nrows * nc):
            ri, ci = divmod(idx, nc)
            fig.add_subplot(outer[ri, ci]).axis("off")

        fig.suptitle(f"{datatype} {prop_y} vs {prop_x} — {cond}", fontsize=11)
        _save(fig, f"{datatype}_{prop_y}_vs_{prop_x}_{cond}.png", fig_dir)


# ---------------------------------------------------------------------------
# PowerPoint export
# ---------------------------------------------------------------------------
def save_pptx(fig_dir, out_path) -> None:
    """Compile every PNG in fig_dir into a 16:9 deck, one image per slide (centered, aspect-fit).
    Imported lazily so the toolkit is usable for plotting even without python-pptx installed."""
    from PIL import Image as PILImage
    from pptx import Presentation
    from pptx.util import Inches, Emu

    fig_dir, out_path = Path(fig_dir), Path(out_path)
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
