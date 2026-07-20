"""
Reference plotting driver for a compiled/annotated biophysics experiment.

This is the TEMPLATE the plot-experiment skill adapts per experiment. When generating a driver:
  - copy biophys_plot_toolkit.py next to this driver in the analysis directory (so the import
    below resolves locally and the analysis stays self-contained / reproducible),
  - set COMPILED_DIR / COULTER_DIR, FIG_DIR and BASELINE_DENSITY for this experiment,
  - review the printed plan and set ROLE_OVERRIDES for any [CONFIRM] columns.

The toolkit reads whatever hand-added annotation columns exist and classifies each into a role
(boolean / categorical / time / ordered / continuous / label / structural). `build_plan` proposes
plots from those roles; `autoplot` runs them. Every column becomes a grouping/comparison axis; a
time column is used to order samples sequentially. Use the explicit combinators (see the bottom of
main) for full control.

Run it with the same Python env that has: numpy, pandas, matplotlib, openpyxl, python-pptx, Pillow.
"""
from pathlib import Path

import biophys_plot_toolkit as tk

# ---------------------------------------------------------------------------
# Per-experiment configuration  (EDIT THESE)
# ---------------------------------------------------------------------------
# The loaders read the raw biophys_helpers outputs directly (no reorg step). Point each at the
# timestamped output dir it needs (set one to None if that half of the experiment doesn't exist):
#   iFXM   -> compile_experiment.py's '<...>_compiled/' dir      (holds experiment_data.xlsx)
#   Coulter-> annotate_coulter_samples.py's '<...>_coulter_sample_annotation/' dir
EXP_NAME     = "<exp>"
COMPILED_DIR = Path(r"C:\path\to\<...>_compiled")                    # iFXM (or None)
COULTER_DIR  = Path(r"C:\path\to\<...>_coulter_sample_annotation")   # Coulter (or None)
# COMPILED_DIR = Path("/Users/you/.../<...>_compiled")              # macOS
# COULTER_DIR  = Path("/Users/you/.../<...>_coulter_sample_annotation")

OUT_ROOT = COMPILED_DIR.parent if COMPILED_DIR else COULTER_DIR.parent
FIG_DIR  = OUT_ROOT / f"{EXP_NAME}_fig"
PPTX_OUT = OUT_ROOT / f"{EXP_NAME}_figures.pptx"

# REQUIRED, per-experiment. Fluid baseline (g/mL) added to measured buoyant density to get the
# absolute density. It is NOT stored in any data file and varies between experiments — set it
# deliberately. (The FL5 reference experiments used 1.008.)
BASELINE_DENSITY = 1.008

# Axis labels for each property (also selects which props to plot). Trim to what exists.
IFXM_LABELS    = dict(tk.IFXM_PROPS)        # mass / density / vol_cal / vol_uncal
COULTER_LABELS = dict(tk.COULTER_PROPS)     # volume

# 2-D per-cell scatters with marginal histograms (iFXM only). (prop_x, prop_y, xlabel, ylabel)
SCATTER_PAIRS = [
    ("mass",    "density", "Buoyant mass (pg)",      "Total density (g/cm^3)"),
    ("vol_cal", "mass",    "Calibrated volume (fL)", "Buoyant mass (pg)"),
    ("vol_cal", "density", "Calibrated volume (fL)", "Total density (g/cm^3)"),
]

# Pin any column's role after reviewing the printed plan (resolves every [CONFIRM] flag).
# e.g. {"dose_uM": "ordered", "plate": "categorical"}. Leave empty to accept the inferred roles.
ROLE_OVERRIDES = {}


def main() -> None:
    print(f"=== plotting {EXP_NAME} ===")
    coulter = tk.load_coulter(COULTER_DIR) if COULTER_DIR else []
    ifxm = ifxm_paired = []
    if COMPILED_DIR:
        ifxm = tk.load_ifxm(COMPILED_DIR, baseline_density=BASELINE_DENSITY)
        # row-aligned variant for scatter (keeps a cell's props paired under one shared gate mask)
        ifxm_paired = tk.load_ifxm_paired(COMPILED_DIR, baseline_density=BASELINE_DENSITY)
    print(f"  coulter samples: {len(coulter)} | ifxm samples: {len(ifxm)}")

    # --- iFXM: inspect the annotation schema, show the plan, run it -----------------------
    if ifxm:
        roles = tk.infer_roles(ifxm, overrides=ROLE_OVERRIDES)
        plan = tk.build_plan(ifxm, "ifxm", roles=roles, props=list(IFXM_LABELS),
                             scatter_pairs=SCATTER_PAIRS)
        print(tk.render_plan(plan))       # review; set ROLE_OVERRIDES for any [CONFIRM] columns
        tk.autoplot(ifxm, plan, "ifxm", FIG_DIR, prop_labels=IFXM_LABELS,
                    paired_records=ifxm_paired)

    # --- Coulter --------------------------------------------------------------------------
    if coulter:
        croles = tk.infer_roles(coulter, overrides=ROLE_OVERRIDES)
        cplan = tk.build_plan(coulter, "coulter", roles=croles, props=list(COULTER_LABELS))
        print(tk.render_plan(cplan))
        tk.autoplot(coulter, cplan, "coulter", FIG_DIR, prop_labels=COULTER_LABELS)

    # --- Explicit combinators (comment out autoplot above and drive it yourself) ----------
    # Each groups/compares/orders by ANY column name; roles= is optional (inferred if omitted).
    #   tk.plot_grouped(ifxm, "mass", "Buoyant mass (pg)", "ifxm", FIG_DIR, group_col="is_activated")
    #   tk.compare_groups(ifxm, "mass", "Buoyant mass (pg)", "ifxm", FIG_DIR, group_col="is_activated")
    #   tk.timecourse_by(ifxm, "mass", "Buoyant mass (pg)", "ifxm", FIG_DIR, series_col="is_activated")
    #   tk.scatter_by(ifxm_paired, "mass", "density", "Mass", "Density", "ifxm", FIG_DIR,
    #                 group_col="is_activated")
    #   tk.facet(ifxm, "mass", "Buoyant mass (pg)", "ifxm", FIG_DIR, facet_col="media")
    #   tk.cross_groups(ifxm, "mass", "Buoyant mass (pg)", "ifxm", FIG_DIR,
    #                   cols=("is_activated", "media"))   # crossing on request

    tk.save_pptx(FIG_DIR, PPTX_OUT)


if __name__ == "__main__":
    main()
