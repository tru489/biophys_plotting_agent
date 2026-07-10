"""
Reference plotting driver for a compiled/annotated biophysics experiment.

This is the TEMPLATE the plot-experiment skill adapts per experiment. When generating a driver:
  - copy biophys_plot_toolkit.py next to this driver in the analysis directory (so the import
    below resolves locally and the analysis stays self-contained / reproducible),
  - set DATA_DIR, FIG_DIR and BASELINE_DENSITY for this experiment,
  - trim the property/condition/scatter lists to what the metadata actually contains.

Run it with the same Python env that has: numpy, pandas, matplotlib, tables, python-pptx, Pillow.
"""
from pathlib import Path

import biophys_plot_toolkit as tk

# ---------------------------------------------------------------------------
# Per-experiment configuration  (EDIT THESE)
# ---------------------------------------------------------------------------
# Compiled/annotated data root: expects <DATA_DIR>/coulter/{metadata.csv,data.h5}
# and <DATA_DIR>/ifxm/experiment_data.h5
DATA_DIR = Path(r"C:\path\to\<exp>_data")          # Windows
# DATA_DIR = Path("/Users/you/.../<exp>_data")     # macOS

FIG_DIR  = DATA_DIR.parent / f"{DATA_DIR.name.replace('_data', '')}_fig"
PPTX_OUT = DATA_DIR.parent / f"{DATA_DIR.name.replace('_data', '')}_figures.pptx"

# REQUIRED, per-experiment. Fluid baseline (g/mL) added to measured buoyant density to get the
# absolute density. It is NOT stored in any data file and varies between experiments — set it
# deliberately. (The FL5 reference experiments used 1.008.)
BASELINE_DENSITY = 1.008

# Which properties to plot. Trim to what this experiment actually has.
COULTER_PROPS = tk.COULTER_PROPS            # [("volume", "Volume (fL)")]
IFXM_PROPS    = tk.IFXM_PROPS               # mass / density / vol_cal / vol_uncal

# 2-D per-cell scatters with marginal histograms (iFXM only), one figure per condition.
# (prop_x, prop_y, xlabel, ylabel, trim_y)  — trim_y='mad' tames density's heavy tails.
SCATTER_PAIRS = [
    ("mass",    "density", "Buoyant mass (pg)",      "Total density (g/cm^3)", "mad"),
    ("vol_cal", "mass",    "Calibrated volume (fL)", "Buoyant mass (pg)",      None),
    ("vol_cal", "density", "Calibrated volume (fL)", "Total density (g/cm^3)", "mad"),
]


def main() -> None:
    print(f"=== plotting {DATA_DIR.name} ===")
    # If this experiment's annotation columns are named differently, pass the mapping, e.g.:
    #   tk.load_ifxm(DATA_DIR, baseline_density=BASELINE_DENSITY,
    #                condition_col="treatment", time_col="t_hours", drug_col="compound")
    # Missing columns / gates degrade gracefully; unknown condition/drug values auto-color.
    coulter = tk.load_coulter(DATA_DIR)
    ifxm    = tk.load_ifxm(DATA_DIR, baseline_density=BASELINE_DENSITY)
    # row-aligned variant for scatter (keeps a cell's props paired; no per-prop outlier drops)
    ifxm_paired = tk.load_ifxm_paired(DATA_DIR, baseline_density=BASELINE_DENSITY)
    print(f"  coulter samples: {len(coulter)} | ifxm samples: {len(ifxm)}")

    for prop, ylabel in COULTER_PROPS:
        tk.ridge_box_by_condition(coulter, prop, ylabel, "coulter", FIG_DIR)
        tk.timecourse(coulter, prop, ylabel, "coulter", FIG_DIR)
        tk.drug_split(coulter, prop, ylabel, "coulter", FIG_DIR)

    for prop, ylabel in IFXM_PROPS:
        tk.ridge_box_by_condition(ifxm, prop, ylabel, "ifxm", FIG_DIR)
        tk.timecourse(ifxm, prop, ylabel, "ifxm", FIG_DIR)
        tk.drug_split(ifxm, prop, ylabel, "ifxm", FIG_DIR)

    for px, py, xl, yl, trim_y in SCATTER_PAIRS:
        tk.scatter_2d(ifxm_paired, px, py, xl, yl, "ifxm", FIG_DIR, trim_y=trim_y)

    tk.save_pptx(FIG_DIR, PPTX_OUT)


if __name__ == "__main__":
    main()
