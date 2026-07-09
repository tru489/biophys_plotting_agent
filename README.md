# biophys_plotting_agent

A Claude Code plugin for automated plotting of compiled biophysics experiments (integrated
FXM / SMR / Coulter — fluorescence-exclusion microscopy, suspended-microchannel-resonator buoyant
mass, and Coulter-counter volume).

The plugin ships a skill, **`plot-experiment`**, that:

1. inspects a compiled/annotated experiment's schema (conditions, drug treatments, time points),
2. generates a short, self-contained plotting driver on top of a bundled toolkit, and
3. runs it to produce the standard figure grid — ridge histograms, box+jitter plots, timecourse
   scatters, and per-cell property-vs-property scatters (e.g. volume vs mass / density) — plus a
   PowerPoint deck.

You then fine-tune the generated driver in normal Claude Code.

## Install

On any machine:

```
/plugin marketplace add tru489/biophys_plotting_agent
/plugin install biophys-plotting@biophys-tools
```

The plugin is copied into the local plugin cache at install (no per-use fetch). Invoke the skill
as `/biophys-plotting:plot-experiment`, or just ask Claude to "plot my experiment" from a data
directory.

Python dependencies (for running generated drivers) are in [environment.yaml](environment.yaml)
(`numpy pandas matplotlib pytables python-pptx pillow`):

```
conda env create -f environment.yaml
conda activate biophys_plotting
```

Or reuse any existing analysis env that has those packages.

## Input data

The skill plots **compiled + annotated** output of the `biophys_helpers` pipeline, one experiment
per `<exp>_data/` directory:

```
<exp>_data/
├── coulter/{metadata.csv, data.h5}      # Coulter single-cell volumes + annotations
└── ifxm/experiment_data.h5              # paired iFXM mass/volume/density + annotations
```

See [skills/plot-experiment/references/data_schema.md](skills/plot-experiment/references/data_schema.md)
for the full schema. Conditions, drugs, and time points are hand-annotated columns and vary per
experiment; the skill reads them from the metadata rather than assuming a fixed set.

> **Note:** absolute density needs a per-experiment `baseline_density` (g/mL) that is *not* stored
> in any data file. The skill asks for it; there is no silent default.

## Layout

```
biophys_plotting_agent/
├── .claude-plugin/{plugin.json, marketplace.json}
├── skills/plot-experiment/
│   ├── SKILL.md                     # the workflow
│   ├── biophys_plot_toolkit.py      # reusable plotting library
│   ├── reference_driver.py          # driver template the skill adapts
│   └── references/data_schema.md    # h5/metadata schema
└── environment.yaml                 # conda env for running generated drivers
```
