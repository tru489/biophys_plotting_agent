# <EXP_NAME> — plotting notes

<!-- Written by the plot-experiment skill next to <EXP_NAME>_figures.pptx. Keep it short: only what a
future run needs to recapitulate this analysis on a new dataset. Delete sections that don't apply. -->

## Inputs
- Data type(s): iFXM / Coulter (which halves existed)
- iFXM: `<...>_compiled/experiment_data.xlsx` — N samples (P paired, M mass-only, V volume-only)
- Coulter: `<...>_coulter_sample_annotation/` — N samples
- Driver: `<driver>.py` (+ copied `biophys_plot_toolkit.py`), plugin version X.Y.Z

## Cells per sample
<!-- Paste the driver's printed render_cell_counts table(s) (same as the deck's first slide). Note
any sample with has_gate = False (compiled ungated) or an unusually large drop. -->
| Sample | Cells from | Before gating | After gating | After outlier rejection: mass | … |
|---|---|---|---|---|---|

## Annotation columns → roles
| column | role | levels / range | notes |
|--------|------|----------------|-------|
| `time_h` | time | 0, 4, 24 h | orders samples; timecourse x-axis |
| `media` | categorical | RPMI, 250 mOsm, 400 mOsm | control = RPMI |
| `dose_uM` | ordered | 0, 0.1, 1, 10 | override: was `[CONFIRM]` |

- `ROLE_OVERRIDES`: `{...}` (and why)
- Ignored / label-only columns worth knowing about:
- Design: which columns are crossed (grid) vs. independent; repeated controls; missing combos

## Plots made
- Properties: mass / density / vol / coulter volume
- Grouping/compare axes: …; timecourse series: …; scatter pairs: …
- `GRID_PAIRS`: `[...]` (show_repeats?) — or none
- Cross-products (`cross_groups`), facets, hand-rolled figures: …
- Props-vs-time deck (`<EXP_NAME>_props_vs_time.pptx`): samples, props, trend — or skipped (why)

## Customizations & data handling
- Outlier rejection: none / `reject_outliers(...)` spec
- Axis limits, palettes, excluded samples, label tweaks, anything non-default

## To reuse on a new dataset
- Columns the new dataset must have (or their equivalents) and how to map renamed ones
- Anything that was experiment-specific and should be re-decided rather than copied
