# Source manifest

The release was finalized on June 25, 2026. Checksums are SHA-256.

## Public release files

| File | SHA-256 |
| --- | --- |
| `main_plot_rsn_bifurcation_distributions_model_based.m` | `41db0d01926a3ea7319d7e5a66d043da00acbbc97c5624e9f06baac44b5d994c` |
| `main_plot_rsn_bifurcation_median_summaries_model_based.m` | `bd5db625f457a83fd18839db72b008a726e5c710e81a16046fc66feba188480d` |
| `main_plot_global_connectivity_distributions_model_free.m` | `63a0560c36bde7c21c866a0c98cf5dd931ee859a4b51eeaaeaa95851cac61c07` |
| `data/model_bifurcation_parameters.mat` | `e065a838d3aed2af35eebc037e36cf52ff397337ace4f7ebe6dfd3cb644a26a2` |
| `data/subject_global_connectivity_aligned.csv` | `7245beb0648efa86afb6474101efb7ed5e9b5dabab5fab0e33b37ed98f679f04` |
| `third_party/violinplot/Violin.m` | `f206af9e23a71f20cd2052df6c3f4d1eccc93b0bd1c58245d36e343163d12ebb` |
| `third_party/violinplot/violinplot.m` | `5d2d4dda7970ffaaff9dc461b3b9084848282cbfb903de0762f6f8802dc127fe` |

## Model-based source mapping

| Release filename | Original working filename |
| --- | --- |
| `main_plot_rsn_bifurcation_distributions_model_based.m` | `GA_plot_RSN_subcortical_boxplots_all_v2.m` |
| `main_plot_rsn_bifurcation_median_summaries_model_based.m` | `aux_figure_paper_subcortical_v2.m` |

The plotting operations are preserved. The release files are runnable scripts
with portable paths to `data/` and `third_party/violinplot/`.

## Model-free provenance

`main_plot_global_connectivity_distributions_model_free.m` is a MATLAB
translation of the Python aligned-paired-test analysis for the two selected
outlier-trimmed global-connectivity figures.

The derived release table was produced from `run-01` using Tian S1 subcortex
followed by Schaefer-200 cortex, linear detrending, order-2 Butterworth
`0.008-0.08 Hz` filtering, removal of 2 leading and 3 trailing filtered
samples, Pearson FC over 234 nodes, and group-to-whole-brain connectivity
excluding node self-connections.

MATLAB model-free statistics match the Python reference tables with maximum
absolute difference below `6e-15`.

## Non-release material

The following local folders are intentionally ignored by Git:

- `si_not_for_release/`
- `deprecated_not_for_release/`
- `reference_not_for_release/`

The reference derivation reproduces the released subject-level connectivity
table from private parcellated time series with maximum absolute difference
`1.33e-15`.
