# Jhana whole-brain modeling paper code

MATLAB code and derived data for selected model-based and model-free figures
from the Jhana whole-brain modeling paper.

Raw participant time series, empirical FC matrices, model-fitting code, and
private analyses are not included in the GitHub release.

## Requirements

- MATLAB R2023b or newer
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Public analysis scripts

The repository root contains three public analysis files:

- `main_plot_rsn_bifurcation_distributions_model_based.m`
- `main_plot_rsn_bifurcation_median_summaries_model_based.m`
- `main_plot_global_connectivity_distributions_model_free.m`

Each file is a plain MATLAB script. Run it by name from MATLAB:

```matlab
main_plot_rsn_bifurcation_distributions_model_based
main_plot_rsn_bifurcation_median_summaries_model_based
main_plot_global_connectivity_distributions_model_free
```

The scripts display the figures without writing generated files. Calculated
statistics remain available in the MATLAB workspace.

## Data included in the release

The `data/` folder contains:

- `model_bifurcation_parameters.mat`
- `subject_global_connectivity_aligned.csv`

The MAT file contains the `40 x 10 x 8` modeled bifurcation-parameter array
used by both model-based scripts. The CSV contains derived subject-level
global network connectivity for 20 subjects, 10 conditions, and 8 model-space
node groups. It does not contain original time series.

The model-free code reproduces complete-case selection, Tukey `1.5 x IQR`
outlier removal, paired t-tests, Cohen's dz, exact two-sided Wilcoxon
signed-rank tests, and Benjamini-Hochberg FDR correction. Its statistics
match the original Python outputs to floating-point precision.

## Supporting folders

- `data/`: releasable modeled and derived data
- `third_party/violinplot/`: bundled violin-plot helper and BSD license
- `SOURCE_MANIFEST.md`: provenance and checksums

The local folders `si_not_for_release/`, `deprecated_not_for_release/`, and
`reference_not_for_release/` are intentionally excluded through `.gitignore`.

## Reproducibility notes

- Numerical values and reported statistics are deterministic.
- Horizontal point jitter is seeded for the model-free plots.
- The plotting and analysis operations are preserved from the working code.
- The public entry points are scripts with portable data and dependency paths.

## Licensing

No project-wide license has been selected. The bundled violin-plot helper
retains its BSD 3-Clause license in `third_party/violinplot/LICENSE`.
