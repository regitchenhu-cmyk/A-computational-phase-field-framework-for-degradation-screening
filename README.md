# Supplementary Software and Data

This archive accompanies the revised numerical study of coupled degradation
screening and contact-aware leakage post-processing in an elastomer seal. It
contains the MATLAB implementation, regression and study QA, portable
processed data, and the Python source used for the nine manuscript figures.

## Scope

The reported evidence is simulation-only. The supplied checks establish
implementation consistency and numerical sensitivity within the tested model
and parameter set. They do not validate a probability of failure, crack-growth
life, leakage onset, a replacement threshold, or absolute service durability.

## Quick start

Run commands from the extracted archive root. Paths in this archive are
relative; no drive-letter mapping is required.

1. Inspect the already generated results without running MATLAB:
   `processed_data/` contains the canonical tables, `figures/manuscript/`
   contains the submitted figures, `module_verification/outputs/` contains
   the isolated analytical checks, and
   `Results_Paper_ScientificallyCorrected_20260823/M7_Strain_Quantification_20260828/`
   contains the endpoint-strain audit.
2. Re-run the fast independent module checks:

   ```text
   matlab -batch "addpath(fullfile(pwd,'module_verification')); run_module_verification"
   ```

3. Run one canonical paper case, for example the 28 x 40 reference grid:

   ```text
   matlab -batch "run_scientifically_corrected_case_20260823('Grid_G28x40')"
   ```

4. On Windows, run the complete 16-case paper campaign with:

   ```powershell
   .\run_scientifically_corrected_campaign_20260823.ps1
   ```

   After all 16 `end.mat` files exist, run
   `postprocess_scientifically_corrected_20260823` in MATLAB. The campaign
   launcher resolves `matlab` from `PATH`; it contains no machine-specific
   MATLAB path.

The canonical paper entry points are the dated
`run_scientifically_corrected_*_20260823` wrappers. Calling `main_seal(1)`
directly uses historical demonstration defaults and **does not reproduce the
paper parameter set**. Paper cases obtain their fixed prescribed parameters
from `calibrated_degradation_options.m` through `run_one_corrected_case.m`.

## Runtime and resource requirements

- Tested with MATLAB R2024b. Recreating solver states requires Symbolic Math
  Toolbox for element construction and Statistics and Machine Learning
  Toolbox for percentile summaries. Parallel Computing Toolbox is optional;
  the canonical case wrapper disables an internal pool.
- The PowerShell launchers require Windows PowerShell 5.1 or PowerShell 7 and
  a `matlab` executable available on `PATH`.
- Figure regeneration requires Python 3.10 or later, NumPy, Matplotlib, and an
  installed Times New Roman font. The frozen figures and source-data tables
  can be inspected without Python.
- The isolated module verification completed in about 11 seconds on the
  authors' machine; allow about one minute on a typical workstation. Solver
  runtime is hardware dependent and the full campaign can take hours.
- The campaign launcher starts four MATLAB processes by default. Reduce
  `$maxConcurrent` in the launcher if memory is limited. A workstation with
  at least 16 GB RAM is recommended for the default parallel launch; serial
  execution needs less memory. Allow at least 1 GB of free disk space for
  regenerated states, logs, tables, and figures.

## Meaning of "calibrated" in legacy identifiers

The filename `calibrated_degradation_options.m` and several historical code
comments are retained for interface compatibility. Here, "calibrated" means
the fixed numerical parameterization selected before the reported comparative
campaign. The values are prescribed/tuned model inputs; they were not fitted
to an experimental seal-life, leakage, fatigue, uptake, or aging dataset.
Accordingly, the manuscript describes the results as calibration-dependent
comparisons, not as experimental calibration or validation.

## Corrected implementation

- Fatigue is evaluated as a trial integration-point state during assembly and
  committed once only after the time step converges.
- The phase residual and tangent include the slope of the residual-stiffness
  degradation function.
- A nonconverged Newton or staggered solve is rejected before any path-dependent
  state is committed.
- Each Newton subproblem is capped at 80 corrections, consistent with the manuscript.
- The phase variable is a diffuse degradation field. It is not interpreted as
  a measured crack opening.
- The output `ScreeningIndex` is the maximum of the 99th-percentile phase,
  aging, and fatigue states and is used only for comparisons within this study.
- The energetic fatigue-rate factor is explicit as `energyRateScale = 1.0` (the manuscript's `r_0`).
- With the reported parameters, every activated fatigue point uses the capped rate; energy and pressure determine activation and localization, while the below-threshold branch retains pressure/contact scaling.
- The fluid and aging subproblems use first-order backward Euler.
- The outer staggered comparison is a native-scaled code-level stopping heuristic, not a dimensionless coupled residual norm.
- The Darcy estimate scales with `outOfPlaneWidth`; the maximum-aperture cubic-law upper-bound proxy does not receive a second width multiplier.

## Compatibility identifiers and manuscript notation

Several public property and result-field names predate the final manuscript
notation. They are retained only so existing scripts and saved interfaces
remain compatible; they do not introduce alternative physical definitions.

| Legacy code identifier | Manuscript term | Meaning |
| --- | --- | --- |
| `DL_crack` | `D_L^deg` | degradation-enhanced diffusivity endpoint; the diffuse field is not treated as a discrete crack |
| `D_age` | `D_a` | effective aging-extent transport coefficient, not a separately identified molecular-oxygen diffusivity |
| `Q_poiseuille` | `Q_P` | maximum-aperture cubic-law upper-bound proxy |
| `leakiness` / `leakinessField` | `Lambda_flow` | dimensionless flow-conductance factor |
| `Q` / `leakage_Q` | `Q_D` | Darcy flow estimate |
| `energyRateScale` | `r_0` | energetic fatigue-rate scale |

## Main entry points

- `tests/test_fatigue_single_commit.m`: regression test for repeated assembly
  and single commit.
- `tests/test_leakage_width_semantics.m`: regression test showing that doubling the extrusion width doubles `Q_D` and leaves `Q_P` unchanged.
- `run_scientifically_corrected_campaign_20260823.ps1`: parallel 16-case
  Windows launcher.
- `run_scientifically_corrected_case_20260823.m`: one named study case.
- `postprocess_scientifically_corrected_20260823.m`: strict study QA and
  canonical table export.
- `run_scientifically_corrected_validation_20260823.m`: fixed/contact baseline
  and validation workflow.
- `run_parallel_scientifically_corrected_validation_20260823.ps1`: isolated
  validation workers.
- `finalize_scientifically_corrected_validation_20260823.m`: read-only
  aggregation of completed validation states.
- `run_aging_kinetics_verification.m`: backward-Euler comparison with the
  uniform analytical aging reaction.
- `module_verification/run_module_verification.m`: independent Fick-diffusion
  and fatigue-point checks with predeclared acceptance criteria.
- `PostProcessing/quantify_m7_endpoint_strain.m`: recovers endpoint
  integration-point strain diagnostics from saved states without rerunning a
  case.
- `figures/source/build_figures_v4.py`: nine-figure publication builder.

## Data layout

- `processed_data/parametric_summary.csv`: canonical 16-case table.
- `processed_data/qa/leakage_width_regression_20260824.csv`: compact human-readable width-semantics regression results.
- `processed_data/l9_main_effect_ranges.csv`: descriptive L9 factor-level
  means and ranges.
- `processed_data/qa/`: study integrity and numerical-sensitivity checks.
- `processed_data/aging_kinetics_verification.csv`: analytical reaction-step
  comparison.
- `processed_data/validation/`: three-day grid, time-step, and fixed/contact
  summaries.
- `Results_Paper_ScientificallyCorrected_20260823/M7_Strain_Quantification_20260828/`:
  M7 endpoint-strain tables, an auditable integration-point export, MATLAB
  tables, and bounded suggested text.
- `figures/source_data/`: portable numeric inputs decoded from the corrected
  MATLAB result objects.
- `figures/manuscript/`: 45 final outputs (nine figures in manuscript PNG, 600-dpi PNG, 600-dpi TIFF, SVG, and PDF formats) plus their figure manifest.
- The canonical campaign manifest contains 16 completed records: three mesh, four regularization-length, and nine L9 cases.

Path fields in the processed QA CSVs are portable provenance labels for source solver
states used during export; the corresponding raw states are not distributed.

Raw serialized solver objects are omitted because they are MATLAB-version and
class-definition dependent. The entry points regenerate them. The M7 audit
script therefore requires regenerated `end.mat` states; its already generated
tables are included for inspection. The study was run with MATLAB R2024b. No
experimental data are included or claimed.
