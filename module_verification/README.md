# Independent module verification for Major Comment 3

This directory contains two independent analytical checks. They reuse the
reported equations and parameter values but do not call or modify the
production solver.

## Checks

### 1. One-dimensional Fick diffusion

The check solves

`dC/dt - D d2C/dx2 = 0`, `0 < x < L`,

with the mixed boundary conditions `C(0,t)=0.65` and `dC/dx(L,t)=0`, and the
initial condition `C(x,0)=0`. The numerical calculation uses consistent-mass
linear finite elements and backward Euler, matching the capacity and diffusion
terms assembled in `Models/@FluidDiffusion/FluidDiffusion.m`.

The reference solution is the Fourier series

`C/Cs = 1 - sum[2/mu_n sin(mu_n x/L) exp(-D mu_n^2 t/L^2)]`,

where `mu_n=(n+1/2)pi`. Its spatial mean is evaluated analytically from the
same series. Six calculations combine two grids (`Nx=14,70`) and three fixed
time steps (`12 h, 3 h, 0.25 h`). The common `0.25 h` step is used to isolate
the spatial-refinement comparison; the larger steps retain a separate temporal
refinement check. The final concentration distribution, mean history,
absolute/relative L2 error, Linf error, and mean error are reported.

Parameters: `L=3.5 mm`, `D=6.0e-13 m^2/s`, `Cs=0.65`, and `t_end=20 days`.

### 2. Constant-Wplus fatigue point update

The check reproduces the single-integration-point update in
`Models/@FatigueDamage/FatigueDamage.m`:

`D_(n+1) = min(D_n + r DeltaN, 1)`.

For constant tensile energy, pressure amplitude, and contact multiplier, the
closed form is `D(N)=min(D0+rN,1)`. Three cases exercise:

- the uncapped below-threshold base-rate branch;
- an activated but uncapped low-pressure-amplitude branch;
- the capped active branch at the main pressure amplitude.

Cycle jumps of 1, 37, and 1000 cycles are compared with the closed form.
The parameter values come from `calibrated_degradation_options.m` and the
FatigueDamage inputs in `main_seal.m`.

## Predeclared acceptance criteria

- Fine (`Nx=70`, `dt=0.25 h`) Fick relative L2 error: `<= 1.0e-2`.
- Fine Fick relative Linf error: `<= 2.0e-2`.
- Fine Fick mean-concentration relative error: `<= 1.0e-2`.
- Fick spatial-refinement L2 error ratio at the common `0.25 h` step: `< 1`.
- Fick temporal-refinement L2 error ratio on `Nx=70`: `< 1`.
- Both capped and uncapped fatigue branches, including activated examples,
  must be present.
- Maximum fatigue discrete-versus-closed-form history error: `<= 5.0e-12`.
  This is an absolute double-precision roundoff tolerance for as many as
  120,000 repeated scalar additions; it is not a model-error allowance.
- Final-fatigue-state range across cycle jumps: `<= 1.0e-12`.

## Run

Open MATLAB or a terminal in the extracted supplementary-archive root. No
drive-letter mapping or machine-specific path is required.

```matlab
addpath(fullfile(pwd, 'module_verification'))
run_module_verification
```

For a command-line batch run:

```text
matlab -batch "addpath(fullfile(pwd,'module_verification')); run_module_verification"
```

The check was run with MATLAB R2024b and uses only small dense/sparse
one-dimensional systems. It completed in about 11 seconds on the authors'
machine; allow one minute and less than 1 GB RAM on a typical workstation.
It neither launches the coupled production solver nor requires its saved
states.

## Generated outputs

The script creates `outputs/` and writes:

- `verification_status.csv`
- `parameter_provenance.csv`
- `fick_summary.csv`
- `fick_profiles.csv`
- `fick_average_history.csv`
- `fick_verification_results.mat`
- `fick_profiles.png`
- `fick_average_history.png`
- `fatigue_point_summary.csv`
- `fatigue_point_rates.csv`
- `fatigue_point_histories.csv`
- `fatigue_verification_results.mat`
- `fatigue_point_verification.png`
- `module_verification_all_results.mat`
- `module_verification.log`

Precomputed copies of all outputs are included. The three MAT files have
different purposes:

- `fick_verification_results.mat` contains the Fick parameters, all six case
  profiles and mean histories, summary table, and acceptance status.
- `fatigue_verification_results.mat` contains the point-test parameters,
  branch definitions, rates, histories, summary table, and acceptance status.
- `module_verification_all_results.mat` combines both result structures with
  the overall status table and scalar `allPassed` flag.

These calculations verify the isolated numerical formulas and their
implementations. They are not material calibration, experimental validation,
or validation of seal service life.
