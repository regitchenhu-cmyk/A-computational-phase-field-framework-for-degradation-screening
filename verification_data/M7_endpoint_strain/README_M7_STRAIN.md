# M7 endpoint-strain quantification

Generated from saved `end.mat` states; no model case was rerun.

## Metric

Infinitesimal von Mises-equivalent strain at all interior Gauss points: `sqrt((2/3) dev(eps):dev(eps))`, with `eps = sym(grad(u))` and plane-strain `eps_zz = 0`.

P99 is the whole-domain percentile. The raw maximum and its coordinates are retained to expose constrained-corner localization. The corner window is `x <= 0.20 mm, y <= 0.20 mm`.

Green-Lagrange comparisons are kinematic diagnostics only; they are not a finite-strain validation of the coupled solution.

`det(F) <= 0` marks a kinematic inversion at an integration point. Such raw maxima are reported for auditability but must not be interpreted as physical material strains.

## Reference result

For `Grid_G28x40`, whole-domain P99 = 21.0959% and the raw maximum = 148.746% at (0.0140877, 0.0140877) mm.

## Files

- `m7_strain_summary.csv`: all endpoint cases.
- `m7_mesh_strain_summary.csv`: three main-study meshes.
- `m7_reference_Grid_G28x40_ip_strain.csv`: auditable integration-point values for the reference grid.
- `m7_strain_quantification.mat`: MATLAB tables.

Three mesh rows written: 3.
