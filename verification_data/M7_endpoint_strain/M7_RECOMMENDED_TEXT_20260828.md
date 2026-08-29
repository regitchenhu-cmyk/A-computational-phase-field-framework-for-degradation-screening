# M7 recommended manuscript and response text

## Definition to state once

The endpoint infinitesimal strain was recovered at every interior Gauss
point from the saved displacement field. The reported scalar is the
three-dimensional plane-strain deviatoric equivalent strain

\[
\varepsilon_{\mathrm{eq}}
=\sqrt{\frac{2}{3}\,\mathrm{dev}(\boldsymbol\varepsilon):
\mathrm{dev}(\boldsymbol\varepsilon)},\qquad
\boldsymbol\varepsilon=\mathrm{sym}\nabla\mathbf u,
\quad \varepsilon_{zz}=0 .
\]

P99 is the unweighted percentile over all interior integration points.
The raw maximum is retained separately and is not treated as a converged
material strain.

## Suggested manuscript paragraph

At the saved 20.063-day endpoint of the 28 x 40 reference case, the
whole-domain unweighted 99th percentile of the integration-point
deviatoric equivalent infinitesimal strain was 21.10%, whereas the raw
integration-point maximum was 148.75% at the lower-left pressure/constraint
corner. Across the 18 x 26, 28 x 40, and 40 x 56 meshes, P99 decreased from
25.74% to 21.10% and 18.06%, while the raw maximum increased from 101.98%
to 148.75% and 200.75%, respectively. The raw maximum is therefore a
corner-localized, mesh-sensitive diagnostic rather than a converged
material strain. Across the nine L9 endpoints, P99 ranged from 21.11% to
39.14%; the upper part of this range confirms that the small-strain model
is a qualitative screening idealization under the most severe prescribed
conditions. These endpoint statistics bound the kinematic limitation but
do not constitute a finite-strain error estimate.

## Suggested response to Major Comment 7

We agree that the small-strain assumption required quantitative context.
Without rerunning the cases, we recovered the displacement gradient from
the saved endpoint states at every interior Gauss point and added a
post-hoc strain audit. We report the whole-domain unweighted P99 and the
raw point maximum separately. For the 28 x 40 reference case, P99 was
21.10%, while the raw maximum was 148.75% at the lower-left constrained
pressure corner. Mesh refinement reduced P99 (25.74%, 21.10%, and 18.06%)
but increased the raw maximum (101.98%, 148.75%, and 200.75%), identifying
the latter as a non-converged corner hotspot. The nine L9 endpoint P99
values ranged from 21.11% to 39.14%. We have therefore restricted the
interpretation to within-study numerical screening and now state explicitly
that the severe-condition results are not quantitatively reliable finite-
deformation predictions. The audit bounds the limitation; it does not
claim a finite-strain accuracy percentage.

## Internal caution before insertion

A diagnostic reconstruction of `F = I + grad(u)` found two integration
points with `det(F) <= 0` in each of O05 and O07--O09, all at the upper-right
corner. The three-day contact-aware verification likewise has two such
points. These are extremely sparse corner outliers, but they confirm that
the raw maxima cannot be interpreted as physical finite-deformation
strains. Do not replace this caveat with a claim such as "Z% accuracy": the
available data do not contain a validated finite-strain constitutive or
solution-error benchmark.
