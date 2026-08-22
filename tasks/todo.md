# ParametricPlot / PolarPlot under-sampling fix

## Problem
Default ParametricPlot/PolarPlot curves look angular: the adaptive refiner
barely engages.

## Root cause (measured)
- `PARAM_FLAT_TOL = 0.0025` (fraction of bbox diagonal) in
  `src/graphics/parametricplot.c` is ~4x looser than the `y=f(x)` sampler's
  `FLAT_TOL = 0.0006` (fraction of y-range) in `src/graphics/sampling.c`.
  A default circle refines only to depth 1 -> 49 pts. At 0.0006 -> depth 2 ->
  97 pts (sub-pixel, smooth). Roses/Lissajous scale up proportionally.
- 1-iterator default seed = 25, below Plot's 50 (anti-alias floor the refiner
  cannot recover).
- PolarPlot forces PlotPoints=75 and delegates to ParametricPlot, so it
  inherits the tolerance fix. Its plain circle is already smooth at 75.

## Plan
- [x] Investigate + measure point counts across tolerances (done)
- [x] `PARAM_FLAT_TOL`: 0.0025 -> 0.0006 (parity with sampling.c); update comment
- [x] 1-iterator default seed: 25 -> 50 (match Plot); update comment + call site
- [x] Fix stale "ParametricPlot's 25" comment in polarplot.c
- [x] Changelog note in docs/spec/changelog/2026-08-17.md
- [x] Rebuild; re-measure (circle smooth, curves reasonable)
- [x] Run test_parametricplot + test_autocompile (no regressions)

## Review

Root cause was the flatness threshold, not the seed count: `PARAM_FLAT_TOL`
(fraction of bbox diagonal) was 0.0025, ~4x looser than the y=f(x) sampler's
FLAT_TOL=0.0006. A default circle stopped at recursion depth 1. Fixed both
levers for parity with Plot:

- `PARAM_FLAT_TOL` 0.0025 -> 0.0006 (primary fix)
- 1-iterator default PlotPoints 25 -> 50 (match Plot; anti-alias floor)

Measured point counts (default options), before -> after:
- ParametricPlot circle     49  -> 99
- PolarPlot circle          75  -> 75  (already sub-pixel; unchanged)
- ParametricPlot rose      177  -> 357
- PolarPlot rose Cos[5t]   149  -> 297
- ParametricPlot Lissajous 153  -> 284
- PolarPlot spiral          93  -> 178

Verification:
- test_parametricplot, test_autocompile, graphics_tests: all pass.
- PDF renders (headless vector path) of circle / 5-petal rose / Lissajous:
  all smooth, including high-curvature petal tips and turning points.
- No API/option change; only default sampling density. Tests pin explicit
  PlotPoints/MaxRecursion + relative counts, so none depended on the old
  defaults.

Considered but deferred: adding sampling.c's MAX_CHORD_FRAC chord-length
backstop to the parametric sampler. It addresses a different artifact (long
on-screen gaps in steep-but-locally-straight stretches), not the reported
angularity, which the tolerance fix fully cures. Left out to keep the change
minimal.
