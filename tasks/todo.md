# StreamPlot: evenly-spaced streamlines

## Problem
`StreamPlot[{-y, x}, {x, -3, 3}, {y, -3, 3}]` (rotation → concentric circles)
rendered as a hairball of ~225 short, overlapping, fragmented arcs. Careless
with stream length and curvature.

## Root cause (src/graphics/streamplot.c)
1. One forward-only stream per grid seed (15×15 = 225) → overlap.
2. StreamScale default 0.08 → every stream a short stub.
3. Raw-field RK4 (step h·|v|) → curvature resolution tied to local speed.

## Fix — Jobard–Lefebvre evenly-spaced streamlines
- [x] Normalized-field RK4 (fixed arc-length step) → uniform spacing/curvature.
- [x] Grow each line both directions from its seed.
- [x] Terminate on: proximity to another line (½·d_sep, hash-grid O(1)),
      boundary, critical point, or closed-orbit return (draws whole circles).
- [x] Even placement: candidate grid finer than d_sep + cull seeds within
      d_sep of an existing line. StreamPoints sets d_sep.
- [x] Periodic direction chevrons along each Line (not one lone mid-arrow).
- [x] StreamScale default → run to natural end (None/Automatic too); s>0 caps.
- [x] Doc (graphics.md) + changelog updated.

## Verification
- Rendered rotation / saddle {x,-y} / {-y Exp[-x^2], x Sin[y]} → all clean,
  evenly-spaced, correct flow direction and speed coloring.
- Point counts: rotation 225 fragments → 19 clean circles.
- streamplot_tests + autocompile_tests pass.
- valgrind: leak total identical to baseline and to the OLD code (13,496 B /
  421 blocks), constant across 1 vs 20 lines → no new leak (pre-existing fixed
  interning + macOS baseline noise).
- Fast: ~0.03s for the field evals.

## Considered / deferred
- Perpendicular offspring seeding (full JL) would fill sparse-region gaps even
  better; the finer candidate grid + culling already gives good coverage, so
  left out to keep the change contained.
- PDF vector export still stretches square domains to the page (pre-existing,
  affects all plotters); windowed/PNG render honors AspectRatio. Out of scope.
