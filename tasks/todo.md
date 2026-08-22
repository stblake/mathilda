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

## Follow-up (closed-orbit overlap)
- [x] Bug: closed orbits drawn TWICE — grow_streamline integrated both
      directions, and each direction walks the whole loop. Fixed: integrate_dir
      reports closure; closed orbit = single forward pass, closed by repeating
      the seed; backward pass skipped. Rings ~halve in point count, first==last.
- [x] Drop sub-resolution closed loops (bbox diag < d_sep) → clean vortex centre
      (removes the tiny "hook" at origin of {-y,x}).
- [x] streamplot/autocompile tests pass; valgrind unchanged (no new leak).

## Density match to Mathematica (mma_example.png)
- [x] "Almost start to overlap" / "compare to Mathematica" = DENSER rings.
      MMA draws ~13 concentric rings nearly touching; ours drew 7 sparse.
- [x] Default StreamPoints 15 -> 25 (d_sep = min(extent)/25) -> ~12 rings that
      nearly touch, matching MMA. Verified rotation + saddle + complex render
      clean (no merge/overlap) at the higher density.
- [x] Book figure regenerated; streamplot_tests pass.
## Dashed-arrow style (exact MMA match)
- [x] emit_dashed_stream: walk arc length, alternate dash (curve-following
      make_line_range + chevron_tip Arrow at the tip) / gap. Dash 1.6*d_sep,
      gap 1.0*d_sep, chevron 0.7*d_sep. Renderer sizes Arrow head from last
      segment -> dash must be Line + separate chevron, not one multi-pt Arrow.
- [x] Animate path unchanged (solid AnimatedStreamline + periodic chevrons).
- [x] rotation/saddle/complex render in MMA dashed style; tests pass; valgrind
      unchanged (no new leak). Book figure regenerated. Doc + changelog updated.

## Fix: chunky arrowheads (dashed style looked worse)
- [x] Root cause: PDF exporter (graphics_export.c, SEPARATE from render.c) used
      a FIXED 8pt arrowhead -> swamped the ~230 short chevrons of a dashed plot.
- [x] graphics_export Arrow head now scales: hl=min(0.55*final_seg, 8), hw=hl*.38.
      Long arrows (VectorPlot) unchanged at cap; short ones shrink.
- [x] Streamplot: lighter lines (Thickness 0.006->0.0035), tuned dash_arrow 0.8.
- [x] rotation/saddle look like MMA now; VectorPlot fine; tests + valgrind OK.
      Book figures (streamplot + vectorplot) regenerated.

## Fix: arrow-with-a-straight-line (user: "looks ridiculous")
- [x] Root cause: chevron_tip emitted Arrow[{a,b}]; BOTH renderers draw every
      Arrow as shaft-polyline + head, so each dash tip grew a ~0.8*d_sep straight
      stick (a chord across the curved streamline) plus a head.
- [x] chevron_tip now emits a bare filled arrowhead Polygon[{b1,b2,tip}] — no
      shaft. Tip at dash end, base dash_arrow upstream, direction = chord over
      the head's own arc length (stays aligned on tight bends). Fills with the
      stream's colour directive on both paths. Animate path unchanged.
- [x] User follow-up: arrowheads 0.8x (dash_arrow 0.42->0.336), streamlines
      1.25x thicker (Thickness 0.0035->0.004375).
- [x] rotation + saddle render clean; streamplot_tests pass; build clean.

## Considered / deferred
- Perpendicular offspring seeding (full JL) would fill sparse-region gaps even
  better; the finer candidate grid + culling already gives good coverage, so
  left out to keep the change contained.
- PDF vector export still stretches square domains to the page (pre-existing,
  affects all plotters); windowed/PNG render honors AspectRatio. Out of scope.
