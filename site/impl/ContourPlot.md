---
references:
  - "W. E. Lorensen and H. E. Cline, *Marching Cubes: A High Resolution 3D Surface Construction Algorithm*, SIGGRAPH Comput. Graph. **21** (1987) 163-169 — marching squares is its 2D case."
source: src/graphics/contourplot.c
---
**Algorithm.** `builtin_contourplot` is `HoldAll`. It evaluates `f` on an
`(N+1) x (N+1)` point grid (`N = PlotPoints`, default **75**, since contours are
plain polylines with no adaptive refinement so the cell pitch is the facet size),
compiling the body once with `grid_compile` → `autocompile_new(body, {x, y}, 2)`
and falling back to the interpreter per sample. For each contour level,
`cell_march` computes a 4-bit corner state and emits the isoline segment(s):
states 0 and 15 give nothing, twelve states give a single segment, and the two
**saddle states (5 and 10)** are disambiguated by the bilinear cell-centre value
`(v00+v10+v11+v01)/4`, which picks the connecting diagonal. `edge_pt` places each
crossing by linear interpolation along the edge, and `push_segment` emits a
`Line[{{x1,y1},{x2,y2}}]`; cells with a non-finite corner are skipped. Levels come
from `Contours -> n` (count), `Contours -> {list}` (explicit), or the default of
**10** interior levels evenly spaced strictly between the data min and max.
`ContourShading` (default Automatic = shading only when a `ColorFunction` is
given) fills each cell with a `Rectangle` coloured by its clamped normalised
corner average; `ContourStyle` cycles an explicit colour list across levels (or
`None` suppresses lines); `ContourLabels -> True` places a `Text[level]` at the
midpoint of each level's first segment. `ContourPlot[lhs == rhs, ...]` rewrites
the body to `lhs - rhs` with a single level at 0; a `List` of equations plots each
as its own zero-contour in a cycling palette. The result is an inert
`Graphics[prims, opts...]`.

**Data structures.** `ContourOpts`; a flat `double* grid` of `(N+1)^2` samples;
`GridCtx {xvar, yvar, body, AutoCompiled* ac}`; the cell state is a transient
per-cell `int` (no persistent cell buffer); `Expr** prims` grown by `ENSURE_CAP`.
The default colour ramp is `default_ramp_rgb` (Viridis); the ColorFunction
callback is tried at arity 1 (`f[t]`) then 3 (`f[x, y, t]`).

**Complexity / limits.** Grid evaluation `O((N+1)^2)` (~5.6k samples at the
default 75); marching squares `O(n_levels * N^2)`. There is no adaptive
refinement, so contour smoothness is set by `PlotPoints`; defaults `Frame -> True`,
`Axes -> False`, `AspectRatio -> 1`.
