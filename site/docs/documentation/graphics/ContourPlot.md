# ContourPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Generates iso-contour lines of f(x,y) using the marching squares algorithm and returns a Graphics\[...\] object (auto-displayed). ContourPlot is HoldAll: f is held unevaluated until x and y are bound to numeric values. Options: Contours         - Integer n (n evenly spaced auto levels, default 10), or {c1, c2, ...} (explicit contour values). ContourStyle     - Style directive(s) for the contour lines. A single directive is applied to all levels; a List cycles through the levels. Automatic (default) colours by height. None/False suppresses lines (leaves only shading). ContourLabels    - True: draw the z value at the midpoint of each level's first visible segment. Default False. ContourShading   - True: fill each grid cell by its z value (via ColorFunction or the built-in thermal gradient). False/None: lines only. Automatic (default): shade when ColorFunction is set, otherwise lines only. ColorFunction    - A function f\[t\] → color (t in \[0,1\] after scaling), or a named ramp string: "Rainbow", "Temperature", "CoolTones", "WarmTones", "Greyscale". Applied to shading and auto line colors. ColorFunctionScaling - True (default): normalise z to \[0,1\] before calling ColorFunction. False: pass raw z. PlotPoints       - Grid resolution per axis (default 25; increase for smoother contours). RegionFunction   - f\[x,y\] mask: cells where the function is False are skipped (neither shaded nor contoured). Standard Graphics options (Axes, AspectRatio, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog, Epilog, ...) pass through to the Graphics\[...\] result.

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[1]= -Graphics-
```

### Options (9)

```mathematica
In[2]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, Contours -> 5]
Out[2]= -Graphics-

In[3]:= ContourPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, Contours -> {-2, -1, 0, 1, 2}]
Out[3]= -Graphics-

In[4]:= ContourPlot[Sin[x + y], {x, -3, 3}, {y, -3, 3}, ColorFunction -> "Rainbow", ContourShading -> True]
Out[4]= -Graphics-

In[5]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, ContourShading -> True, Contours -> 8]
Out[5]= -Graphics-

In[6]:= ContourPlot[Sin[x] Cos[y], {x, -Pi, Pi}, {y, -Pi, Pi}, ContourShading -> False, ContourStyle -> {Thickness[0.006]}, PlotPoints -> 40]
Out[6]= -Graphics-

In[7]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, ContourLabels -> True, Contours -> 5]
Out[7]= -Graphics-

In[8]:= ContourPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 4], ContourShading -> True]
Out[8]= -Graphics-

In[9]:= ContourPlot[Sin[x + y], {x, -3, 3}, {y, -3, 3}, ContourStyle -> {Red, Blue, Green}, Contours -> 6]
Out[9]= -Graphics-

In[10]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3}, ContourStyle -> None, ContourShading -> True, ColorFunction -> "Temperature"]
Out[10]= -Graphics-
```

### Applications (5)

```mathematica
In[11]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[11]= -Graphics-

In[12]:= ContourPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, Contours -> {-2, -1, 0, 1, 2}]
Out[12]= -Graphics-

In[13]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, ContourShading -> True, Contours -> 8]
Out[13]= -Graphics-

In[14]:= Head[ContourPlot[x^2 + y^2 == 1, {x, -2, 2}, {y, -2, 2}]]
Out[14]= Graphics

In[15]:= Length[Cases[ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, Contours -> 5, PlotPoints -> 25], _Line, Infinity]] > 0
Out[15]= True
```

## Algorithm

contourplot.c — ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]

Generates iso-contour lines of a 2D function f(x, y) using the marching squares algorithm and returns a Graphics[...] object auto-displayed by the

```text
REPL.  ContourPlot is HoldAll: f and the iterator specs are held unevaluated
```

until x and y are bound to numeric values.

Algorithm:

```text
  1. Evaluate f on a (PlotPoints+1) × (PlotPoints+1) grid.
  2. Choose contour levels: explicit list, or N evenly spaced levels.
  3. Optionally shade each grid cell by its z value (ContourShading).
  4. For each level: run marching squares over the grid, emitting Line[]
     segments.  Saddle cells (states 5 and 10) use the bilinear centre
     value to pick the correct one of the two possible pairings.
  5. Optionally label each level at the midpoint of its first segment.
  6. Wrap everything in Graphics[...] with Plot's default options. 
```

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- W. E. Lorensen and H. E. Cline, *Marching Cubes: A High Resolution 3D Surface Construction Algorithm*, SIGGRAPH Comput. Graph. **21** (1987) 163-169 — marching squares is its 2D case.
- Source: [`src/graphics/contourplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/contourplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

### Notes

`ContourPlot` is `HoldAll`. It samples `f` on an `(N+1) x (N+1)` grid (`PlotPoints`
default 75) and runs **marching squares** per contour level, disambiguating the
two saddle cell configurations by the bilinear cell-centre value. Each isoline
segment is emitted as a `Line[...]`. Levels come from `Contours -> n` (a count),
`Contours -> {...}` (explicit values), or the default of 10 interior levels evenly
spaced between the data min and max.

By default (Automatic) only the lines are drawn; shading turns on when a
`ColorFunction` is given, or explicitly with `ContourShading -> True`, filling
each cell with a `Rectangle` coloured by its normalised average.
`ContourPlot[lhs == rhs, ...]` rewrites the body to `lhs - rhs` and draws the
single zero-contour; a list of equations plots each as its own curve. There is no
adaptive refinement, so contour smoothness is set by `PlotPoints`.
