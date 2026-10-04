# ParametricPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ParametricPlot[{fx, fy}, {t, tmin, tmax}, opts...]`**

Adaptively samples the parametric curve (fx(t), fy(t)) over \[tmin, tmax\] and returns a Graphics\[...\] object (auto-displayed). The body may be any expression that evaluates to a 2-element {x,y} list (not just a literal {fx,fy}). Multiple curves: ParametricPlot\[{{fx1,fy1}, ...}, {t,...}\]. Two-iterator (filled region) form: ParametricPlot\[body, {t,...}, {r,...}\] samples a PlotPoints x PlotPoints grid and emits Polygon\[\] quads. Default AspectRatio -\> 1 (both axes equally important). Options: PlotPoints (default 25), MaxRecursion (default 6), MaxPlotPoints, Mesh (All: dots for curves, grid lines for regions), PlotLegends (Automatic/"Expressions"/{labels...}: draws a legend), ColorFunction ("Rainbow" or f\[t\] / f\[t,r\]), ColorFunctionScaling (default True), RegionFunction (f\[x,y\] mask), PlotStyle, AspectRatio, Axes, PlotRange, PlotRangePadding, AxesLabel, AxesOrigin, Frame, FrameLabel, GridLines, Prolog, Epilog, PlotLabel, Background, ImageSize (all passed through to Graphics).

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[1]= -Graphics-

In[2]:= ParametricPlot[{Sin[2 t], Sin[3 t]}, {t, 0, 2 Pi}]
Out[2]= -Graphics-

In[3]:= ParametricPlot[{{Cos[t], Sin[t]}, {2 Cos[t], Sin[t]}}, {t, 0, 2 Pi}]
Out[3]= -Graphics-

In[4]:= ParametricPlot[2 {Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[4]= -Graphics-

In[5]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}]
Out[5]= -Graphics-

In[6]:= ParametricPlot[r^2 {Sqrt[t] Cos[t], Sin[t]}, {t, 0, 3 Pi/2}, {r, 1, 2}]
Out[6]= -Graphics-
```

### Options (3)

```mathematica
In[7]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}, ColorFunction -> (Hue[#] &)]
Out[7]= -Graphics-

In[8]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}, RegionFunction -> Function[{x, y}, x > 0]]
Out[8]= -Graphics-

In[9]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}, Mesh -> All]
Out[9]= -Graphics-
```

### Applications (5)

```mathematica
In[10]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[10]= -Graphics-

In[11]:= ParametricPlot[{Sin[2 t], Sin[3 t]}, {t, 0, 2 Pi}]
Out[11]= -Graphics-

In[12]:= ParametricPlot[2 {Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[12]= -Graphics-

In[13]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}]
Out[13]= -Graphics-

In[14]:= Attributes[ParametricPlot]
Out[14]= {HoldAll, Protected}
```

## Algorithm

parametricplot.c — ParametricPlot[body, {t, tmin, tmax}, opts...]

```text
                 — ParametricPlot[body, {t, tmin, tmax}, {r, rmin, rmax}, opts...]
```

HoldAll: the body and all iterator specs are unevaluated when received.

One-iterator form:

```text
  body must evaluate to {x, y} for each t. Single body or
  {body1, body2,...} (list of bodies) for multi-curve.
  Adaptive 2D sampling with a Euclidean chord-deviation flatness test.
  Output: Graphics[{Line[...], ...}, opts].
```

Two-iterator form:

```text
  body evaluates to {x, y} for each (t, r) pair.
  Samples a PlotPoints x PlotPoints grid and builds Polygon[] quads,
  like Plot3D but mapped back into the xy-plane.
  Output: Graphics[{Polygon[...], ...}, opts].
```

In both forms the body can be any expression that evaluates to a 2-element numeric list: a literal {fx, fy}, or a computed form such as r^2 * {Sqrt[t] Cos[t], Sin[t]}, as long as the result has head List and exactly two finite-real elements.

## Implementation notes

**Algorithm.** `builtin_parametricplot` is `HoldAll`. `is_iterator` tests each
trailing `{var, min, max}` spec; one iterator is the **curve** form, two is the
**filled region** form. The body is any expression that must evaluate to a
2-element list of finite reals — a literal `{fx, fy}` or a computed form such as
`r {Cos[t], Sin[t]}`. Two evaluation paths live in `ParamEvalCtx`: a compiled
fast path (`param_ctx_compile`, all-or-nothing, firing only for a *literal*
2-element `List` whose two coordinates each compile with `autocompile_new`) and
an interpreter fallback (`param_eval_at` → `eval_body_xy`, which binds the
iterator via an OwnValue and checks the result is a finite real pair). The
one-iterator form (`build_param_curve` → `param_sample`) runs its **own**
adaptive bisection (`param_subdivide`), the same three-probe refinement as `Plot`
but measuring **Euclidean chord deviation normalised by the bounding-box
diagonal** (neither axis is privileged for a parametric curve), deliberately kept
at the same tolerance `PARAM_FLAT_TOL = 0.0006`; invalid or region-rejected
samples break the polyline into separate `Line[...]` runs (one coloured two-point
`Line` per segment under `ColorFunction`). The two-iterator form
(`build_param_region`) samples a **uniform `n x n`** `(t, r)` grid, maps each pair
to `(x, y)`, and emits one filled `Polygon[{p00, p10, p11, p01}]` per cell whose
four corners are all valid — no adaptivity — with `Mesh -> All` overlaying the
grid edges as `Line`s. A `List`-of-`List`s body is the multi-curve form, each
sub-body drawn in a `palette_color` with non-colour directives given their own
`List` scope. The result is an inert `Graphics[prims, opts...]`.

**Data structures.** `ParamEvalCtx` (iterator vars, body, `acx`/`acy` compiled
coordinates, `RegionFunction`); `ParamPt`/`ParamBuf` growing sample buffers; the
option bundle from `split_options_param`.

**Complexity / limits.** `PlotPoints` initial samples (default **50** for the
curve, **75** for the region) plus refinement to `MaxRecursion = 6` for the
curve; the region grid is uniform with no refinement and drops any cell with an
invalid corner (no partial-triangle handling). Default `AspectRatio -> 1`,
`Axes -> True`, `PlotStyle` colour `RGBColor[0.2, 0.4, 0.8]`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/parametricplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/parametricplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_parametricplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parametricplot.c)

## Notes & additional examples

### Notes

The one-iterator form draws a curve `{fx(t), fy(t)}`; the two-iterator form fills
a region `{fx(t, r), fy(t, r)}` with `Polygon` quads over a uniform grid. The body
is any expression that evaluates to a 2-element list of finite reals — a literal
`{fx, fy}` takes a compiled fast path, a computed form such as `2 {Cos[t],
Sin[t]}` takes the interpreter path.

The curve form runs its **own** adaptive sampler (parallel to `Plot`'s), refining
by Euclidean deviation from the chord normalised by the bounding-box diagonal,
since neither axis is privileged for a parametric curve; invalid samples break the
curve into separate `Line[...]` runs. The region form is a uniform grid with no
refinement (`Mesh -> All` overlays its grid lines). Default `AspectRatio -> 1`.
