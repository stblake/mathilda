# ParametricPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ParametricPlot[{fx, fy}, {t, tmin, tmax}, opts...]`**

Adaptively samples the parametric curve (fx(t), fy(t)) over \[tmin, tmax\] and returns a Graphics\[...\] object (auto-displayed). The body may be any expression that evaluates to a 2-element {x,y} list (not just a literal {fx,fy}). Multiple curves: ParametricPlot\[{{fx1,fy1}, ...}, {t,...}\]. Two-iterator (filled region) form: ParametricPlot\[body, {t,...}, {r,...}\] samples a PlotPoints x PlotPoints grid and emits Polygon\[\] quads. Default AspectRatio -\> 1 (both axes equally important). Options: PlotPoints (default 25), MaxRecursion (default 6), MaxPlotPoints, Mesh (All: dots for curves, grid lines for regions), PlotLegends (Automatic/"Expressions"/{labels...}: draws a legend), ColorFunction ("Rainbow" or f\[t\] / f\[t,r\]), ColorFunctionScaling (default True), RegionFunction (f\[x,y\] mask), PlotStyle, AspectRatio, Axes, PlotRange, PlotRangePadding, AxesLabel, AxesOrigin, Frame, FrameLabel, GridLines, Prolog, Epilog, PlotLabel, Background, ImageSize (all passed through to Graphics).

## Examples (8)

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

### Options (2)

```mathematica
In[7]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}, ColorFunction -> (Hue[#] &)]
Out[7]= -Graphics-

In[8]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}, RegionFunction -> Function[{x, y}, x > 0]]
Out[8]= -Graphics-
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

**Attributes:** `HoldAll`, `Protected`.

## See also

[HoldAll](../../expression-information/HoldAll/)

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_parametricplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parametricplot.c)
