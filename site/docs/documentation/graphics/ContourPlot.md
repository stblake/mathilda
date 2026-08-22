# ContourPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Generates iso-contour lines of f(x,y) using the marching squares algorithm and returns a Graphics\[...\] object (auto-displayed). ContourPlot is HoldAll: f is held unevaluated until x and y are bound to numeric values. Options: Contours         - Integer n (n evenly spaced auto levels, default 10), or {c1, c2, ...} (explicit contour values). ContourStyle     - Style directive(s) for the contour lines. A single directive is applied to all levels; a List cycles through the levels. Automatic (default) colours by height. None/False suppresses lines (leaves only shading). ContourLabels    - True: draw the z value at the midpoint of each level's first visible segment. Default False. ContourShading   - True: fill each grid cell by its z value (via ColorFunction or the built-in thermal gradient). False/None: lines only. Automatic (default): shade when ColorFunction is set, otherwise lines only. ColorFunction    - A function f\[t\] → color (t in \[0,1\] after scaling), or a named ramp string: "Rainbow", "Temperature", "CoolTones", "WarmTones", "Greyscale". Applied to shading and auto line colors. ColorFunctionScaling - True (default): normalise z to \[0,1\] before calling ColorFunction. False: pass raw z. PlotPoints       - Grid resolution per axis (default 25; increase for smoother contours). RegionFunction   - f\[x,y\] mask: cells where the function is False are skipped (neither shaded nor contoured). Standard Graphics options (Axes, AspectRatio, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog, Epilog, ...) pass through to the Graphics\[...\] result.

## Examples (10)

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
