# Graphics

A 2D graphics engine modeled on Mathematica's: symbolic primitives compose
into a `Graphics[...]` object, which is rendered in an interactive window
via a Raylib backend. Rendering is owned by the REPL "front end": **any
top-level result whose head is `Graphics` is auto-displayed**, just as
Mathematica's notebook front end does. `Show[]` and `Plot[]` therefore
simply build and return a `Graphics[...]` object -- the window is opened by
the REPL when that object is the output, so `g // Graphics`, `Show[g]` and
`Plot[...]` all render through one path. A trailing `;` (whose result is
`Null`) suppresses the window. The textual `Out[n]=` line still shows the
`-Graphics-` placeholder, matching Mathematica's behavior outside a
notebook.

Requires Raylib (a system library, autodetected via `pkg-config` at build
time -- `brew install raylib` on macOS, `apt install libraylib-dev` on
Debian/Ubuntu). When Raylib isn't available, the build still succeeds;
`Show`/`Plot` print a one-line message instead of opening a window.

A small 3D extension lives alongside this: `Plot3D[]` builds a
`Graphics3D[...]` object (rendered with an orbit camera) the same way `Plot[]`
builds a `Graphics[...]`, reusing the same primitive heads (`Polygon`, `Line`,
color directives) with 3-coordinate `{x,y,z}` points instead of 2-coordinate
`{x,y}` ones -- see `Graphics3D` and `Plot3D` below.

---

## Feature summary

### 2D plotters

| Function | Purpose | Key options |
|----------|---------|-------------|
| [`Plot`](#plot) | Adaptive function curve(s) of one variable | `PlotPoints`, `MaxRecursion`, `Filling`, `ColorFunction`, `RegionFunction`, `Exclusions`, `ScalingFunctions` |
| [`ListPlot`](#listplot) | Scatter / line plot of discrete data | `Joined`, `Filling`, `DataRange`, `PlotStyle`, `ScalingFunctions` |
| [`ParametricPlot`](#parametricplot) | Parametric curves `{fx(t), fy(t)}` (1-iter) or surfaces `{fx(t,u), fy(t,u)}` (2-iter) | `ColorFunction`, `RegionFunction`, `Mesh` |
| [`PolarPlot`](#polarplot) | Polar curve `r(θ)` | `ColorFunction`, `PlotStyle`, `RegionFunction` |
| [`ContourPlot`](#contourplot) | Iso-contour lines / filled contours of `f(x,y)` | `Contours`, `ContourStyle`, `ContourShading`, `ColorFunction`, `ContourLabels`, `ScalingFunctions` |
| [`DensityPlot`](#densityplot) | Heatmap of `f(x,y)` | `ColorFunction`, `ColorFunctionScaling`, `RegionFunction`, `PlotLegends`, `ScalingFunctions` |
| [`VectorPlot`](#vectorplot) | Arrow grid for a 2D vector field | `VectorPoints`, `VectorScale`, `VectorStyle`, `ColorFunction`, `RegionFunction`, `ScalingFunctions` |
| [`StreamPlot`](#streamplot) | RK4-integrated streamlines of a 2D vector field | `StreamPoints`, `StreamScale`, `StreamStyle`, `ColorFunction`, `RegionFunction`, `ScalingFunctions` |
| [`BarChart`](#barchart) | Vertical bar chart for categorical / grouped data | `ChartStyle`, `ChartLabels`, `BarSpacing` |
| [`Histogram`](#histogram) | Frequency histogram with flexible bin specs | `ChartStyle`, `BarSpacing`; bin: `k`, `{step}`, `{min,max,step}` |

### 3D plotters

| Function | Purpose | Key options |
|----------|---------|-------------|
| [`Plot3D`](#plot3d) | Surface plot of `f(x,y)` with orbit camera | `PlotPoints`, `Mesh`, `ColorFunction`, `RegionFunction`, `Lighting` |
| [`ParametricPlot3D`](#parametricplot3d) | Parametric space curves or surfaces in 3D | `ColorFunction`, `RegionFunction`, `Mesh` |

### Shared option cross-reference

| Option | Supported by |
|--------|-------------|
| `PlotRange` | all plotters (auto-embedded; user can override) |
| `PlotPoints` | all function plotters |
| `PlotStyle` | `Plot`, `ListPlot`, `ParametricPlot`, `PolarPlot`, `Plot3D`, `ParametricPlot3D` |
| `PlotLegends` | `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `StreamPlot` |
| `PlotLabel` | all (pass-through to `Graphics`/`Graphics3D`) |
| `ColorFunction` | `Plot`, `ParametricPlot`, `PolarPlot`, `Plot3D`, `ParametricPlot3D`, `DensityPlot`, `ContourPlot`, `VectorPlot`, `StreamPlot` |
| `ColorFunctionScaling` | same set as `ColorFunction` |
| `RegionFunction` | `Plot`, `ParametricPlot`, `Plot3D`, `ParametricPlot3D`, `DensityPlot`, `ContourPlot`, `VectorPlot`, `StreamPlot` |
| `Filling` / `FillingStyle` | `Plot`, `ListPlot` |
| `ScalingFunctions` | `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot`, `StreamPlot` |
| `AspectRatio` | all 2D plotters (default `1/GoldenRatio` for Plot/ListPlot, `1` for field plots) |
| `Axes` / `AxesLabel` / `AxesStyle` | all plotters (pass-through) |
| `Frame` / `FrameLabel` / `FrameStyle` | all plotters (pass-through) |
| `GridLines` / `GridLinesStyle` | all plotters (pass-through) |
| `Background` / `ImageSize` | all plotters (pass-through) |

### Primitives

| Primitive | Description |
|-----------|-------------|
| `Point[{x,y}]` / `Point[{{x1,y1},…}]` | Point(s) |
| `Line[{{x1,y1},…}]` | Polyline |
| `Arrow[{{x1,y1},{x2,y2}}]` | Directed line segment with arrowhead |
| `Circle[{x,y}, r]` | Circle outline |
| `Disk[{x,y}, r]` | Filled disk |
| `Rectangle[{x1,y1},{x2,y2}]` | Axis-aligned filled rectangle |
| `Polygon[{{x1,y1},…}]` | Filled polygon (any winding order) |
| `Text[expr, {x,y}]` | String / expression label (Hershey vector font) |

### Style directives

| Directive | Effect |
|-----------|--------|
| `RGBColor[r,g,b]` / `RGBColor[r,g,b,a]` | Set fill/stroke color |
| `GrayLevel[g]` / `GrayLevel[g,a]` | Grayscale color |
| `Hue[h]` / `Hue[h,s,b]` / `Hue[h,s,b,a]` | HSB color |
| `CMYKColor[c,m,y,k]` | CMYK color |
| `Opacity[a]` | Fill opacity `[0,1]` |
| `Thickness[t]` | Line thickness (plot coords) |
| `PointSize[s]` | Point radius (plot coords) |
| 24 named constants | `Red`, `Blue`, `Green`, `Black`, `White`, `Gray`, `Orange`, `Purple`, … |

### ColorFunction ramps (built-in named strings)

| Name | Aliases | Description |
|------|---------|-------------|
| `"Rainbow"` | — | Full HSV sweep (hue 0→1) |
| `"Temperature"` | `"Thermal"` | Dark blue-purple → orange → bright yellow |
| `"CoolTones"` | `"Cool"` | Ice blue → cornflower → deep navy |
| `"WarmTones"` | `"Warm"` | Pale cream → amber → deep crimson |
| `"Greyscale"` | `"Grayscale"`, `"Grey"`, `"Gray"` | White → black |

### ScalingFunctions (axis transforms)

Available for `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot`, `StreamPlot`:

| Value | Transform | Tick labels |
|-------|-----------|-------------|
| `"Log"` | natural log (world = ln data) | decade + 2×, 5× sub-ticks |
| `"Log10"` | log base 10 | decade + 2×, 5× sub-ticks |
| `"Log2"` | log base 2 | decade + 2×, 5× sub-ticks |
| `"Reverse"` | negate (world = –data) | linear, reversed |
| `None` / `Automatic` | identity | standard linear |
| `{"sfx", "sfy"}` | per-axis | mixed |

---

## Graphics
A symbolic 2D graphics object.
- `Graphics[primitives]`: wraps a (possibly nested) list of primitives and
  style directives.
- `Graphics[primitives, opt1 -> val1, ...]`: as above, with rendering
  options (see `Show` below for the option list).
- Like `Point`/`Line`/etc., `Graphics` is an inert structural head -- it
  does not evaluate or validate its contents; rendering happens only when
  `Show[]`/`Plot[]` is called on it.

**Features**:
- `Protected`.
- Prints as `-Graphics-` (use `FullForm[]` to inspect the underlying
  expression).
- `Options[Graphics]` reports the rendering options the engine honours (the
  set read by the renderer's `gfx_options_parse`): `AspectRatio`, `Axes`,
  `AxesLabel`, `AxesOrigin`, `AxesStyle`, `Background`, `Epilog`, `Frame`,
  `FrameLabel`, `FrameStyle`, `FrameTicks`, `GridLines`, `GridLinesStyle`,
  `ImageSize`, `LabelStyle`, `PlotLabel`, `PlotRange`, `PlotRangePadding`,
  `PlotStyle`, `Prolog`, `RotateLabel`, `TicksStyle`. See `Show` below for
  each option's effect.

```mathematica
In[1]:= Graphics[{Point[{0, 0}]}]
Out[1]= -Graphics-

In[2]:= Options[Graphics, {Axes, PlotRange}]
Out[2]= {Axes -> False, PlotRange -> Automatic}

In[2]:= FullForm[Graphics[{Point[{0, 0}]}]]
Out[2]= Graphics[List[Point[List[0, 0]]]]
```

## Point
A graphics primitive.
- `Point[{x, y}]`: a single point.
- `Point[{{x1,y1}, {x2,y2}, ...}]`: a list of points.

## Line
A graphics primitive.
- `Line[{{x1,y1}, {x2,y2}, ...}]`: a polyline through the given points.
- `Line[{{{...}}, {{...}}}]`: several disjoint polylines.

## Rectangle
A graphics primitive: `Rectangle[{xmin,ymin}, {xmax,ymax}]` is an
axis-aligned filled rectangle.

## Circle
A graphics primitive: `Circle[{x,y}, r]` is a circle outline of radius `r`
centered at `{x,y}`. The radius and centre carry Mathematica's defaults:
`Circle[{x,y}]` is radius 1, and `Circle[]` is the unit circle at the
origin.

## Disk
A graphics primitive: `Disk[{x,y}, r]` is a filled disk of radius `r`
centered at `{x,y}`. As with `Circle`, `Disk[{x,y}]` is radius 1 and
`Disk[]` is the unit disk at the origin.

## Polygon
A graphics primitive: `Polygon[{{x1,y1}, {x2,y2}, ...}]` is a filled,
closed polygon through the given vertices, in either winding order --
clockwise and counter-clockwise vertex lists both fill correctly (the
renderer detects the winding via the shoelace formula and corrects it
before drawing, since the underlying triangle-fan fill is
winding-sensitive but `Polygon[]` itself, like Mathematica's, is not).

## Text
A graphics primitive: `Text[expr, {x,y}]` renders `expr` (a string,
number, or any expression -- non-strings are stringified the same way
`ToString[]` would) as text centered at `{x,y}`, drawn with the classic
Hershey *Roman Simplex* single-stroke vector font (full printable ASCII,
upper- and lowercase, with proportional advance widths). The glyph data is
transcribed from the historical Hershey dataset (`tools/hershey.dat`) by
`tools/gen_hershey.py` into `src/graphics/hershey_glyphs.inc`.

## RGBColor / GrayLevel / Opacity / Thickness / PointSize
Style directives. Placed alongside primitives in a `Graphics[]` primitive
list, each affects every primitive that follows it (left to right), exactly
like Mathematica's directive semantics. See also `Hue` and `CMYKColor` below,
and the **Named color constants** section above for the 24 pre-defined color
symbols.

- `RGBColor[r, g, b]` / `RGBColor[r, g, b, a]`: color, components in `[0,1]`.
- `GrayLevel[g]` / `GrayLevel[g, a]`: a shade of gray (0 = black, 1 = white).
- `Opacity[a]`: fill opacity in `[0,1]` (does not affect line alpha; use the
  4-argument `RGBColor`/`GrayLevel`/`Hue` forms for full RGBA control).
- `Thickness[t]`: line thickness (plot coordinates) for `Line`/`Circle`.
- `PointSize[s]`: point radius (plot coordinates) for `Point`.

```mathematica
In[1]:= Graphics[{RGBColor[1, 0, 0], Point[{0, 0}], Line[{{0,0},{1,1}}]}]
Out[1]= -Graphics-

In[2]:= Graphics[{GrayLevel[0.3], Disk[{0, 0}], GrayLevel[0.8], Disk[{2, 0}]}]
Out[2]= -Graphics-

In[3]:= Graphics[{Blue, Opacity[0.5], Polygon[{{0,0},{1,0},{0.5,1}}],
                  Red,  Opacity[0.5], Polygon[{{0.5,0},{1.5,0},{1,1}}]}]
Out[3]= -Graphics-  (* two overlapping semi-transparent triangles *)
```

## Named color constants

All 24 named color constants are protected `OwnValue`s that evaluate to an
`RGBColor[r,g,b]` or `GrayLevel[g]` literal. They are usable anywhere a
color literal is accepted: as a style directive inside a primitive list,
as the value of `PlotStyle`, `Background`, `FrameStyle`, `AxesStyle`,
`TicksStyle`, `GridLinesStyle`, `ColorFunction` return values, and so on.

Because `Graphics[...]` is not held, its arguments evaluate at construction
time — a bare `Red` inside a primitive list is immediately replaced by
`RGBColor[1, 0, 0]`. `Plot`/`ParametricPlot`/etc. are `HoldAll` but
evaluate each option value once before storing it on the result, so a
named color in any option (e.g. `PlotStyle -> Orange` or
`Epilog -> {Blue, Disk[{0,0}]}`) resolves identically.

All four color forms — `RGBColor`, `GrayLevel`, `Hue`, `CMYKColor` — are
recognised as style directives in both the Raylib CLI renderer and the
Plotly notebook frontend.

**Color table** (24 constants):

| Name | Evaluates to | Approx. swatch |
|---|---|---|
| `Black` | `GrayLevel[0]` | black |
| `White` | `GrayLevel[1]` | white |
| `Gray` | `GrayLevel[0.5]` | mid-gray |
| `LightGray` | `GrayLevel[0.85]` | pale gray |
| `Red` | `RGBColor[1, 0, 0]` | pure red |
| `Green` | `RGBColor[0, 1, 0]` | lime green |
| `Blue` | `RGBColor[0, 0, 1]` | pure blue |
| `Cyan` | `RGBColor[0, 1, 1]` | cyan |
| `Magenta` | `RGBColor[1, 0, 1]` | magenta |
| `Yellow` | `RGBColor[1, 1, 0]` | yellow |
| `Brown` | `RGBColor[0.6, 0.4, 0.2]` | brown |
| `Orange` | `RGBColor[1, 0.5, 0]` | orange |
| `Pink` | `RGBColor[1, 0.5, 0.5]` | pink |
| `Purple` | `RGBColor[0.5, 0, 0.5]` | purple |
| `LightRed` | `RGBColor[1, 0.85, 0.85]` | pale red |
| `LightGreen` | `RGBColor[0.88, 1, 0.88]` | pale green |
| `LightBlue` | `RGBColor[0.87, 0.94, 1]` | pale blue |
| `LightCyan` | `RGBColor[0.9, 1, 1]` | pale cyan |
| `LightMagenta` | `RGBColor[1, 0.9, 1]` | pale magenta |
| `LightYellow` | `RGBColor[1, 1, 0.85]` | pale yellow |
| `LightBrown` | `RGBColor[0.94, 0.91, 0.88]` | pale brown |
| `LightOrange` | `RGBColor[1, 0.9, 0.8]` | pale orange |
| `LightPink` | `RGBColor[1, 0.925, 0.925]` | pale pink |
| `LightPurple` | `RGBColor[0.94, 0.88, 0.94]` | pale purple |

```mathematica
(* Named colors evaluate to their literal form *)
In[1]:= Red
Out[1]= RGBColor[1, 0, 0]

In[2]:= Black
Out[2]= GrayLevel[0]

In[3]:= {Red, Green, Blue, Black, White, Gray}
Out[3]= {RGBColor[1, 0, 0], RGBColor[0, 1, 0], RGBColor[0, 0, 1],
         GrayLevel[0], GrayLevel[1], GrayLevel[0.5]}

(* Style directive inside a primitive list *)
In[4]:= Graphics[{Red, Disk[{0, 0}], Blue, Disk[{2, 0}], Green, Disk[{4, 0}]}]
Out[4]= -Graphics-

(* As PlotStyle / Background options *)
In[5]:= Plot[Sin[x], {x, 0, 2 Pi}, PlotStyle -> Orange, Background -> LightGray]
Out[5]= -Graphics-

(* Per-color Disk grid — all 24 named colors *)
In[6]:= With[{cols = {Red, Green, Blue, Black, White, Gray,
                       Cyan, Magenta, Yellow, Brown, Orange, Pink,
                       Purple, LightRed, LightGreen, LightBlue, LightGray,
                       LightCyan, LightMagenta, LightYellow, LightBrown,
                       LightOrange, LightPink, LightPurple}},
         Table[{cols[[n]], Disk[{Mod[n - 1, 6] * 2, -Floor[(n - 1) / 6] * 2}]},
               {n, Length[cols]}]] // Graphics
Out[6]= -Graphics-

(* Named colors in multi-curve PlotStyle *)
In[7]:= Plot[{Sin[x], Cos[x], Sin[2 x]}, {x, 0, 2 Pi},PlotStyle -> {Red, Blue, Green}]
Out[7]= -Graphics-

(* Named colors in Epilog *)
In[8]:= Plot[Sin[x], {x, 0, 2 Pi},
          Epilog -> {Orange, PointSize[0.02], Point[{Pi/2, 1}],
                     Purple, PointSize[0.02], Point[{3 Pi/2, -1}]}]
Out[8]= -Graphics-

(* AxesStyle / GridLinesStyle *)
In[9]:= Plot[x^2, {x, -2, 2}, GridLines -> Automatic,
          GridLinesStyle -> LightBlue, AxesStyle -> Gray]
Out[9]= -Graphics-

(* Combining named colors with Opacity *)
In[10]:= Graphics[{Red, Opacity[0.4], Disk[{0, 0}],
                   Blue, Opacity[0.4], Disk[{1, 0}],
                   Green, Opacity[0.4], Disk[{0.5, 0.866}]}]
Out[10]= -Graphics-  (* three overlapping translucent disks *)

(* ParametricPlot with named color *)
In[11]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}, PlotStyle -> Purple]
Out[11]= -Graphics-

(* Plot3D with named color surface *)
In[12]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}, PlotStyle -> Cyan, Mesh -> None]
Out[12]= -Graphics3D-
```

## Named ColorFunction ramps

All gradient-based plotters — `DensityPlot`, `ContourPlot`, `VectorPlot`,
`StreamPlot`, `Plot`, `Plot3D`, `ParametricPlot`, `ParametricPlot3D` — accept
`ColorFunction -> "name"` where `name` is one of the following built-in
string ramps.  Each ramp is a continuous 5-stop RGB gradient parameterised
by `t ∈ [0,1]` (normalised to the data range when
`ColorFunctionScaling -> True`, the default).

| Name | Aliases | Appearance |
|------|---------|------------|
| `"Rainbow"` | — | `Hue` sweep red → violet (stops at 0.8 to avoid wrapping back to red) |
| `"Temperature"` | `"Thermal"` | dark blue-purple → purple → red → orange → bright yellow (Mathilda's default stream/contour ramp) |
| `"CoolTones"` | `"Cool"` | near-white ice blue → sky blue → cornflower → deep navy/indigo |
| `"WarmTones"` | `"Warm"` | pale cream → amber → orange → deep crimson |
| `"Greyscale"` | `"Grayscale"`, `"Grey"`, `"Gray"` | white (t=0) → black (t=1) |

**Calling convention by plotter:**

| Plotter | `t` is derived from… |
|---------|----------------------|
| `Plot`, `ParametricPlot`, `PolarPlot` | position along curve: `(x − xmin)/(xmax − xmin)` |
| `Plot3D`, `ParametricPlot3D` | z-height: `(z − zmin)/(zmax − zmin)` |
| `DensityPlot`, `ContourPlot` | cell's normalised function value |
| `VectorPlot`, `StreamPlot` | field magnitude (speed) |

For custom function forms see the individual plotter's `ColorFunction` option
table.

```mathematica
(* DensityPlot with each named ramp *)
In[1]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3},
          ColorFunction -> "CoolTones"]
Out[1]= -Graphics-

In[2]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3},
          ColorFunction -> "WarmTones"]
Out[2]= -Graphics-

In[3]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3},
          ColorFunction -> "Greyscale"]
Out[3]= -Graphics-

In[4]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3},
          ColorFunction -> "Temperature"]
Out[4]= -Graphics-

(* ContourPlot with a warm ramp *)
In[5]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2},
          ContourShading -> True, ColorFunction -> "WarmTones"]
Out[5]= -Graphics-

(* VectorPlot with cool tones *)
In[6]:= VectorPlot[{-y, x}, {x, -2, 2}, {y, -2, 2},
          ColorFunction -> "CoolTones"]
Out[6]= -Graphics-

(* StreamPlot with greyscale speed coloring *)
In[7]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2},
          StreamColorFunction -> "Greyscale", PlotLegends -> Automatic]
Out[7]= -Graphics-
```

## Hue
A style directive: `Hue[h]` (fully saturated, `h` in `[0,1]`, wrapping),
`Hue[h, s, b]` (saturation, brightness), `Hue[h, s, b, a]` (with opacity) --
standard HSB-to-RGB, recognized everywhere `RGBColor`/`GrayLevel` are
(style directives, `PlotStyle`, `Background`, `FrameStyle`, `AxesStyle`,
`TicksStyle`, `GridLinesStyle`, `ColorFunction`'s return value).

## CMYKColor
A style directive for the subtractive (print) color model used by inks.

- `CMYKColor[c, m, y, k]`: cyan, magenta, yellow, black components in `[0,1]`.
- `CMYKColor[c, m, y]`: equivalent to `k = 0`.
- `CMYKColor[c, m, y, k, a]`: with opacity `a`.
- `CMYKColor[{c, m, y, k}]` / `CMYKColor[{c, m, y, k, a}]`: list forms.

Components and opacity outside `[0,1]` are clipped. At render time the color
is converted to RGB by `r = (1-c)(1-k)`, `g = (1-m)(1-k)`, `b = (1-y)(1-k)`,
so `CMYKColor` is recognized everywhere `RGBColor`/`GrayLevel`/`Hue` are.

```mathematica
In[1]:= Graphics[{CMYKColor[1, 0, 0, 0], Disk[]}]   (* a cyan disk *)
```

## Show
Normalizes a `Graphics[...]` object for display; the REPL front end opens
the window when the result is a top-level `Graphics` (see the auto-display
note at the top). The window blocks the REPL until closed.
- `Show[graphics]`: returns `graphics` (a copy).
- `Show[graphics, opt -> val, ...]`: merges the given options into
  `graphics`'s option list (a later/explicit option overrides one already
  present) and returns the merged `Graphics[...]`.

Rendering options (read by `Show`/`Plot`'s renderer, stored as trailing
`Rule[...]` arguments on the `Graphics[...]` object):

| Option | Default (bare `Graphics`) | Default (via `Plot`) | Meaning |
|---|---|---|---|
| `Axes` | `False` | `True` (unless `Frame` on) | draw coordinate axes with tick labels |
| `Frame` | `False` | `False` | draw a tick-and-label frame around the plot: `True`/`None`/`False`, or per-edge `{{left,right},{bottom,top}}` of `True`/`False` |
| `FrameTicks` | `Automatic` | `Automatic` | `Automatic`/`None`, or per-edge `{{left,right},{bottom,top}}` of `Automatic`/`None`; ticks appear only on drawn frame edges |
| `FrameStyle` | gray | gray | `RGBColor`/`GrayLevel` for the frame box, ticks and labels |
| `AspectRatio` | `Automatic` (true geometry) | `1/GoldenRatio` | window height-to-width ratio: `Automatic` shapes it from the data extent, `Full` keeps the `ImageSize` box and stretches the data to fill it, or an explicit ratio `a` (incl. symbolic like `1/GoldenRatio`) |
| `PlotRange` | `Automatic` (data bbox, y spike-clipped) | same | `Automatic` (default), `All` (no clipping), `{{xmin,xmax},{ymin,ymax}}` (fix both axes), or `{ymin,ymax}` (fix y only, x stays automatic) |
| `PlotStyle` | a default blue | same | style directive(s) used as the initial draw color |
| `Background` | white | white | window background color |
| `ImageSize` | width `800` | same | `w` sets the width (height follows from `AspectRatio`), or `{w, h}` fixes both in pixels (then `AspectRatio` shapes the data inside the fixed box) |
| `AxesLabel` | none | none | `{xlabel, ylabel}` |
| `PlotLabel` | none | none | title drawn above the plot |
| `AxesOrigin` | `Automatic` (0 if in range, else clamp to the near edge) | same | `{x, y}` pins exactly where the axes cross |
| `AxesStyle` | gray | gray | `RGBColor`/`GrayLevel` for the axis lines and tick marks |
| `TicksStyle` | gray | gray | `RGBColor`/`GrayLevel` for the axis tick *labels* (independent of `AxesStyle`) |
| `FrameLabel` | none | none | `{xlabel, ylabel}`, the `Frame`-mode equivalent of `AxesLabel`: xlabel below the frame, ylabel beside it |
| `RotateLabel` | `True` | `True` | whether `FrameLabel`'s y-label is rotated 90° (`True`) or horizontal (`False`) |
| `PlotRangePadding` | `0.08` (8%) | same | extra fraction of the data extent added on each side during auto-fit: a number (both axes), `{xfrac, yfrac}`, or `None` for no padding |
| `GridLines` | `None` | `None` | `None` (no grid), `Automatic` (light lines at the same major ticks as `Axes`/`Frame`), or `{xspec, yspec}` with each independently `None`/`Automatic`/an explicit `List[...]` of positions |
| `GridLinesStyle` | light gray | light gray | `RGBColor`/`GrayLevel` for the grid lines |
| `Prolog` | none | none | extra primitive(s), drawn first in data coordinates (under the grid/axes/curve) |
| `Epilog` | none | none | extra primitive(s), drawn last in data coordinates (on top of the curve; still pans/zooms with it) |
| `LabelStyle` | none | none | `RGBColor`/`GrayLevel`/`Hue` fallback seeded into `AxesStyle`/`TicksStyle`/`FrameStyle` before any of those are individually resolved -- a specific one still overrides it, regardless of argument order |

**Automatic vertical range** — with `PlotRange -> Automatic` (the default),
the y-extent is chosen to frame the visible body of the curve rather than its
raw min/max. A function with an asymptote (e.g. `Tan`, `Sec`, `1/x`, `Gamma`)
climbs to enormous finite values near each pole; taking those at face value
would crush the interesting part of the curve to a flat line (`Plot[Tan[x],
{x, -2 Pi, 2 Pi}]` would span ±10^16). The renderer detects such a curve by a
near-vertical *runaway* segment — one whose slope is ~1000× the curve's median
slope, the signature of a pole — and, only then, clips the y-range to a robust
median/MAD band so the spikes run off-screen. Smooth curves (including steep
ones like `x^5` or `Exp`, and curves with a vertical tangent like `Sqrt` or a
removable singularity like `Sinc`) have no such runaway and keep their full,
legitimate extent. `PlotRange -> All` disables the clipping and shows every
sampled point; an explicit `{ymin, ymax}` overrides it outright.

**Window shape (`AspectRatio`)** — `AspectRatio` is the plot's height-to-width
ratio and so sizes the *window itself*, not merely the y-scale inside a fixed
box. A default `Plot` (`AspectRatio -> 1/GoldenRatio`) opens a short, wide
window (e.g. `800 × 494`) with the curve filling it edge to edge. The window
height is `round(width * ratio)`, clamped to a screen-friendly `[100, 2000]`
band. `Automatic` takes the ratio from the data extent (so `Graphics[Circle[]]`
is square); `Full` instead keeps the `ImageSize` box and stretches the data to
fill it; and `ImageSize -> {w, h}` pins both dimensions, with `AspectRatio` then
shaping the data within the fixed frame. The value may be any real, including a
symbolic constant such as `1/GoldenRatio` or `GoldenRatio`.

**Frame (`Frame`)** — `Frame -> True` rules the plot with a rectangle along its
edges, replacing the through-the-origin `Axes` cross (in `Plot`, supplying
`Frame -> True` withholds the default `Axes -> True`; pass `Axes -> True`
explicitly to keep both). Each edge can be toggled independently with
`Frame -> {{left, right}, {bottom, top}}`. The frame reserves a margin (~5% of
each window dimension; the bottom and left a little more to hold the numeric
labels and the help line) and the curve is **fitted to and clipped against** the
interior, so it never spills past the box. Frame ticks come in two tiers: major
ticks land on the same "nice" values as the axes, and minor sub-ticks subdivide
each major interval — the number of subdivisions is read off the major step's
leading digit so the minors always fall on round values (a step of `1` splits
into 5, `2` into 4, `5` into 5). Minor ticks are half-length and unlabeled. Ticks
point inward; the major-tick numeric labels sit just *outside* the frame in the
reserved margin, on the bottom and left edges (or the top/right when those are
the only ones drawn). The frame box, its ticks and its labels render at 1.5 px;
the `Axes` lines, ticks and labels render at the same 1.5 px, so the plot looks
identical in weight whether it is framed or axed. `FrameTicks -> None` keeps the
box but drops all ticks and labels; `FrameStyle` colors the box, ticks and labels.

**Interactive controls** (within the window): mouse wheel to zoom *about
the cursor* (the point under the mouse stays fixed, so you magnify whatever
you point at rather than always diving toward the plot centre), right-drag
(or middle-drag) to pan, `Q`/`E` to rotate the view, `R` to reset, `Esc` or
the window's close button to return to the REPL.

**Toolbar** — a Plotly-style modebar of gray icon buttons sits in the
top-right corner (hover for a tooltip):

| Icon | Action |
|---|---|
| Camera | Download the plot as `mathilda_plot.png` (UI chrome excluded) in the current directory |
| Magnifier | Zoom tool: left-drag a box to zoom to that region |
| Move arrows | Pan tool: left-drag to move (the default active tool) |
| Magnifier `+` | Zoom in about the center |
| Magnifier `−` | Zoom out about the center |
| Expand arrows | Autoscale to fit the data |
| House | Reset axes to the initial view |
| ✕ | Close the window (red hover tint; same as `Esc`) |

The active drag tool (Pan or Zoom) is highlighted. Left-drag follows the
selected tool; right/middle-drag always pans regardless of tool. All icons
are anti-aliased hand-drawn vector glyphs (4× MSAA).

**Live re-sampling** — for a window created by `Plot`, zooming or panning
re-runs the adaptive sampler over the newly visible x-range, so a magnified
curve keeps full resolution (e.g. `Plot[Sin[1/x^2], {x, 1, 3}]` stays
smooth when zoomed in) rather than exposing the home view's coarse grid.
Re-sampling lands on drag release to keep gestures smooth. Because the
mouse wheel zooms about the cursor, you can drive the view into an
off-centre region of interest and watch the sampler refine it — handy when
the plot's centre is its worst point (the `x=0` essential singularity of
`Sin[1/x^2]`, say, which no sampling can resolve).

**Features**:
- `Protected`.
- Declines to evaluate (stays unevaluated) if its first argument isn't a
  `Graphics[...]` expression, or any trailing argument isn't a `Rule`.
- When Raylib isn't compiled in, prints a one-line message instead of
  opening a window; the option merge still happens and the merged
  `Graphics[...]` is still returned.

```mathematica
In[1]:= Show[Graphics[{Point[{0,0}]}], Axes -> True]
Out[1]= -Graphics-
```

## Plot
Adaptively samples and displays a function of one real variable.
- `Plot[f, {x, xmin, xmax}]`: samples `f` over `[xmin, xmax]`, opens an
  interactive window showing the curve, and returns it as a
  `Graphics[{Line[...], ...}, opts]` object. A gap in `f` (a singularity,
  or a non-real result) breaks the curve into separate `Line[...]`
  segments rather than bridging it.
- `Plot[{f1, f2, ...}, {x, xmin, xmax}]`: plots several curves on the same
  axes. Each is sampled independently (its own gaps/singularities) and
  drawn in a distinct colour drawn from the default plot palette
  (Mathematica's `ColorData[97]`), chosen so the series read as distinct
  yet harmonious; the palette cycles when there are more curves than
  entries. A single curve keeps the `PlotStyle` colour.
- `Plot[f, {x, xmin, xmax}, opts...]`: as above, with sampling and/or
  rendering options (see below and the `Show` table above).
- `HoldAll`: `f` (or the function list) and the iterator spec are not
  pre-evaluated (`x` has no value until sampling binds it).

Options: see **Feature summary** above. `HoldAll`, `Protected`. Unknown options pass through to `Graphics[...]`.

```mathematica
In[1]:= Plot[Sin[x], {x, 0, 2 Pi}]
Out[1]= -Graphics-

In[2]:= Plot[1/x, {x, -2, 2}]
Out[2]= -Graphics-

In[3]:= Plot[Sin[x], {x, a, b}]
Out[3]= Plot[Sin[x], {x, a, b}]

In[4]:= Plot[{Sin[1/x], Cos[1/x]}, {x, -Pi, Pi}]
Out[4]= -Graphics-

In[5]:= Plot[Sin[x] + Sin[7 x], {x, -2, 2}, Mesh -> All]
Out[5]= -Graphics-

In[6]:= Plot[Sin[x], {x, 0, 2 Pi}, GridLines -> Automatic,
            Epilog -> {Red, Point[{0, 0}]}, AxesOrigin -> {0, 0}]
Out[6]= -Graphics-

In[7]:= Plot[Sin[x], {x, 0, 2 Pi}, ColorFunction -> "Rainbow", Filling -> Axis]
Out[7]= -Graphics-

In[8]:= Plot[{Sin[x], Cos[x]}, {x, 0, 2 Pi}, PlotLegends -> "Expressions"]
Out[8]= -Graphics-
```

## Graphics3D
A symbolic 3D graphics object, as built by `Plot3D[]`.
- `Graphics3D[primitives]`: wraps a (possibly nested) list of primitives and
  style directives -- the *same* heads `Graphics[]` uses (`Polygon`, `Line`,
  `RGBColor`/`GrayLevel`/`Hue`, `Opacity`), but with 3-coordinate `{x,y,z}`
  points instead of 2-coordinate `{x,y}` ones, exactly mirroring how
  Mathematica itself reuses primitive heads between `Graphics`/`Graphics3D`.
- `Graphics3D[primitives, opt1 -> val1, ...]`: as above, with options (see
  `Plot3D` below).
- Like `Graphics`, it is an inert structural head -- nothing renders until
  `Show[]`/`Plot3D[]` is called on it.
- Rendered in an orbit-camera window: drag to rotate, scroll to zoom,
  right-drag (or middle-drag) to pan, `R` to reset the view, `S` to save a
  screenshot, `Esc` to close. There is no live re-sampling on
  rotate/zoom -- unlike a 2D pan/zoom, orbiting the camera never changes
  which `(x,y)` domain was sampled, so there is nothing to re-sample.

**Features**:
- `Protected`.
- Prints as `-Graphics3D-` (use `FullForm[]` to inspect the underlying
  expression).

```mathematica
In[1]:= Graphics3D[{Polygon[{{0,0,0}, {1,0,0}, {1,1,1}, {0,1,1}}]}]
Out[1]= -Graphics3D-
```

## Plot3D
Samples and displays a function of two real variables as a surface.
- `Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}]`: samples `f` over a uniform
  grid on `[xmin,xmax] x [ymin,ymax]`, opens an orbit-camera window showing
  the surface, and returns it as a `Graphics3D[{Polygon[...], ...}, opts]`
  object. Each grid cell becomes one quad `Polygon[]`; a cell with any
  corner that doesn't evaluate to a finite real (or fails `RegionFunction`)
  is simply omitted, leaving a hole in the surface.
- `Plot3D[{f1, f2, ...}, {x,...}, {y,...}]`: plots several surfaces over the
  same domain, each in a distinct palette colour (the same `ColorData[97]`
  palette `Plot` uses).
- `Plot3D[f, {x,...}, {y,...}, opts...]`: as above, with options below.
- `HoldAll`: `f` and both iterator specs are not pre-evaluated.

Options: see **Feature summary** above. `HoldAll`, `Protected`. `ExclusionStyle` (default `GrayLevel[0.35]`) styles boundary edges when `RegionFunction` is active.

```mathematica
In[1]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[1]= -Graphics3D-

In[2]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
Out[2]= -Graphics3D-  (* Rainbow maps z-height to hue *)

In[3]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
Out[3]= -Graphics3D-  (* exclusion boundary drawn in default GrayLevel[0.35] *)

In[4]:= Plot3D[x + y, {x,-2,2}, {y,-2,2},RegionFunction -> Function[{x,y,z}, x^2+y^2 <4],ExclusionStyle -> RGBColor[1, 0.3, 0]]
Out[4]= -Graphics3D-  (* exclusion boundary in orange *)

In[5]:= Plot3D[{Sin[x + y], Cos[x - y]}, {x, -2, 2}, {y, -2, 2}]
Out[5]= -Graphics3D-  (* palette colors for each surface *)

In[6]:= Plot3D[{x^2, x^2 + 1}, {x,-2,2}, {y,-2,2},PlotStyle -> {Blue, Red}]
Out[6]= -Graphics3D-  (* explicit per-surface colors *)
```

## ParametricPlot
Samples and displays a parametric planar curve or filled region.

**One-iterator form** — parametric curve:
- `ParametricPlot[body, {t, tmin, tmax}]`: evaluates `body` (which must
  produce a 2-element numeric list `{x, y}`) for each sampled `t` value;
  returns `Graphics[{Line[...], ...}, opts]`. `body` may be a literal `{fx,
  fy}` or any expression that evaluates to `{x, y}` (e.g. `r {Cos[t], Sin[t]}`
  with `r` set to a number elsewhere).
- `ParametricPlot[{{body1}, {body2}, ...}, {t, tmin, tmax}]`: plots multiple
  curves over the same parameter range in distinct palette colours.

**Two-iterator form** — filled parametric region:
- `ParametricPlot[body, {t, tmin, tmax}, {r, rmin, rmax}]`: samples a
  `PlotPoints × PlotPoints` grid of `(t, r)` pairs, maps each pair to `(x, y)`
  via `body`, and emits filled `Polygon[{p00,p10,p11,p01}]` quads for each
  valid grid cell.  Produces `Graphics[{Polygon[...], ...}, opts]`.  `Mesh ->
  All` overlays the grid lines.  The fill is **solid** by default; use
  `PlotStyle -> {color, Opacity[a]}` to get a transparent region.

`HoldAll`, `Protected`. Options: see **Feature summary** above.

```mathematica
(* --- One-iterator: curves --- *)
In[1]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[1]= -Graphics-  (* unit circle, AspectRatio -> 1 *)

In[2]:= ParametricPlot[{Sin[2 t], Sin[3 t]}, {t, 0, 2 Pi}]
Out[2]= -Graphics-  (* Lissajous figure *)

In[3]:= ParametricPlot[{{Cos[t], Sin[t]}, {2 Cos[t], Sin[t]}}, {t, 0, 2 Pi}]
Out[3]= -Graphics-  (* two curves in palette colours *)

In[4]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},
          ColorFunction -> (Hue[#] &)]
Out[4]= -Graphics-  (* rainbow-colored circle *)

In[5]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},
          RegionFunction -> Function[{x, y}, x > 0]]
Out[5]= -Graphics-  (* right semicircle only *)

In[6]:= ParametricPlot[2 {Cos[t], Sin[t]}, {t, 0, 2 Pi}]
Out[6]= -Graphics-  (* computed body (not a literal List) -- works fine *)

(* --- Two-iterator: filled regions --- *)
In[7]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}]
Out[7]= -Graphics-  (* annular region, r from 1 to 2 *)

In[8]:= ParametricPlot[r^2 {Sqrt[t] Cos[t], Sin[t]},
          {t, 0, 3 Pi/2}, {r, 1, 2}]
Out[8]= -Graphics-  (* weighted spiral region *)

In[9]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2},
          Mesh -> All]
Out[9]= -Graphics-  (* with grid lines overlaid *)
```

## PolarPlot
Plots one or more curves given in polar coordinates `r(theta)`, converting
to Cartesian `{r*Cos[theta], r*Sin[theta]}` and sampling adaptively.  Returns a
`Graphics[...]` object (auto-displayed).

- `PolarPlot[r, {theta, tmin, tmax}]` — single polar curve.
- `PolarPlot[{r1, r2, ...}, {theta, tmin, tmax}]` — multiple curves in
  distinct palette colours.

`HoldAll`, `Protected`. Negative `r` values plot in the opposite direction. Default `AspectRatio -> 1`. Options: same as `ParametricPlot`; see **Feature summary** above. `PolarAxes -> True` is accepted but polar grid overlay is not yet rendered.

```mathematica
In[1]:= PolarPlot[1, {t, 0, 2Pi}]
Out[1]= -Graphics-  (* unit circle *)

In[2]:= PolarPlot[Sin[2t], {t, 0, 2Pi}]
Out[2]= -Graphics-  (* four-petal rose *)

In[3]:= PolarPlot[{1, 2}, {t, 0, 2Pi}, PlotStyle -> {Blue, Red}]
Out[3]= -Graphics-  (* two concentric circles *)

In[4]:= PolarPlot[t, {t, 0, 4Pi}, ColorFunction -> "Rainbow"]
Out[4]= -Graphics-  (* Archimedean spiral, rainbow-coloured *)

In[5]:= PolarPlot[Sin[2t], {t, 0, 2Pi}, Mesh -> All, PlotLabel -> "Rose"]
Out[5]= -Graphics-  (* rose with sample-point overlay and title *)
```

## PolarAxes
`PolarPlot` option: `True` requests a polar grid overlay (radial circles at
regular intervals plus angular degree/radian labels). Currently accepted and
documented but not yet rendered; Cartesian axes are drawn instead.

## ParametricPlot3D
Samples and displays a parametric 3D space curve or surface patch, returning
a `Graphics3D[...]` object rendered in an orbit-camera window.

**One-iterator form** — parametric space curve:
- `ParametricPlot3D[body, {t, tmin, tmax}]`: evaluates `body` (which must
  produce a 3-element numeric list `{x, y, z}`) for each sampled `t` value;
  returns `Graphics3D[{Line[...], ...}, opts]`. `body` may be a literal
  `{fx, fy, fz}` or any expression that evaluates to `{x, y, z}`.
- `ParametricPlot3D[{{body1}, {body2}, ...}, {t, tmin, tmax}]`: plots
  multiple space curves over the same parameter range in distinct palette
  colours.

**Two-iterator form** — parametric 3D surface patch:
- `ParametricPlot3D[body, {t, tmin, tmax}, {u, umin, umax}]`: samples a
  `PlotPoints × PlotPoints` grid of `(t, u)` pairs, maps each to `{x, y, z}`
  via `body`, and emits filled `Polygon[{p00,p10,p11,p01}]` quads for each
  valid grid cell. Produces `Graphics3D[{Polygon[...], ...}, opts]`.

`HoldAll`, `Protected`. Options: see **Feature summary** above. `ColorFunction` receives scaled spatial `{xs,ys,zs}` coordinates (not parameter values). `RegionFunction` is tried as `f[x,y,z]` first, then `f[x,y]`. Interactive controls same as `Plot3D`.

```mathematica
(* --- One-iterator: space curves --- *)
In[1]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}]
Out[1]= -Graphics3D-  (* helix *)

In[2]:= ParametricPlot3D[{Sin[2 t] Cos[t], Sin[2 t] Sin[t], Cos[t]}, {t, 0, 2 Pi}]
Out[2]= -Graphics3D-  (* spherical Lissajous *)

In[3]:= ParametricPlot3D[{{Cos[t], Sin[t], 0}, {Cos[t], 0, Sin[t]}}, {t, 0, 2 Pi}]
Out[3]= -Graphics3D-  (* two circles in different planes *)

In[4]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi},
          ColorFunction -> "Rainbow"]
Out[4]= -Graphics3D-  (* hue sweeps with z-height *)

In[5]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}, Mesh -> All]
Out[5]= -Graphics3D-  (* with sample dots *)

(* --- Two-iterator: surface patches --- *)
In[6]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]},
          {u, 0, 2 Pi}, {v, 0, Pi}]
Out[6]= -Graphics3D-  (* unit sphere *)

In[7]:= ParametricPlot3D[{(2 + Cos[v]) Cos[u], (2 + Cos[v]) Sin[u], Sin[v]},
          {u, 0, 2 Pi}, {v, 0, 2 Pi}]
Out[7]= -Graphics3D-  (* torus, R=2, r=1 *)

In[8]:= ParametricPlot3D[{u Cos[v], u Sin[v], u},
          {u, 0, 2}, {v, 0, 2 Pi},
          ColorFunction -> "Rainbow", Mesh -> All]
Out[8]= -Graphics3D-  (* cone with rainbow coloring and mesh *)

In[9]:= ParametricPlot3D[{{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]},
          {2 Cos[u] Sin[v], 2 Sin[u] Sin[v], 2 Cos[v]}},
          {u, 0, 2 Pi}, {v, 0, Pi}]
Out[9]= -Graphics3D-  (* two concentric spheres in palette colours *)

In[10]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]},
           {u, 0, 2 Pi}, {v, 0, Pi}, Lighting -> None]
Out[10]= -Graphics3D-  (* flat color, no shading *)
```

## ListPlot
Plots explicit data as a point (scatter) plot, returning a
`Graphics[{Point[...], ...}, opts]` object rendered through the same engine
as `Plot`/`Show`. Unlike `Plot` (which is `HoldAll` so its function body stays
symbolic), `ListPlot`'s data is concrete and **is** evaluated, so
`ListPlot[Range[10]]` and `ListPlot[Table[i^2, {i, 5}]]` work.

Data is classified into one or more datasets:

- `ListPlot[{y1, ..., yn}]` — a flat list of **heights**, plotted as the
  points `{i, yi}` (`x` running over `DataRange`, default `1` to `n`).
  Non-numeric entries are treated as missing and skipped (their index slot is
  still consumed).
- `ListPlot[{{x1, y1}, ..., {xn, yn}}]` — a list of **coordinate pairs**,
  plotted as the given points (a scatter plot).
- `ListPlot[{data1, data2, ...}]` — several **datasets** (any element that is
  itself a sublist and not a flat pair list), each drawn in a distinct palette
  colour (`ColorData[97]`, as for multi-curve `Plot`).

A non-list argument, an empty list, or all-missing data leaves `ListPlot`
unevaluated.

Options: see **Feature summary** above. Defaults: `Axes -> True`, `AspectRatio -> 1/GoldenRatio`. Marker size adapts to scatter density; an explicit `PointSize` overrides it.

```mathematica
In[1]:= ListPlot[{1, 4, 9, 16, 25}]
Out[1]= -Graphics-

In[2]:= ListPlot[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}]
Out[2]= -Graphics-

In[3]:= ListPlot[{{1, 1}, {2, 4}, {3, 9}}, Joined -> True]
Out[3]= -Graphics-

In[4]:= ListPlot[{Table[Sin[n], {n, 20}], Table[Cos[n], {n, 20}]},
            PlotLegends -> {"sin", "cos"}]
Out[4]= -Graphics-

In[5]:= ListPlot[{1, 4, 9, 16}, Filling -> Axis]
Out[5]= -Graphics-
```

## StreamPlot

Traces streamlines of a 2-D vector field `{vx, vy}` by RK4 integration
from a uniform grid of seed points, renders each stream as an
`Arrow[...]` primitive (a directed polyline with an arrowhead at its
end), and returns a `Graphics[...]` object (auto-displayed).

`HoldAll`, `Protected`. Options: see **Feature summary** above. Defaults: `Axes -> True`, `AspectRatio -> 1`. `ColorFunction` is an alias for `StreamColorFunction`.

```mathematica
StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}]
StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

```mathematica
(* --- Basic stream plots --- *)
In[1]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[1]= -Graphics-   (* circular rotation *)

In[2]:= StreamPlot[{1 - y^2, x}, {x, -3, 3}, {y, -2, 2}]
Out[2]= -Graphics-   (* nonlinear field *)

In[3]:= StreamPlot[{Sin[x + y], Cos[x - y]}, {x, 0, 2 Pi}, {y, 0, 2 Pi}]
Out[3]= -Graphics-

(* --- Denser seeding --- *)
In[4]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamPoints -> 25]
Out[4]= -Graphics-

(* --- StreamScale: let streams run freely --- *)
In[5]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamScale -> None]
Out[5]= -Graphics-

(* --- Speed-colored streams --- *)
In[6]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamColorFunction -> "Rainbow"]
Out[6]= -Graphics-

(* --- Custom style --- *)
In[7]:= StreamPlot[{x, -y}, {x, -2, 2}, {y, -2, 2}, StreamStyle -> {Thickness[0.004], RGBColor[0.8, 0.2, 0.1]}]
Out[7]= -Graphics-

(* --- RegionFunction: mask to a disk --- *)
In[8]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 1.5^2]]
Out[8]= -Graphics-

(* --- Arrow primitive directly --- *)
In[9]:= Graphics[{Blue, Arrow[{{0,0}, {1,0}, {1,1}}]}]
Out[9]= -Graphics-
```

## ContourPlot
Generates iso-contour lines of a 2-D function `f(x, y)` using the marching
squares algorithm and returns a `Graphics[...]` object (auto-displayed).

```mathematica
ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}]
ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

`HoldAll`, `Protected`. Uses marching squares with bilinear saddle-cell disambiguation. Options: see **Feature summary** above. Defaults: `Axes -> True`, `AspectRatio -> 1`.

```mathematica
(* --- Basic contour plots --- *)
In[1]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[1]= -Graphics-  (* 10 auto-levels, coloured by height *)

In[2]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, Contours -> 5]
Out[2]= -Graphics-  (* 5 circular contours *)

(* --- Explicit contour values --- *)
In[3]:= ContourPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, Contours -> {-2, -1, 0, 1, 2}]
Out[3]= -Graphics-  (* hyperbolas at specified levels *)

(* --- Shading with ColorFunction --- *)
In[4]:= ContourPlot[Sin[x + y], {x, -3, 3}, {y, -3, 3},
          ColorFunction -> "Rainbow", ContourShading -> True]
Out[4]= -Graphics-  (* rainbow-filled density with contour lines *)

In[5]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2},
          ContourShading -> True, Contours -> 8]
Out[5]= -Graphics-  (* thermal-gradient fill, 8 levels *)

(* --- Lines only --- *)
In[6]:= ContourPlot[Sin[x] Cos[y], {x, -Pi, Pi}, {y, -Pi, Pi},
          ContourShading -> False, ContourStyle -> {Thickness[0.006]},
          PlotPoints -> 40]
Out[6]= -Graphics-

(* --- Labels --- *)
In[7]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2},
          ContourLabels -> True, Contours -> 5]
Out[7]= -Graphics-  (* level values annotated at first segment *)

(* --- RegionFunction: circular mask --- *)
In[8]:= ContourPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3},
          RegionFunction -> Function[{x, y}, x^2 + y^2 < 4],
          ContourShading -> True]
Out[8]= -Graphics-

(* --- ContourStyle: cycle explicit colours --- *)
In[9]:= ContourPlot[Sin[x + y], {x, -3, 3}, {y, -3, 3},
          ContourStyle -> {Red, Blue, Green}, Contours -> 6]
Out[9]= -Graphics-  (* cycles Red, Blue, Green across 6 levels *)

(* --- Suppress lines, shading only --- *)
In[10]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3},
           ContourStyle -> None, ContourShading -> True,
           ColorFunction -> "Temperature"]
Out[10]= -Graphics-
```

---

## DensityPlot

```
DensityPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

Renders `f(x,y)` as a heatmap by shading each grid cell with the colour corresponding to its average `f` value. `HoldAll`, `Protected`. Default `ColorFunction`: thermal blue→yellow ramp. Options: see **Feature summary** above.

**Examples:**

```mathematica
(* Basic heatmap *)
In[1]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3}]
Out[1]= -Graphics-

(* Rainbow colour scheme *)
In[2]:= DensityPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2},
          ColorFunction -> "Rainbow", PlotPoints -> 60]
Out[2]= -Graphics-

(* Custom ColorFunction + legend *)
In[3]:= DensityPlot[Sin[x + y], {x, 0, 6}, {y, 0, 6},
          ColorFunction -> (GrayLevel[#]&), PlotLegends -> Automatic]
Out[3]= -Graphics-

(* RegionFunction: circular mask *)
In[4]:= DensityPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3},
          RegionFunction -> Function[{x,y}, x^2 + y^2 < 4]]
Out[4]= -Graphics-
```

---

## BarChart

```
BarChart[{v1, v2, ..., vn}, opts...]
BarChart[{{v1,...}, {w1,...}, ...}, opts...]
```

Draws a vertical bar chart. Single-dataset: `n` bars at `x = 1..n`. Multi-dataset: grouped bars in distinct palette colours. Options: `ChartStyle`, `ChartLabels`, `BarSpacing` (gap fraction, default `0.2`), standard `Graphics` options. Defaults: `Axes→True`, `AspectRatio→0.618`.

**Examples:**

```mathematica
(* Simple bar chart *)
In[1]:= BarChart[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[1]= -Graphics-

(* Custom colours and labels *)
In[2]:= BarChart[{2.5, 4.1, 3.3, 5.7},
          ChartStyle -> {Red, Blue, Green, Orange},
          ChartLabels -> {"Q1", "Q2", "Q3", "Q4"}]
Out[2]= -Graphics-

(* Grouped datasets *)
In[3]:= BarChart[{{1, 3, 2}, {4, 2, 5}}, BarSpacing -> 0.3]
Out[3]= -Graphics-

(* Negative values *)
In[4]:= BarChart[{3, -1, 4, -1, 5}]
Out[4]= -Graphics-
```

---

## Histogram

```
Histogram[data, opts...]
Histogram[data, k, opts...]
Histogram[data, {step}, opts...]
Histogram[data, {min, max, step}, opts...]
```

Bins numeric `data` and draws a frequency histogram. Default bin count: Sturges' rule (`⌈log₂(n)⌉ + 1`). Options: `ChartStyle`, `BarSpacing` (default `0.2`), standard `Graphics` options. Defaults: `Axes→True`, `AspectRatio→0.618`.

**Examples:**

```mathematica
(* Auto-binned histogram *)
In[1]:= Histogram[Table[RandomReal[], {200}]]
Out[1]= -Graphics-

(* Bind data first, then vary bin specs *)
In[2]:= data = Table[RandomReal[], {200}]

(* Explicit bin count *)
In[3]:= Histogram[data, 20]
Out[3]= -Graphics-

(* Fixed bin width *)
In[4]:= Histogram[data, {0.1}]
Out[4]= -Graphics-

(* Explicit range and width *)
In[5]:= Histogram[data, {0, 1, 0.05}]
Out[5]= -Graphics-
```

---

## VectorPlot

```
VectorPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

Draws a grid of arrows for the vector field `{vx, vy}`. `HoldAll`, `Protected`. Default `ColorFunction`: thermal ramp keyed to speed. Arrow sizing is screen-normalized so arrows stay visible across mixed-scale axes. Options: see **Feature summary** above. `VectorScale -> None` makes lengths proportional to magnitude; a real value scales relative to grid spacing.

**Examples:**

```mathematica
(* Simple rotation field *)
In[1]:= VectorPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[1]= -Graphics-

(* Gradient field with Rainbow colouring *)
In[2]:= VectorPlot[{x, y}, {x, -3, 3}, {y, -3, 3},
          ColorFunction -> "Rainbow", VectorPoints -> 20]
Out[2]= -Graphics-

(* Magnitude-proportional arrows *)
In[3]:= VectorPlot[{Sin[y], Cos[x]}, {x, 0, 2Pi}, {y, 0, 2Pi},
          VectorScale -> None]
Out[3]= -Graphics-

(* RegionFunction: unit disk only *)
In[4]:= VectorPlot[{-y, x}, {x, -1.5, 1.5}, {y, -1.5, 1.5},
          RegionFunction -> Function[{x,y}, x^2 + y^2 < 1]]
Out[4]= -Graphics-
```
