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

Sampling options (consumed by `Plot` itself; not stored on the resulting
`Graphics[...]`), implemented by one shared adaptive sampler so the same
algorithm and option semantics will back future plotting functions
(`ListPlot`, `ParametricPlot`, ...):

| Option | Default | Meaning |
|---|---|---|
| `PlotPoints` | `50` | initial uniform sample count across `[xmin,xmax]` |
| `MaxRecursion` | `6` | max bisection depth per initial interval when curvature/a singularity demands it |
| `MaxPlotPoints` | `Infinity` | overall cap on stored sample points |
| `Mesh` | `None` | `All` overlays a dot at every evaluation point (in the curve's colour); `None` draws the line only. `True`/`False` are accepted as synonyms |
| `RegionFunction` | none | `f[x,y]` (or `f[x]`) -- a sample where this isn't `True` is dropped exactly like a singularity (the curve gaps there) |
| `Exclusions` | `None` | `{x1, x2, ...}` or `{x == a, ...}`: forces a discontinuity at each given x, independent of curvature/singularity detection |
| `ColorFunction` | none | a function `f[xscaled,yscaled]`/`f[xscaled]` returning a color literal, applied per sampled segment; or a named ramp string: `"Rainbow"`, `"Temperature"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"` (see **Named ColorFunction ramps** section). Overrides `PlotStyle` for the curve itself |
| `ColorFunctionScaling` | `True` | whether `ColorFunction` receives `x`/`y` scaled to `[0,1]` over the plot's range, or raw data values |
| `Filling` | `None` | `Axis` (to `y=0`), `Bottom`/`Top` (to the run's own min/max y), or a number (to that y) -- a strip of small quad `Polygon[]`s (one per consecutive sampled-point pair) hugging the curve down to the baseline, under the curve's outline |
| `FillingStyle` | none | color for the fill; default is the curve's own color at `Opacity[0.3]` |
| `PlotLegends` | none | `Automatic`/`"Expressions"` (label each curve with its own expression, via `ToString`) or an explicit list of label strings -- drawn as a small swatch+label box, top-right below the toolbar |

The adaptive sampler judges each candidate segment by the maximum **vertical**
gap between the curve and the straight chord at three interior probe points
(1/4, 1/2, 3/4), normalized by the displayed y-extent. Probing more than the
midpoint alone makes it robust to oscillatory functions (e.g. a sum of sines
of different frequencies), where a single-midpoint test would resonate with the
grid and miss the wiggle entirely no matter how high `MaxRecursion` is set.

Every probe is first **clamped to the displayed y-band** — an explicit numeric
`PlotRange` y, or (during zoom re-sampling) the current visible extent. This
keeps refinement on the part of the curve the window actually shows: a function
that dives far out of frame (a truncated Taylor series past its radius of
convergence, a steep asymptote) would otherwise pour its whole detail budget
into the off-screen plunge and leave the visible body and its on-screen crossing
coarsely faceted. Clamped, the off-band stretch collapses onto the clip line and
reads as flat, so the recursion concentrates where it is seen. With no explicit
`PlotRange` y the curve is sampled over its full extent, exactly as before.

This vertical-gap test is the **primary** refinement driver, and because a
chord's sagitta grows with curvature, sample density tracks curvature: points
cluster where the curve bends (the peaks of `Sin`, the ends of `x^3`, the body
of a spike) and thin out where it is straight — including steep but locally
linear stretches such as `Sin`'s zero-crossings. A secondary, deliberately
loose chord-**length** cap (≈1/12 of the frame) acts only as a backstop, so a
near-vertical asymptotic approach (e.g. `Log[x]` near `0`) cannot leave a
conspicuous on-screen gap; it is kept generous on purpose, since a tight cap
keys on slope rather than curvature and would invert the density, over-sampling
straight-steep regions at the expense of the curvy ones.

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if the iterator spec isn't `{x, xmin, xmax}` with
  numeric (possibly symbolic-but-numeric, e.g. `2 Pi`) bounds, or any
  option has a malformed value (e.g. `PlotPoints -> 1`).
- Auto-displays (like Mathematica's notebook front end auto-rendering a
  `Graphics[]` result) -- no separate `Show[]` call is needed to see the
  plot, though `Show[Plot[...], opt -> val]` works for re-styling.

**Not yet implemented**: `ColorFunction`'s companion `MeshFunctions`/
`MeshShading` (finer control over an already-working `Mesh`), curve-to-curve
`Filling -> {1 -> {2}}` for multi-curve plots, `Method` (one sampler, not
multiple selectable algorithms), `ClippingStyle`, `PlotTheme`, and
`BaseStyle` (broader than `LabelStyle`, overlapping confusingly with
`PlotStyle`) are silently accepted on the option list (`Plot` copies any
option it doesn't specifically recognize through onto the returned
`Graphics[...]` rather than erroring) but have no rendering effect yet.
Likewise, notebook-embedding-only options (`BaselinePosition`,
`AlignmentPoint`, `CoordinatesToolOptions`, `ImageSizeRaw`, `PlotRegion`,
`ImageMargins`, `ImagePadding`) and `WorkingPrecision` (the sampler is
architecturally machine-precision `double` throughout) are accepted but
inherently inert in a standalone window.

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

`Plot3D` reuses `Plot`'s option *semantics* wherever a 3D analogue makes
sense, sharing the actual evaluation code (`src/graphics/plot_common.c`) for
`RegionFunction`, `ColorFunction`, and the multi-curve/multi-surface palette:

| Option | Default | Meaning |
|---|---|---|
| `PlotPoints` | `25` | initial per-axis grid resolution (an `n x n` grid) |
| `MaxRecursion` | `2` | doubles the *whole* grid's resolution (up to this many times, capped at 200 points/axis) while a cell-center-vs-bilinear-interpolant flatness check fails -- a global, crack-free analogue of `Plot`'s per-interval adaptive bisection (a per-cell quadtree would leave T-junction cracks where differently-refined cells meet) |
| `Mesh` | `True` | overlays the grid wireframe on the surface (unlike `Plot`'s default `None` -- Mathematica's `Plot3D` shows mesh lines out of the box too); `None`/`False` draws the filled surface only. **Only interior grid lines are drawn** -- the perimeter edges of the domain are omitted, giving a smooth boundary silhouette |
| `PlotStyle` | `RGBColor[0.2,0.4,0.8]` | the surface's fill color. For a single surface, a direct color literal; for multi-surface plots (`{f1,f2,...}`), a `List` of colors is indexed per surface (cycling if shorter than the surface count), or a single literal applies to all. Ignored per-cell where `ColorFunction` overrides |
| `ColorFunction` | none | a function returning a color literal, evaluated per grid cell. **Calling convention** (tried in order): `f[xs,ys,zs]`, then `f[xs,zs]`, then `f[zs]`. Named ramp strings: `"Rainbow"` (inverted height hue sweep), `"Temperature"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"` (see **Named ColorFunction ramps**) |
| `ColorFunctionScaling` | `True` | when `True`, the arguments passed to `ColorFunction` are scaled to `[0,1]` within the data range; `False` passes raw coordinate values |
| `RegionFunction` | none | `f[x,y,z]` (Mathematica's `Plot3D` convention) tried first; falls back to `Plot`'s `f[x,y]`/`f[x]` forms if that doesn't resolve to `True`/`False`, so a `RegionFunction` written for `Plot` keeps working |
| `ExclusionStyle` | `GrayLevel[0.35]` | style directive (a color literal or `Thickness[...]` etc.) used to draw the boundary edges between included and excluded grid cells when `RegionFunction` is active |
| `PlotRange` | `Automatic` | an explicit `{zmin,zmax}` (or the last `{zmin,zmax}` of a longer nested list) feeds the same flatness check `MaxRecursion` uses, and bounds the rendered box |
| `Axes` | `True` | draws a wireframe bounding box with tick labels |
| `Lighting` | `Automatic` | surface shading model. `Automatic` (default): per-face Lambertian flat shading — each polygon's face normal is dotted with a directional light that is **fixed relative to the camera** (upper-right-front in view space: `0.3 right + 1.0 up + 0.5 forward`, normalised), so the lit/dark pattern updates as you orbit. Intensity `I = 0.3 + 0.7 × |n·L|` (ambient 0.3 + diffuse 0.7, two-sided so back faces are lit too). `None` or `False`: disables shading and draws surfaces in their raw `PlotStyle`/`ColorFunction` color |

Other options (`PlotLabel`, `Background`, `ImageSize`, `AxesLabel`, ...) are
evaluated once (since `Plot3D` is `HoldAll`) and copied through onto the
`Graphics3D[...]` result, exactly like `Plot`'s fallback branch. 2D-only
chrome with no 3D analogue in this engine (`Filling`, `FillingStyle`,
`Exclusions`, `Frame*`, `GridLines*`, `PlotLegends`) is likewise not
specially recognized -- it passes through inertly rather than erroring, so
copying a `Plot[]` option list onto `Plot3D[]` by habit doesn't break
anything, it just has no effect for those names.

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if either iterator spec isn't `{var, min, max}` with
  numeric (possibly symbolic-but-numeric) bounds, if `PlotPoints -> 1` (or
  any other malformed sampling option), or if no grid cell is valid (e.g. a
  `RegionFunction` that excludes every sampled point).
- Auto-displays, exactly like `Graphics`/`Plot`.

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

`HoldAll`: the body and all iterator specs are not pre-evaluated.

Adaptive sampling (1-iterator form) uses a 2D chord-deviation test in `(x, y)`
space: the midpoint of each candidate segment is checked for Euclidean
displacement from the linear interpolant of the two endpoints, normalised by the
bounding-box diagonal.  Three interior probes per interval (quarter-point and
midpoint) avoid aliasing against periodic curves.

| Option | Default | Meaning |
|---|---|---|
| `PlotPoints` | `25` | initial uniform sample count (1-iter) or grid size (2-iter) |
| `MaxRecursion` | `6` | maximum adaptive-subdivision depth (1-iter only) |
| `MaxPlotPoints` | `Infinity` | overall point cap (1-iter only) |
| `Mesh` | `None` | `All` overlays sample dots (1-iter) or grid lines (2-iter) |
| `PlotLegends` | none | `Automatic` or `"Expressions"` labels each curve with its body expression; an explicit `{l1, l2, ...}` uses those labels. Embeds `$PlotLegendData` in the Graphics result for the renderer to draw a legend box. |
| `ColorFunction` | none | named ramp string (`"Rainbow"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"`, `"Temperature"`) or a function `f[t_scaled]` (1-iter) / `f[t_scaled, r_scaled]` (2-iter). `t` is scaled to `[0,1]` when `ColorFunctionScaling -> True` |
| `ColorFunctionScaling` | `True` | whether to scale parameters before `ColorFunction` |
| `RegionFunction` | none | `f[x,y]` mask; points where it returns `False` are excluded |
| `PlotStyle` | `RGBColor[0.2,0.4,0.8]` | curve/polygon color; a `List` of directives is accepted, so `PlotStyle -> {Blue, Opacity[0.4]}` gives a semi-transparent fill for the two-iterator form |
| `AspectRatio` | `1` | square by default (both axes equally important) |
| `Axes` | `True` | draws coordinate axes |

All other display options (`PlotRange`, `PlotRangePadding`, `AxesLabel`,
`AxesOrigin`, `AxesStyle`, `Frame`, `FrameLabel`, `FrameStyle`, `FrameTicks`,
`GridLines`, `GridLinesStyle`, `Prolog`, `Epilog`, `PlotLabel`, `Background`,
`ImageSize`, `TicksStyle`, `LabelStyle`, `RotateLabel`) are evaluated and
passed through to the `Graphics[...]` result, where the renderer interprets
them identically to how it does for bare `Graphics[]` objects.

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or `PlotPoints < 2`.
- Auto-displays exactly like `Graphics`/`Plot`.

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

`PolarPlot` is `HoldAll`: `r` and the iterator spec are held unevaluated
until sampling begins. Negative `r` values are allowed (they plot in the
opposite direction, matching standard polar convention). Default
`AspectRatio -> 1` (equal axes, so circles look round).

**Options** (same as `ParametricPlot` — see that section for details):

| Option | Default | Description |
|---|---|---|
| `PlotPoints` | `75` | Initial sample count per curve |
| `MaxRecursion` | `6` | Adaptive refinement depth |
| `MaxPlotPoints` | `Infinity` | Total point cap |
| `Mesh` | `None` | `All` overlays evaluation dots |
| `ColorFunction` | `None` | named ramp string or `f[t]`; `t` is scaled theta ∈ [0,1] (see **Named ColorFunction ramps**) |
| `ColorFunctionScaling` | `True` | Normalise theta to `[0,1]` |
| `RegionFunction` | `None` | `f[x,y]` mask |
| `PlotStyle` | palette | Color/style directives |
| `PlotLegends` | `None` | `Automatic` / `"Expressions"` / label list |
| `PolarAxes` | `False` | Accepted; polar grid overlay not yet rendered |
| `AspectRatio` | `1` | Equal axes (unlike `Plot`'s 1/GoldenRatio) |
| `Axes`, `Frame`, `PlotRange`, `GridLines`, `PlotLabel`, … | — | Pass through to `Graphics[...]` |

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

`HoldAll`: the body and all iterator specs are not pre-evaluated.

Adaptive sampling (1-iterator form) uses a 3D Euclidean chord-deviation test
in `(x, y, z)` space, normalized by the curve's bounding-box diagonal. Three
interior probes per interval prevent aliasing against periodic curves, exactly
mirroring `ParametricPlot`'s 2D algorithm extended to 3D.

`ColorFunction` receives **scaled spatial coordinates** `{xs, ys, zs}` (in
`[0,1]` over the sampled range when `ColorFunctionScaling -> True`), not
parameter values — this matches `Plot3D`'s convention and means a function
like `Function[{x,y,z}, Hue[z]]` colors by height. `"Rainbow"` maps hue to
the z-extent for surfaces and to z-range for curves. `RegionFunction` is
tried as `f[x,y,z]` first, then falls back to `f[x,y]` forms.

| Option | Default | Meaning |
|---|---|---|
| `PlotPoints` | `25` | initial uniform sample count (1-iter) or grid size (2-iter) |
| `MaxRecursion` | `6` | maximum adaptive-subdivision depth (1-iter only) |
| `MaxPlotPoints` | `Infinity` | overall point cap (1-iter only) |
| `Mesh` | `None` | `All`/`True` overlays sample dots (1-iter) or grid lines (2-iter) |
| `PlotLegends` | none | `Automatic` or `"Expressions"` labels each curve; an explicit `{l1,...}` uses those labels |
| `ColorFunction` | none | named ramp string or `f[x,y,z]`/`f[x,z]`/`f[z]` (scaled spatial); colors per segment (1-iter) or per cell (2-iter). Ramps: `"Rainbow"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"`, `"Temperature"` |
| `ColorFunctionScaling` | `True` | whether spatial coords are scaled to `[0,1]` before `ColorFunction` |
| `RegionFunction` | none | `f[x,y,z]` mask; points where it returns `False` are excluded |
| `PlotStyle` | `RGBColor[0.2,0.4,0.8]` | curve/surface color; a `List` of directives is accepted, so `PlotStyle -> {Blue, Opacity[0.4]}` gives a semi-transparent surface for the two-iterator form |
| `Axes` | `True` | draws the 3D bounding-box wireframe with tick labels |
| `Lighting` | `Automatic` | `Automatic`: per-face Lambertian shading (same model as `Plot3D`); `None`/`False`: raw flat color |

All other options (`PlotRange`, `AxesLabel`, `PlotLabel`, `Background`,
`ImageSize`, ...) are evaluated and passed through to the `Graphics3D[...]`
result. `AspectRatio` is silently ignored (the orbit camera has no fixed
aspect ratio). Interactive controls are the same as `Plot3D` (drag to rotate,
scroll to zoom, `R` to reset, `S` for screenshot, `Esc` to close).

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or `PlotPoints < 2`.
- Auto-displays exactly like `Graphics3D`/`Plot3D`.

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

`ListPlot`-specific options (the rest pass through to the `Graphics[...]`
result and are interpreted by the renderer — see the `Show` table above):

| Option | Default | Meaning |
|---|---|---|
| `Joined` | `False` | `True` connects the points with a `Line` polyline instead of drawing point markers |
| `DataRange` | `Automatic` | the x-range `{xmin, xmax}` to assume for a heights list (`Automatic`/`All` = `1` to `n`); `DataRange -> All` forces a flat list of pairs to be read as several height datasets |
| `Filling` | `None` | `Axis` / `Bottom` / `Top` / a number — fills down to that baseline. With `Joined -> True` the whole region under the connecting curve is filled continuously (quad/triangle `Polygon[]` strip, shared with `Plot`); otherwise a vertical stem is drawn from each point |
| `FillingStyle` | `Automatic` | colour/opacity directive for the fill (default `Opacity[0.3]`) |
| `PlotStyle` | `Automatic` | colour for the points/line (default blue `RGBColor[0.2, 0.4, 0.8]`) |
| `PlotLegends` | `None` | `{lbl1, ...}` labels each dataset; `Automatic` labels by index; emits the internal `$PlotLegendData` the renderer's legend box reads |
| `PlotMarkers` | `None` | accepted; markers currently render with the default glyph |

Plot-style defaults distinct from a bare `Graphics[]`: `Axes -> True`,
`AspectRatio -> 1/GoldenRatio`. A supplied `Frame -> True` suppresses the
default `Axes -> True`, as for `Plot`.

Marker size is **adaptive**: when no explicit `PointSize` is given, the
renderer shrinks the dots as the scatter grows denser, so a large point cloud
stays legible instead of merging into an ink blob. The home-zoom radius is held
to a small fraction of the mean inter-point spacing (≈ `sqrt(area / N)` for `N`
points over the plot region), capped at the sparse default and floored at a
sub-pixel minimum. It therefore depends on the point count, the window size, and
the pixel resolution. An explicit `PointSize` overrides this entirely.

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

`StreamPlot` is `HoldAll`: the field components `vx`, `vy` and the
iterator specs are held unevaluated — `x` and `y` have no values until
the sampler substitutes numeric coordinates, exactly like `Plot`'s
function body.

```mathematica
StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}]
StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

**Options**

| Option | Default | Meaning |
|---|---|---|
| `StreamPoints` | `Automatic` (15×15 grid) | Integer `n` for an n×n seed grid; `Automatic` uses the default 15×15 |
| `StreamScale` | `Automatic` | `Automatic` limits each stream to ≈ 8% of the domain diagonal; `None` lets streams run until they leave the domain; a positive real sets the fraction of the diagonal |
| `StreamStyle` | *(thin steel-blue)* | A style directive or list applied to every stream (e.g. `Thickness[0.003]`, `RGBColor[...]`) |
| `StreamColorFunction` | `None` | `f[x,y,vx,vy,speed]` (or fewer args) → color directive per stream at midpoint; named ramp strings: `"Rainbow"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"`, `"Temperature"` (keyed to scaled speed) |
| `ColorFunction` | `None` | Alias for `StreamColorFunction` |
| `RegionFunction` | `None` | `f[x,y]` mask; seeds and integration steps outside the region are skipped |
| `PlotLegends` | `None` | `Automatic` / `"Expressions"` / explicit list |
| `PlotPoints` | — | Alias for `StreamPoints` (integer only) |

All other options (`PlotRange`, `Axes`, `AspectRatio`, `Frame`, `AxesLabel`,
`GridLines`, `PlotLabel`, `Background`, `ImageSize`, …) pass through to the
`Graphics[...]` result.

Default style: `Axes -> True`, `AspectRatio -> 1` (square domain).

**Arrow primitive** — `Arrow[{{x1,y1}, ..., {xn,yn}}]` draws a directed
polyline with a filled arrowhead triangle at the final point. The
arrowhead size scales with line thickness and viewport size so it remains
visible at any zoom level. `Arrow` is an inert protected head (like `Line`
or `Polygon`) and can be used directly inside `Graphics[...]`.

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

**Features**:
- `HoldAll`, `Protected`.
- RK4 integration; step size adapts to seed density and domain size.
- Declines to evaluate if the field arg is not a 2-element List, or if
  bounds aren't numeric.
- `StreamColorFunction` is evaluated at each stream's midpoint; tries
  `f[x,y,vx,vy,speed]` → `f[x,y,vx,vy]` → `f[x,y]` → `f[speed]` in
  order, using the first form that returns a recognized color.
- Arrow arrowhead size is viewport-relative: it stays visible regardless
  of `PlotRange` scale or interactive zoom.

## ContourPlot
Generates iso-contour lines of a 2-D function `f(x, y)` using the marching
squares algorithm and returns a `Graphics[...]` object (auto-displayed).

```mathematica
ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}]
ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
```

`ContourPlot` is `HoldAll`: `f` is held unevaluated until `x` and `y` are
given numeric values, exactly like `Plot`'s function body.

**Algorithm**: evaluates `f` on a `(PlotPoints+1) × (PlotPoints+1)` grid,
then for each contour level runs marching squares over the `PlotPoints ×
PlotPoints` grid of cells, linearly interpolating the exact crossing points on
cell edges. Saddle cells (where both pairs of opposite corners straddle the
level) use the bilinear cell-centre value to choose the correct diagonal. If
`ContourShading` is active, each grid cell is filled with a coloured
`Rectangle[]` primitive whose colour corresponds to the cell's average `z`
value via `ColorFunction` or the built-in blue-cyan-yellow-red thermal ramp.

**Options**

| Option | Default | Meaning |
|---|---|---|
| `Contours` | `10` | Integer `n`: `n` evenly spaced auto levels in `[zmin, zmax]`; or `{c1, c2, ...}` to set levels explicitly |
| `ContourStyle` | `Automatic` | Style directive(s) for the contour lines. `Automatic`: each level is coloured by its height using the thermal ramp. A single directive is applied to all levels; a `List` of directives cycles through the levels. `None`/`False` suppresses lines entirely |
| `ContourLabels` | `False` | `True`: draws the level's `z` value as a `Text[]` at the midpoint of each level's first visible grid-segment |
| `ContourShading` | `Automatic` | `True`: fill each grid cell with a colour derived from its average `z`. `False`/`None`: lines only. `Automatic`: enable shading when `ColorFunction` is set, otherwise lines only |
| `ColorFunction` | `None` | named ramp string (`"Rainbow"`, `"Temperature"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"`) or a function `f[t]` with `t ∈ [0,1]` after scaling. Applied to both cell shading and auto contour-line colours |
| `ColorFunctionScaling` | `True` | `True`: normalise `z` to `[0,1]` before calling `ColorFunction`; `False`: pass raw `z` |
| `PlotPoints` | `25` | Grid resolution per axis; increase for smoother contours |
| `RegionFunction` | `None` | `f[x,y]` mask; cells whose centre lies outside return `False` are neither shaded nor contoured |
| `Axes` | `True` | Coordinate axes (same Plot/StreamPlot default) |
| `AspectRatio` | `1` | Default square aspect (both axes equally important) |

All other `Graphics` options (`Frame`, `PlotRange`, `AxesLabel`, `GridLines`,
`PlotLabel`, `Background`, `ImageSize`, `Prolog`, `Epilog`, …) pass through
to the `Graphics[...]` result.

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or the function argument
  does not have exactly two iterator specs.
- Marching squares handles all 16 cell states, including the two saddle
  cases (5 and 10), with bilinear centre disambiguation.
- An explicit `PlotRange` passed in options suppresses the auto-range
  injection (which otherwise pins the range to `{xmin,xmax}` × `{ymin,ymax}`
  so the axes match the sampling domain).

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

Renders `f(x,y)` as a heatmap by shading each grid cell with the colour
corresponding to the cell's average `f` value.  Returns a `Graphics[...]`
object (auto-displayed).  `DensityPlot` is `HoldAll`: `f` is held
unevaluated until `x` and `y` are bound to numeric values.

**Default ColourFunction:** thermal blue→yellow ramp (same as `StreamPlot`).

| Option | Default | Meaning |
|--------|---------|---------|
| `PlotPoints` | `50` | Grid resolution per axis |
| `ColorFunction` | thermal | named ramp string or `f[t]→color` (`t∈[0,1]`). Ramps: `"Rainbow"`, `"Temperature"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"` |
| `ColorFunctionScaling` | `True` | Normalise `z` to `[0,1]` before calling `ColorFunction` |
| `RegionFunction` | none | `f[x,y]` mask; excluded cells are not drawn |
| `PlotLegends` | none | `Automatic`: attach a vertical color-scale bar |
| Standard `Graphics` options | | `Axes→True`, `AspectRatio→1`, `Frame`, `PlotRange`, … |

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

Draws a vertical bar chart.  Single-dataset form: `n` bars at `x = 1..n`
with heights `v1..vn`.  Multi-dataset form: grouped bars — each inner list
is one group, each group sub-bar uses a distinct palette colour.
Returns a `Graphics[...]` object (auto-displayed).

| Option | Default | Meaning |
|--------|---------|---------|
| `ChartStyle` | palette | Color or `{color,...}` cycling through bars |
| `ChartLabels` | none | `{label,...}` drawn below each bar on the x-axis |
| `BarSpacing` | `0.2` | Gap as a fraction of bar width; `0` = touching |
| Standard `Graphics` options | | `Axes→True`, `AspectRatio→0.618`, `Frame`, `PlotRange`, … |

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

Bins numeric `data` and draws a frequency histogram.  The first form uses
Sturges' rule (`⌈log₂(n)⌉ + 1` bins) by default.
Returns a `Graphics[...]` object (auto-displayed).

| Option | Default | Meaning |
|--------|---------|---------|
| `ChartStyle` | palette | Color or `{color,...}` cycling through bins |
| `BarSpacing` | `0.2` | Gap as a fraction of bin width |
| Standard `Graphics` options | | `Axes→True`, `AspectRatio→0.618`, `Frame`, `PlotRange`, … |

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

Draws a grid of arrows showing the direction (and optionally magnitude) of
the vector field `{vx, vy}` at each grid point.  Returns a `Graphics[...]`
object (auto-displayed).  `VectorPlot` is `HoldAll`: `vx`, `vy` are held
unevaluated until `x` and `y` are bound to numeric values.

**Default colour:** thermal ramp keyed to local field magnitude (speed).

| Option | Default | Meaning |
|--------|---------|---------|
| `VectorPoints` | `15` | Grid density per axis (`n` → `n×n` grid); `Automatic` = 15 |
| `VectorScale` | `Automatic` | `Automatic`: equal-length arrows (direction only); `None`: proportional to magnitude; real `f`: arrow length = `f × grid_spacing` |
| `VectorStyle` | — | Style directive(s) applied globally to all arrows |
| `ColorFunction` | thermal | named ramp string (keyed to speed) or `f[vx,vy,speed]`/`f[speed]` → color. Ramps: `"Rainbow"`, `"CoolTones"`, `"WarmTones"`, `"Greyscale"`, `"Temperature"` |
| `ColorFunctionScaling` | `True` | Normalise speed to `[0,1]` before calling `ColorFunction` |
| `RegionFunction` | none | `f[x,y]` mask: skip grid points outside the region |
| Standard `Graphics` options | | `Axes→True`, `AspectRatio→1`, `Frame`, `PlotRange`, … |

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
