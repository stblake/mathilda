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
`ToString[]` would) as text centered at `{x,y}`, drawn with a built-in
single-stroke vector font.

## RGBColor / GrayLevel / Opacity / Thickness / PointSize
Style directives. Placed alongside primitives in a `Graphics[]` primitive
list, each affects every primitive that follows it (left to right), exactly
like Mathematica's directive semantics.
- `RGBColor[r, g, b]` / `RGBColor[r, g, b, a]`: color, components in `[0,1]`.
- `GrayLevel[g]` / `GrayLevel[g, a]`: a shade of gray.
- `Opacity[a]`: opacity in `[0,1]`.
- `Thickness[t]`: line thickness (plot coordinates) for `Line`/`Circle`.
- `PointSize[s]`: point radius (plot coordinates) for `Point`.

```mathematica
In[1]:= Graphics[{RGBColor[1, 0, 0], Point[{0, 0}], Line[{{0,0},{1,1}}]}]
Out[1]= -Graphics-
```

## Named color constants
`Red`, `Green`, `Blue`, `Black`, `White`, `Gray`, `Cyan`, `Magenta`,
`Yellow`, `Brown`, `Orange`, `Pink`, `Purple` and their light variants
`LightRed`, `LightGreen`, `LightBlue`, `LightGray`, `LightCyan`,
`LightMagenta`, `LightYellow`, `LightBrown`, `LightOrange`, `LightPink`,
`LightPurple` are ordinary protected `OwnValue`s evaluating to an
`RGBColor[...]`/`GrayLevel[...]` literal (e.g. `Red` is exactly
`RGBColor[1, 0, 0]`, `Black` is `GrayLevel[0]`) -- usable anywhere a
color literal is, including `PlotStyle`, `Background`, `FrameStyle`,
`AxesStyle`, `TicksStyle`, `GridLinesStyle`, and as a style directive
directly inside a primitive list (`Graphics[{Blue, Point[{0,0}]}]`).

Because `Plot`/`Show` are `HoldAll`/inert respectively, a non-`Held`
`Graphics[...]`'s arguments evaluate normally (so `Blue` resolves the
instant it's parsed), while `Plot` evaluates each of its own option values
once before storing them, specifically so a named color used inside e.g.
`Epilog -> {Red, ...}` resolves the same way.

```mathematica
In[1]:= Red
Out[1]= RGBColor[1, 0, 0]

In[2]:= Plot[Sin[x], {x, 0, 2 Pi}, PlotStyle -> Green, Background -> LightGray]
Out[2]= -Graphics-
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
| `ColorFunction` | none | a function `f[xscaled,yscaled]`/`f[xscaled]` returning a color literal (`RGBColor`/`GrayLevel`/`Hue`), applied per sampled segment -- or the string `"Rainbow"` for a built-in `Hue` sweep over x. Overrides `PlotStyle` for the curve itself |
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
