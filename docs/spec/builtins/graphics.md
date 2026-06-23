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

```mathematica
In[1]:= Graphics[{Point[{0, 0}]}]
Out[1]= -Graphics-

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
closed polygon through the given vertices.

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
| `Axes` | `False` | `True` | draw coordinate axes with tick labels |
| `AspectRatio` | `Automatic` (true geometry) | `1/GoldenRatio` | window height-to-width ratio: `Automatic` shapes it from the data extent, `Full` keeps the `ImageSize` box and stretches the data to fill it, or an explicit ratio `a` (incl. symbolic like `1/GoldenRatio`) |
| `PlotRange` | `Automatic` (data bbox, y spike-clipped) | same | `Automatic` (default), `All` (no clipping), `{{xmin,xmax},{ymin,ymax}}` (fix both axes), or `{ymin,ymax}` (fix y only, x stays automatic) |
| `PlotStyle` | a default blue | same | style directive(s) used as the initial draw color |
| `Background` | white | white | window background color |
| `ImageSize` | width `800` | same | `w` sets the width (height follows from `AspectRatio`), or `{w, h}` fixes both in pixels (then `AspectRatio` shapes the data inside the fixed box) |
| `AxesLabel` | none | none | `{xlabel, ylabel}` |
| `PlotLabel` | none | none | title drawn above the plot |

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

The adaptive sampler judges each candidate segment by the maximum **vertical**
gap between the curve and the straight chord at three interior probe points
(1/4, 1/2, 3/4), normalized by the displayed y-extent. Probing more than the
midpoint alone makes it robust to oscillatory functions (e.g. a sum of sines
of different frequencies), where a single-midpoint test would resonate with the
grid and miss the wiggle entirely no matter how high `MaxRecursion` is set.

**Features**:
- `HoldAll`, `Protected`.
- Declines to evaluate if the iterator spec isn't `{x, xmin, xmax}` with
  numeric (possibly symbolic-but-numeric, e.g. `2 Pi`) bounds, or any
  option has a malformed value (e.g. `PlotPoints -> 1`).
- Auto-displays (like Mathematica's notebook front end auto-rendering a
  `Graphics[]` result) -- no separate `Show[]` call is needed to see the
  plot, though `Show[Plot[...], opt -> val]` works for re-styling.

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
```
