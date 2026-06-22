# Graphics

A 2D graphics engine modeled on Mathematica's: symbolic primitives compose
into a `Graphics[...]` object, which `Show[]` (or `Plot[]`, which
auto-displays) renders in an interactive window via a Raylib backend.
Renders a `-Graphics-` placeholder in the textual REPL, matching
Mathematica's own behavior outside a notebook front end.

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
centered at `{x,y}`.

## Disk
A graphics primitive: `Disk[{x,y}, r]` is a filled disk of radius `r`
centered at `{x,y}`.

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
Displays a `Graphics[...]` object.
- `Show[graphics]`: opens an interactive window rendering `graphics` and
  returns it. The window blocks the REPL until closed.
- `Show[graphics, opt -> val, ...]`: as above, merging the given options
  into `graphics`'s option list (a later/explicit option overrides one
  already present; the merged `Graphics[...]` is what's returned).

Rendering options (read by `Show`/`Plot`'s renderer, stored as trailing
`Rule[...]` arguments on the `Graphics[...]` object):

| Option | Default (bare `Graphics`) | Default (via `Plot`) | Meaning |
|---|---|---|---|
| `Axes` | `False` | `True` | draw coordinate axes with tick labels |
| `AspectRatio` | `Automatic` (fit window) | `1/GoldenRatio` | plot-box height/width ratio |
| `PlotRange` | `Automatic` (fit data bounding box) | same | `{{xmin,xmax},{ymin,ymax}}` |
| `PlotStyle` | a default blue | same | style directive(s) used as the initial draw color |
| `Background` | white | white | window background color |
| `ImageSize` | `{800, 600}` | same | window size in pixels, or a single width |
| `AxesLabel` | none | none | `{xlabel, ylabel}` |
| `PlotLabel` | none | none | title drawn above the plot |

**Interactive controls** (within the window): mouse wheel to zoom,
right-drag (or middle-drag) to pan, `Q`/`E` to rotate the view, `R` to
reset, `Esc` or the window's close button to return to the REPL.

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
- `Plot[f, {x, xmin, xmax}, opts...]`: as above, with sampling and/or
  rendering options (see below and the `Show` table above).
- `HoldAll`: `f` and the iterator spec are not pre-evaluated (`x` has no
  value until sampling binds it).

Sampling options (consumed by `Plot` itself; not stored on the resulting
`Graphics[...]`), implemented by one shared adaptive sampler so the same
algorithm and option semantics will back future plotting functions
(`ListPlot`, `ParametricPlot`, ...):

| Option | Default | Meaning |
|---|---|---|
| `PlotPoints` | `25` | initial uniform sample count across `[xmin,xmax]` |
| `MaxRecursion` | `6` | max bisection depth per initial interval when curvature/a singularity demands it |
| `MaxPlotPoints` | `Infinity` | overall cap on stored sample points |

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
```
