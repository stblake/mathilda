# StreamPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Traces streamlines of the 2-D vector field {vx, vy} by RK4 integration from a grid of seed points, and returns a Graphics\[{Arrow\[...\], ...}, opts\] object (auto-displayed). StreamPlot is HoldAll: vx, vy, and the iterator specs are held unevaluated until x and y are given numeric values. Options: StreamPoints  - Integer n (n x n seed grid) or Automatic (default 15 x 15). StreamScale   – Automatic (8%% of domain diagonal, default), None (full run), or a real fraction of the domain diagonal. StreamStyle   – Style directive(s) applied to all streams. StreamColorFunction / ColorFunction – f\[x,y,vx,vy,speed\] (or fewer args) returning a color, or a named ramp: "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" (all keyed to scaled speed). RegionFunction – f\[x,y\] mask; seeds outside the region are skipped. PlotLegends   – Automatic / "Expressions" / explicit label list. StreamAnimate  – True | False (default). When True, each streamline is emitted as AnimatedStreamline\[...\] instead of Line\[...\]: the shape is drawn identically, but an interactive window (Show, or embedded in Animate/Manipulate) also draws particles flowing along it in real time. Standard Graphics options (PlotRange, Axes, AspectRatio, Frame, …) pass through to the Graphics\[...\] result.

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[1]= -Graphics-

In[2]:= StreamPlot[{1 - y^2, x}, {x, -3, 3}, {y, -2, 2}]
Out[2]= -Graphics-

In[3]:= StreamPlot[{Sin[x + y], Cos[x - y]}, {x, 0, 2 Pi}, {y, 0, 2 Pi}]
Out[3]= -Graphics-

In[4]:= Graphics[{Blue, Arrow[{{0,0}, {1,0}, {1,1}}]}]
Out[4]= -Graphics-
```

### Options (6)

```mathematica
In[5]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamPoints -> 25]
Out[5]= -Graphics-

In[6]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamScale -> None]
Out[6]= -Graphics-

In[7]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamColorFunction -> "Rainbow"]
Out[7]= -Graphics-

In[8]:= StreamPlot[{x, -y}, {x, -2, 2}, {y, -2, 2}, StreamStyle -> {Thickness[0.004], RGBColor[0.8, 0.2, 0.1]}]
Out[8]= -Graphics-

In[9]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 1.5^2]]
Out[9]= -Graphics-

In[10]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamAnimate -> True]
Out[10]= -Graphics-
```

### Applications (5)

```mathematica
In[11]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[11]= -Graphics-

In[12]:= StreamPlot[{1 - y^2, x}, {x, -3, 3}, {y, -2, 2}]
Out[12]= -Graphics-

In[13]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamPoints -> 30]
Out[13]= -Graphics-

In[14]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamAnimate -> True]
Out[14]= -Graphics-

In[15]:= Attributes[StreamPlot]
Out[15]= {HoldAll, Protected}
```

## Algorithm

streamplot.c — StreamPlot[{vx,vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]

Traces streamlines of a 2-D vector field by RK4 integration from a grid

```text
of seed points, emitting one Arrow[...] primitive per stream.  The result
```

is returned as a Graphics[...] object (auto-displayed by the REPL).

## Implementation notes

**Algorithm.** `builtin_streamplot` is `HoldAll`. It traces **evenly-spaced**
streamlines of the field `{vx, vy}` by the Jobard–Lefebvre scheme. Each line is
grown in **both directions** from a seed (`grow_streamline`/`integrate_dir`) by
`rk4_step_unit`, a 4th-order Runge–Kutta step of the **normalised (unit) field**
`û = v/|v|` with signed step `s = ±h` — integrating the unit field gives uniform
arc-length point spacing independent of local speed. With `ext = min(dx, dy)` and
`np = StreamPoints` (default **25**) the data-space parameters are: target line
separation `d_sep = ext/np`; RK4 step `h = d_sep/8`; proximity-stop radius²
`d_test2 = 0.25*d_sep^2` (a line halts ~`0.5*d_sep` from an existing line);
closed-orbit return radius² `loop2 = (1.6 h)^2` after a minimum arc `6 h`. Even
spacing is enforced by a uniform spatial hash `SGrid` with **cell size exactly
`d_sep`**: every accepted streamline point is registered, and
`sgrid_within2(x, y, r2)` searches the query cell's 3×3 neighbourhood (correct
because `r2 <= cell^2`), giving `O(1)` proximity tests. Seeds are drawn from a
candidate grid twice as fine as `d_sep` (`2*np` per axis), sorted by distance from
the domain centre so placement grows outward. A line stops at the domain
boundary, a `RegionFunction` rejection, a critical point (`|v| <= crit`),
proximity to an existing line, a closed-orbit return, or the step cap
(`max_arc = StreamScale*diag` or `4*diag`, `max_steps` capped at 20000). Default
rendering (`emit_dashed_stream`) walks arc length alternating drawn `dash_len`
and blank `gap_len`, each dash a curve-following `Line` capped by a filled
`Polygon` arrowhead triangle (`chevron_tip`) — a bare triangle, not `Arrow`, to
avoid a straight chord shaft. `StreamAnimate -> True` instead emits one
`AnimatedStreamline[...]` per line plus periodic `Arrow` chevrons, whose particle
dots flow in the interactive renderer. Colour priority is
`StreamColorFunction`/`ColorFunction` at the line midpoint, else a Viridis
speed ramp, else a single `StreamStyle` directive. The result is an inert
`Graphics[prims, opts...]`.

**Data structures.** `SGrid` (the `d_sep`-cell spatial hash: `nx*ny` buckets of
`Point2` with per-cell `cnt`/`cap`); `SeedCand {x, y, d2}` sorted candidate array;
`PBuf` growing per-half point buffers; `all_streams`/`all_lengths`/`all_speed`
parallel arrays. Integration runs in data space; finished polylines map to world
space only for rendering.

**Complexity / limits.** `(2*np)^2` seed candidates (default 2500), `O(1)`
proximity per point via the hash, per-line steps capped at 20000. Returns `NULL`
if no streamline is produced. Defaults `Frame -> True`, `Axes -> False`,
`AspectRatio -> 1`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Show](../../graphics/Show/), [Animate](../../graphics/Animate/), [Manipulate](../../graphics/Manipulate/)

- B. Jobard and W. Lefebvre, *Creating Evenly-Spaced Streamlines of Arbitrary Density*, Proc. 8th Eurographics Workshop on Visualization in Scientific Computing (1997) 43-55.
- Source: [`src/graphics/streamplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/streamplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_streamplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streamplot.c)

## Notes & additional examples

### Notes

`StreamPlot` is `HoldAll`. It traces **evenly-spaced** streamlines by the
Jobard–Lefebvre scheme: each line is grown in both directions from a seed by RK4
integration of the *normalised* field (a fixed arc-length step, so point spacing and
rendered curvature stay uniform regardless of local speed), and stops when it nears
an existing line, leaves the domain, reaches a critical point, or closes on itself.
A uniform spatial hash whose cell size equals the target line separation keeps the
lines evenly spaced with `O(1)` proximity tests, and seeds are placed outward from
the domain centre.

Each line renders in Mathematica's dashed-arrow style — short curve-following `Line`
dashes capped by filled `Polygon` arrowheads, so the flow direction reads
everywhere. `StreamPoints` sets the density (default 25); `StreamScale` caps a line's
arc length, while the default lets every line run to its natural end.
`StreamAnimate -> True` instead emits one `AnimatedStreamline[...]` per line whose
particle dots flow downstream in the interactive renderer.
