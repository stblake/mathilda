# VectorPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`VectorPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Draws a grid of arrows showing the direction (and optionally magnitude) of the vector field {vx, vy} at each grid point. VectorPlot is HoldAll: vx, vy are held unevaluated until x and y are bound to numeric values. Returns a Graphics\[...\] object. Options: VectorPoints   integer n → n×n grid (default 15); Automatic = 15 VectorScale    Automatic: equal-length arrows (direction only) None: proportional to magnitude real f: arrow length = f × grid spacing VectorStyle    style directive(s) applied to all arrows ColorFunction  named ramp string (keyed to speed) or f\[vx,vy,speed\]/f\[speed\]→color. Ramps: "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" ColorFunctionScaling  True (default): normalise speed to \[0,1\] RegionFunction f\[x,y\] mask: skip grid points outside the region Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= VectorPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[1]= -Graphics-
```

### Options (4)

```mathematica
In[2]:= VectorPlot[{x, y}, {x, -3, 3}, {y, -3, 3}, ColorFunction -> "Rainbow", VectorPoints -> 20]
Out[2]= -Graphics-

In[3]:= VectorPlot[{Sin[y], Cos[x]}, {x, 0, 2Pi}, {y, 0, 2Pi}, VectorScale -> None]
Out[3]= -Graphics-

In[4]:= VectorPlot[{-y, x}, {x, -1.5, 1.5}, {y, -1.5, 1.5}, RegionFunction -> Function[{x,y}, x^2 + y^2 < 1]]
Out[4]= -Graphics-

In[5]:= Head[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 5]]
Out[5]= Graphics
```

### Applications (5)

```mathematica
In[6]:= VectorPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
Out[6]= -Graphics-

In[7]:= VectorPlot[{x, y}, {x, -2, 2}, {y, -2, 2}, VectorScale -> None]
Out[7]= -Graphics-

In[8]:= VectorPlot[{-y, x}, {x, -1.5, 1.5}, {y, -1.5, 1.5}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 1]]
Out[8]= -Graphics-

In[9]:= Length[Cases[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 6], _Arrow, Infinity]]
Out[9]= 36

In[10]:= Attributes[VectorPlot]
Out[10]= {HoldAll, Protected}
```

## Algorithm

vectorplot.c — VectorPlot[{vx, vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]

Draws a grid of arrows showing the direction (and optionally magnitude) of the 2-D vector field {vx, vy}. Returns a Graphics[...] object auto- displayed by the REPL.

VectorPlot is HoldAll: vx, vy, and the iterator specs are held unevaluated until x and y are bound to numeric values (same semantics as ContourPlot).

Options:

```text
  VectorPoints    integer n → n × n grid (default 15); Automatic = 15
  VectorScale     Automatic: normalise all arrows to equal display length
                  None:      proportional to magnitude
                  real f:    arrow length = f × grid_spacing
  VectorStyle     style directive(s) applied to all arrows
  ColorFunction   f[vx,vy,speed,x,y] (or fewer args) → color, or "Rainbow"
  ColorFunctionScaling  True (default): normalise speed to [0,1]
  RegionFunction  f[x,y] mask: skip grid points outside the region
  Standard Graphics options pass through (Axes, AspectRatio→1, Frame, …) 
```

## Implementation notes

**Algorithm.** `builtin_vectorplot` is `HoldAll`; `args[0]` must be a 2-element
list `{vx, vy}`. The field is compiled all-or-nothing by `vp_compile` into two
`AutoCompiled` kernels `f(x, y)` (`VpField`), with `vp_eval` falling back to the
interpreter per sample. It samples an **`N x N` grid** (`N = VectorPoints`,
default **15**) in scaled world space, inverting each point to data space
(`scale_invert`) before evaluation and masking it with `RegionFunction`; when a
`ScalingFunctions` transform is active the arrow direction is corrected by a
finite-difference Jacobian. Arrow sizing is **screen-normalised** so arrows stay
legible across mixed-scale axes: directions are divided by the per-axis world
units, and `half_len_screen = 0.5/(N-1)` (half a grid cell) sets the base length.
`VectorScale` selects the length rule — `Automatic` gives a fixed length,
`None` makes length proportional to magnitude (`len ∝ speed/speed_max`), a real
`f` scales relative to the grid spacing. Each arrow is **centred** on its grid
point (`tail = grid - d`, `head = grid + d`) and emitted as
`Arrow[{{tail}, {head}}]` preceded by its colour directive; `vp_color` defaults to
the `default_ramp_rgb` (Viridis) ramp keyed to the normalised speed (a user
function is tried at arity 3 `f[vx, vy, speed]` then arity 1 `f[speed]`). The
result is an inert `Graphics[prims, opts..., PlotRange -> {...}]`.

**Data structures.** `VpField` (the two compiled kernels); a local
`Arrow {wx, wy, vx, vy}` array of `N*N` world-space entries; the option bundle
from `split_vector_options`. The primitive list is a leading `VectorStyle` (or the
default `Thickness[0.005]`) then a colour directive and `Arrow` per grid point.

**Complexity / limits.** `O(N^2)` field evaluations and arrows (default 225). The
grid is always square (no per-axis point count), and only `Arrow` primitives are
emitted. Defaults `Frame -> True`, `Axes -> False`, `AspectRatio -> 1`,
`ColorFunctionScaling -> True`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/vectorplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/vectorplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

### Notes

`VectorPlot` is `HoldAll`. It samples `{vx, vy}` on an `N x N` grid (`VectorPoints`
default 15) and draws one `Arrow` per grid point, **centred** on that point. Arrow
sizing is screen-normalised so arrows stay legible across mixed-scale axes:
`VectorScale -> Automatic` (default) gives a fixed length, `None` makes length
proportional to magnitude, and a real value scales relative to the grid spacing.

The default `ColorFunction` is the Viridis ramp keyed to speed; a custom function
is tried as `f[vx, vy, speed]` then `f[speed]`. When a `ScalingFunctions` transform
is active the arrow direction is corrected by a finite-difference Jacobian. The grid
is always square (there is no per-axis point count).
