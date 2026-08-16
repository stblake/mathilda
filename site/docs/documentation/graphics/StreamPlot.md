# StreamPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Traces streamlines of the 2-D vector field {vx, vy} by RK4 integration from a grid of seed points, and returns a Graphics\[{Arrow\[...\], ...}, opts\] object (auto-displayed). StreamPlot is HoldAll: vx, vy, and the iterator specs are held unevaluated until x and y are given numeric values. Options: StreamPoints  - Integer n (n x n seed grid) or Automatic (default 15 x 15). StreamScale   – Automatic (8%% of domain diagonal, default), None (full run), or a real fraction of the domain diagonal. StreamStyle   – Style directive(s) applied to all streams. StreamColorFunction / ColorFunction – f\[x,y,vx,vy,speed\] (or fewer args) returning a color, or a named ramp: "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" (all keyed to scaled speed). RegionFunction – f\[x,y\] mask; seeds outside the region are skipped. PlotLegends   – Automatic / "Expressions" / explicit label list. StreamAnimate  – True | False (default). When True, each streamline is emitted as AnimatedStreamline\[...\] instead of Line\[...\]: the shape is drawn identically, but an interactive window (Show, or embedded in Animate/Manipulate) also draws particles flowing along it in real time. Standard Graphics options (PlotRange, Axes, AspectRatio, Frame, …) pass through to the Graphics\[...\] result.

## Examples (10)

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

## Algorithm

streamplot.c — StreamPlot[{vx,vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]

Traces streamlines of a 2-D vector field by RK4 integration from a grid

```text
of seed points, emitting one Arrow[...] primitive per stream.  The result
```

is returned as a Graphics[...] object (auto-displayed by the REPL).

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Show](../../graphics/Show/), [Animate](../../graphics/Animate/), [Manipulate](../../graphics/Manipulate/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_streamplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streamplot.c)
