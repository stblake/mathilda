# BarChart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BarChart[{v1, v2, ..., vn}, opts...]`**

Draws a vertical bar chart: n bars at x = 1..n with heights v1..vn.

**`BarChart[{{v1,...}, {w1,...}, ...}, opts...]`**

Multiple grouped datasets, each in a distinct palette colour. Options: ChartStyle    color/style list cycling through bars (default: palette) ChartLabels   list of x-axis tick labels BarSpacing    gap fraction of bar width (default 0.2) Standard Graphics options (Axes, AspectRatio, Frame, PlotRange, PlotLabel, Background, ImageSize, …) pass through.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BarChart[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[1]= -Graphics-

In[2]:= BarChart[{3, -1, 4, -1, 5}]
Out[2]= -Graphics-
```

### Options (2)

```mathematica
In[3]:= BarChart[{2.5, 4.1, 3.3, 5.7}, ChartStyle -> {Red, Blue, Green, Orange}, ChartLabels -> {"Q1", "Q2", "Q3", "Q4"}]
Out[3]= -Graphics-

In[4]:= BarChart[{{1, 3, 2}, {4, 2, 5}}, BarSpacing -> 0.3]
Out[4]= -Graphics-
```

## Algorithm

barchart.c — BarChart[data, opts...] and Histogram[data, opts...]

BarChart renders a vertical bar chart from explicit heights. Histogram bins numeric data and renders a frequency histogram. Both return Graphics[...] objects auto-displayed by the REPL.

BarChart[{v1,...,vn}, opts...]

```text
  n bars at x = 1..n with heights v1..vn, coloured via ChartStyle or
  the default palette.
```

BarChart[{{v1,...}, {w1,...}, ...}, opts...]

```text
  Multiple grouped datasets, each dataset in a distinct palette colour.
```

Histogram[data, opts...]

```text
Histogram[data, k, opts...]          k equal-width bins
Histogram[data, {step}, opts...]     bins of given width
```

Histogram[data, {min,max,step}, opts...] explicit range + width

Options (both):

```text
  ChartStyle   — color/style list cycling through bars
  ChartLabels  — label list for x-axis ticks
  BarSpacing   — gap fraction of bar width (default 0.2)
  PlotLabel, standard Graphics options pass through 
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [ChartLabels](../../other-advanced/ChartLabels/), [BarSpacing](../../other-advanced/BarSpacing/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
