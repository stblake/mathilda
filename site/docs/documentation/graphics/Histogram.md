# Histogram

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Histogram[data, opts...]`**

Bins the numeric values in data and draws a frequency histogram. Bin count defaults to Sturges' rule: ceil(Log2\[n\]) + 1.

**`Histogram[data, k, opts...]`**

k equal-width bins.

**`Histogram[data, {step}, opts...]`**

Bins of width step.

**`Histogram[data, {min, max, step}, opts...]`**

Explicit range and width. Options: ChartStyle   color/style list cycling through bins BarSpacing   gap fraction of bin width (default 0.2) Standard Graphics options pass through.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Histogram[Table[RandomReal[], {200}]]
Out[1]= -Graphics-

In[2]:= data = Table[RandomReal[], {200}]
Out[2]= {0.381811, 0.565411, 0.517667, 0.958092, 0.341145, 0.575315, 0.854509, 0.172816, 0.277656, 0.291537, 0.585821, 0.350394, 0.808445, 0.408548, 0.562146, 0.0376629, 0.822992, 0.658259, 0.606677, 0.613409, 0.53319, 0.534529, 0.411775, 0.563055, 0.590864, 0.0338609, 0.955115, 0.627521, 0.367303, 0.406122, 0.281652, 0.620763, 0.857902, 0.470922, 0.354334, 0.378084, 0.164892, 0.988439, 0.31415, 0.303961, 0.998969, 0.0621575, 0.265327, 0.642811, 0.527975, 0.343933, 0.717521, 0.174583, 0.150473, 0.153289, 0.230089, 0.36056, 0.10068, 0.934111, 0.465915, 0.338566, 0.968959, 0.9678, 0.0154806, 0.686801, 0.315441, 0.469542, 0.758977, 0.9859, 0.400411, 0.921949, 0.542131, 0.364484, 0.345168, 0.445939, 0.775437, 0.757864, 0.405717, 0.319942, 0.71222, 0.479292, 0.194553, 0.988801, 0.487749, 0.0940605, 0.509245, 0.0981105, 0.171662, 0.665692, 0.794151, 0.859102, 0.902043, 0.854773, 0.771742, 0.123439, 0.32708, 0.730482, 0.430564, 0.818876, 0.292121, 0.693662, 0.649036, 0.295219, 0.767811, 0.543109, 0.0880276, 0.368267, 0.958104, 0.797389, 0.65281, 0.701204, 0.169588, 0.825446, 0.162986, 0.951631, 0.514855, 0.922763, 0.285462, 0.696404, 0.442955, 0.743627, 0.16725, 0.418959, 0.495348, 0.674684, 0.525265, 0.886604, 0.577384, 0.786498, 0.445681, 0.301862, 0.87412, 0.355917, 0.72178, 0.417294, 0.548693, 0.240635, 0.0737463, 0.316564, 0.157612, 0.62727, 0.570123, 0.946624, 0.694695, 0.658579, 0.826288, 0.167716, 0.876017, 0.62743, 0.705158, 0.445195, 0.252827, 0.728106, 0.918171, 0.307353, 0.556157, 0.333418, 0.595037, 0.381433, 0.76755, 0.827423, 0.801922, 0.530996, 0.751559, 0.876634, 0.103755, 0.609448, 0.181726, 0.263179, 0.110045, 0.484638, 0.275366, 0.809255, 0.16446, 0.544327, 0.660907, 0.210212, 0.947442, 0.580215, 0.354538, 0.919552, 0.425326, 0.202524, 0.0783208, 0.3343, 0.785656, 0.257343, 0.641495, 0.0443357, 0.447715, 0.664395, 0.688345, 0.76826, 0.733151, 0.798665, 0.561677, 0.144196, 0.25003, 0.583128, 0.0625727, 0.402078, 0.897964, 0.699562, 0.51601, 0.669173}

In[3]:= Histogram[data, 20]
Out[3]= -Graphics-

In[4]:= Histogram[data, {0.1}]
Out[4]= -Graphics-

In[5]:= Histogram[data, {0, 1, 0.05}]
Out[5]= -Graphics-
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

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [BarSpacing](../../other-advanced/BarSpacing/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
