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
Out[2]= {0.6388, 0.19753, 0.259622, 0.927308, 0.982672, 0.53752, 0.613184, 0.314762, 0.952705, 0.159625, 0.550893, 0.719175, 0.566538, 0.0558559, 0.970589, 0.196486, 0.190327, 0.506011, 0.00917659, 0.413478, 0.0728979, 0.69732, 0.456341, 0.794324, 0.624654, 0.56591, 0.638684, 0.424455, 0.960227, 0.398391, 0.793824, 0.999916, 0.354205, 0.448606, 0.312365, 0.0737417, 0.299581, 0.0692005, 0.501132, 0.920691, 0.059684, 0.437271, 0.332981, 0.817597, 0.31809, 0.716554, 0.39215, 0.784283, 0.514209, 0.229246, 0.853981, 0.607652, 0.830069, 0.964919, 0.242843, 0.809771, 0.473568, 0.616931, 0.00964103, 0.628178, 0.904607, 0.555132, 0.667616, 0.776532, 0.849907, 0.417223, 0.443795, 0.686178, 0.0873288, 0.354487, 0.194582, 0.122919, 0.348752, 0.627783, 0.723264, 0.892867, 0.663712, 0.666417, 0.302888, 0.814082, 0.134663, 0.303144, 0.356942, 0.9386, 0.679463, 0.817475, 0.0961075, 0.732634, 0.489657, 0.694455, 0.877461, 0.672443, 0.462053, 0.689886, 0.0781771, 0.139573, 0.0172763, 0.300521, 0.158352, 0.0643939, 0.277327, 0.92007, 0.0621353, 0.986212, 0.586585, 0.165351, 0.387966, 0.703436, 0.0628383, 0.185628, 0.872476, 0.678958, 0.466521, 0.0551069, 0.273423, 0.00967196, 0.348305, 0.141463, 0.637832, 0.779823, 0.656071, 0.278478, 0.212413, 0.295241, 0.935449, 0.777237, 0.600295, 0.16894, 0.297807, 0.902459, 0.68654, 0.999017, 0.641403, 0.510027, 0.196824, 0.0636516, 0.0931354, 0.930469, 0.528662, 0.241469, 0.117527, 0.516602, 0.241963, 0.0733072, 0.788241, 0.998809, 0.815579, 0.998205, 0.748668, 0.411133, 0.540389, 0.123337, 0.486672, 0.992346, 0.842496, 0.813884, 0.697966, 0.22221, 0.681753, 0.419355, 0.162449, 0.64675, 0.178569, 0.112894, 0.225653, 0.18889, 0.56818, 0.401918, 0.378958, 0.610042, 0.778432, 0.635925, 0.303168, 0.57785, 0.597052, 0.400409, 0.00502551, 0.13298, 0.291574, 0.304121, 0.969868, 0.757533, 0.457181, 0.321258, 0.960284, 0.990684, 0.759726, 0.630396, 0.479985, 0.589363, 0.0999078, 0.560652, 0.216007, 0.655303, 0.439905, 0.485233, 0.75561, 0.880164, 0.708224, 0.723673}

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
