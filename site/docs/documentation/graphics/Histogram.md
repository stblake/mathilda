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
Out[2]= {0.312279, 0.0644816, 0.278431, 0.0973643, 0.578796, 0.802526, 0.596102, 0.969221, 0.458642, 0.388538, 0.0485295, 0.402131, 0.294593, 0.465992, 0.255709, 0.19374, 0.493489, 0.95429, 0.00454312, 0.281325, 0.790953, 0.223021, 0.961685, 0.0708533, 0.273226, 0.423613, 0.483569, 0.0763389, 0.916729, 0.568841, 0.619876, 0.441259, 0.417793, 0.955896, 0.744325, 0.915519, 0.833213, 0.39709, 0.071091, 0.331111, 0.995099, 0.4961, 0.270289, 0.902504, 0.927902, 0.666405, 0.352448, 0.226185, 0.855236, 0.204564, 0.349799, 0.237089, 0.23802, 0.741803, 0.561523, 0.907668, 0.21117, 0.581708, 0.741152, 0.235531, 0.963891, 0.326093, 0.0565399, 0.812333, 0.561723, 0.222077, 0.702557, 0.0263261, 0.279356, 0.584519, 0.869087, 0.120552, 0.664597, 0.752338, 0.585846, 0.430432, 0.0234868, 0.231097, 0.199855, 0.757547, 0.3361, 0.484743, 0.558989, 0.802054, 0.267541, 0.930581, 0.114471, 0.320718, 0.869414, 0.97631, 0.399969, 0.232712, 0.745537, 0.67135, 0.472496, 0.884052, 0.290778, 0.215283, 0.673366, 0.99302, 0.817226, 0.364785, 0.73505, 0.412998, 0.483558, 0.520168, 0.483019, 0.354052, 0.0285501, 0.495102, 0.423781, 0.867705, 0.257575, 0.999019, 0.362843, 0.390625, 0.703903, 0.194072, 0.367243, 0.354545, 0.450235, 0.501545, 0.689988, 0.0411663, 0.545247, 0.90909, 0.460224, 0.503369, 0.199828, 0.64877, 0.328311, 0.707156, 0.352909, 0.635605, 0.027416, 0.820036, 0.635003, 0.714314, 0.271235, 0.378361, 0.412376, 0.26072, 0.330885, 0.502162, 0.278839, 0.593518, 0.0885562, 0.610705, 0.266956, 0.774154, 0.093885, 0.654604, 0.812832, 0.824252, 0.385939, 0.185986, 0.735021, 0.54046, 0.577888, 0.0558794, 0.645577, 0.690518, 0.44838, 0.590332, 0.0500516, 0.0649599, 0.220754, 0.384908, 0.642317, 0.508318, 0.570748, 0.372219, 0.937311, 0.554069, 0.299197, 0.0478162, 0.448979, 0.815344, 0.0240904, 0.272927, 0.525415, 0.514603, 0.894152, 0.929939, 0.748466, 0.738843, 0.232639, 0.78562, 0.0920279, 0.134826, 0.380744, 0.805587, 0.929174, 0.375384, 0.441936, 0.290989, 0.0379912, 0.393053, 0.111966, 0.479356}

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
