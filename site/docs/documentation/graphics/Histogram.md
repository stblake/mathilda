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
Out[2]= {0.629586, 0.942751, 0.605166, 0.600698, 0.984104, 0.736091, 0.362343, 0.0536784, 0.659348, 0.343264, 0.990965, 0.259023, 0.658899, 0.598801, 0.376843, 0.159839, 0.824332, 0.275455, 0.819215, 0.644934, 0.327485, 0.544069, 0.728272, 0.254794, 0.809249, 0.913643, 0.274282, 0.233561, 0.88067, 0.537212, 0.000305119, 0.374491, 0.108083, 0.879927, 0.928964, 0.347022, 0.661534, 0.664651, 0.118519, 0.305813, 0.065999, 0.677769, 0.756209, 0.688014, 0.634646, 0.230511, 0.0475002, 0.77173, 0.51776, 0.831674, 0.970246, 0.64772, 0.251662, 0.879048, 0.066056, 0.0177714, 0.984993, 0.105977, 0.521904, 0.158469, 0.30772, 0.24338, 0.282198, 0.061153, 0.611993, 0.661499, 0.252708, 0.419708, 0.181371, 0.532849, 0.0350252, 0.950103, 0.127665, 0.813139, 0.601906, 0.50915, 0.51508, 0.992201, 0.391273, 0.782088, 0.314878, 0.772962, 0.810324, 0.227198, 0.494865, 0.134716, 0.867561, 0.783321, 0.293802, 0.0529513, 0.485413, 0.737276, 0.478325, 0.0369635, 0.258522, 0.291693, 0.488339, 0.0857743, 0.979248, 0.171417, 0.485501, 0.157999, 0.737085, 0.414572, 0.885217, 0.795795, 0.00136579, 0.770361, 0.32317, 0.0637698, 0.427699, 0.242044, 0.391411, 0.751366, 0.356184, 0.811324, 0.575526, 0.724431, 0.308428, 0.551691, 0.0245983, 0.370749, 0.130533, 0.578542, 0.153674, 0.989946, 0.45492, 0.893995, 0.295944, 0.211936, 0.837487, 0.227676, 0.685531, 0.351306, 0.0630368, 0.740355, 0.823377, 0.689718, 0.546103, 0.722287, 0.47965, 0.181514, 0.156719, 0.762385, 0.276024, 0.0961886, 0.878752, 0.734601, 0.137373, 0.421853, 0.458729, 0.720099, 0.10634, 0.232878, 0.682204, 0.684563, 0.792702, 0.0513911, 0.814226, 0.692252, 0.609212, 0.801752, 0.469385, 0.564076, 0.928077, 0.325362, 0.0441925, 0.772785, 0.282269, 0.00356752, 0.814443, 0.430127, 0.0194799, 0.315888, 0.483764, 0.228432, 0.255796, 0.0231939, 0.587596, 0.429218, 0.297711, 0.329079, 0.754798, 0.251364, 0.324147, 0.978771, 0.873615, 0.198789, 0.464512, 0.957097, 0.21785, 0.492796, 0.152478, 0.591044, 0.429971, 0.927662, 0.753379, 0.675894, 0.603158, 0.545611}

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
