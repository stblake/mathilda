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
Out[2]= {0.424917, 0.184099, 0.0464921, 0.324958, 0.310935, 0.481177, 0.926552, 0.566478, 0.467271, 0.735405, 0.22519, 0.956382, 0.980894, 0.182866, 0.553349, 0.27341, 0.080304, 0.333882, 0.748301, 0.845821, 0.0875133, 0.199432, 0.220597, 0.883364, 0.944079, 0.132786, 0.981592, 0.933896, 0.91734, 0.373154, 0.44537, 0.312319, 0.771963, 0.905782, 0.246715, 0.643095, 0.508406, 0.712687, 0.539356, 0.201415, 0.712895, 0.844905, 0.651616, 0.531204, 0.662956, 0.835888, 0.691243, 0.638271, 0.882954, 0.527103, 0.570114, 0.961764, 0.539654, 0.966791, 0.383766, 0.621199, 0.0392646, 0.395594, 0.797025, 0.389332, 0.743853, 0.414341, 0.196954, 0.938111, 0.605745, 0.603165, 0.249204, 0.93443, 0.125843, 0.43182, 0.952917, 0.00794006, 0.577875, 0.420254, 0.745059, 0.482895, 0.315229, 0.849201, 0.219382, 0.104468, 0.622053, 0.330358, 0.242903, 0.326314, 0.338509, 0.844052, 0.121413, 0.671154, 0.0787382, 0.403091, 0.929477, 0.935876, 0.724983, 0.784173, 0.972601, 0.631471, 0.117792, 0.42791, 0.0316047, 0.918928, 0.922752, 0.610133, 0.518691, 0.897481, 0.867212, 0.588986, 0.510219, 0.320482, 0.862074, 0.325684, 0.180259, 0.810277, 0.317588, 0.373461, 0.796287, 0.307068, 0.784014, 0.929455, 0.404778, 0.00302756, 0.311454, 0.684833, 0.238986, 0.470889, 0.424028, 0.785906, 0.403548, 0.711376, 0.289967, 0.0763657, 0.938757, 0.970635, 0.494489, 0.788583, 0.189573, 0.165797, 0.0347499, 0.6163, 0.805765, 0.428065, 0.74841, 0.238969, 0.251565, 0.720248, 0.48168, 0.265513, 0.0798035, 0.839964, 0.914215, 0.36687, 0.608901, 0.103418, 0.309308, 0.47008, 0.649501, 0.428167, 0.186296, 0.856755, 0.596008, 0.785123, 0.698075, 0.0740527, 0.740437, 0.666479, 0.121563, 0.656195, 0.571025, 0.880979, 0.128063, 0.865381, 0.0139424, 0.0599949, 0.568663, 0.0587311, 0.942903, 0.370418, 0.855579, 0.0852149, 0.206909, 0.709839, 0.450389, 0.129512, 0.240619, 0.307779, 0.579937, 0.822378, 0.26431, 0.34313, 0.243926, 0.447062, 0.360335, 0.164046, 0.186445, 0.284328, 0.606725, 0.066421, 0.464345, 0.852019, 0.368492, 0.971551}

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
