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
Out[2]= {0.110154, 0.536385, 0.209726, 0.117104, 0.217996, 0.250765, 0.693947, 0.406657, 0.99533, 0.918872, 0.0259855, 0.274381, 0.72438, 0.816681, 0.0189133, 0.516249, 0.691094, 0.10584, 0.245907, 0.546078, 0.818805, 0.718511, 0.525134, 0.493215, 0.560871, 0.296591, 0.0109602, 0.335789, 0.851715, 0.195852, 0.631146, 0.638094, 0.494272, 0.389946, 0.241286, 0.439083, 0.617964, 0.0732778, 0.0254844, 0.478379, 0.212681, 0.124594, 0.937763, 0.433333, 0.835698, 0.971876, 0.33857, 0.994073, 0.421053, 0.943266, 0.848405, 0.258118, 0.421526, 0.399865, 0.530578, 0.860723, 0.437725, 0.364636, 0.240185, 0.30862, 0.113194, 0.611934, 0.0507744, 0.605423, 0.195244, 0.964672, 0.593293, 0.264964, 0.62414, 0.64938, 0.359053, 0.738315, 0.132249, 0.538033, 0.383538, 0.898199, 0.285637, 0.5888, 0.544785, 0.879029, 0.924187, 0.233848, 0.626732, 0.156664, 0.540331, 0.791892, 0.010229, 0.569963, 0.207901, 0.0601243, 0.421945, 0.308857, 0.561248, 0.311915, 0.363324, 0.0728998, 0.263361, 0.584162, 0.0695768, 0.526357, 0.0240933, 0.345001, 0.150012, 0.75595, 0.216702, 0.730842, 0.90744, 0.499862, 0.426239, 0.0534498, 0.894984, 0.716412, 0.843553, 0.0577056, 0.841553, 0.346863, 0.612339, 0.813428, 0.111374, 0.100759, 0.00556426, 0.222562, 0.439042, 0.846834, 0.200445, 0.530478, 0.911489, 0.102426, 0.090909, 0.639977, 0.601994, 0.521627, 0.95291, 0.673776, 0.268722, 0.690841, 0.614664, 0.572302, 0.193269, 0.107276, 0.711278, 0.118106, 0.313389, 0.506721, 0.171901, 0.216412, 0.0463297, 0.941649, 0.999092, 0.0745381, 0.0808771, 0.528814, 0.278071, 0.783544, 0.706967, 0.350348, 0.766091, 0.562589, 0.597498, 0.705887, 0.787569, 0.183986, 0.329786, 0.646151, 0.287432, 0.0273972, 0.618626, 0.00844734, 0.757675, 0.667286, 0.146083, 0.860337, 0.933992, 0.640971, 0.531583, 0.952841, 0.663577, 0.167872, 0.789876, 0.926664, 0.781689, 0.849932, 0.64447, 0.877813, 0.0980876, 0.563882, 0.582613, 0.264647, 0.278893, 0.325666, 0.12457, 0.905553, 0.364737, 0.128283, 0.0840733, 0.683563, 0.045503, 0.701063, 0.154039, 0.481762}

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
