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
Out[2]= {0.196111, 0.245354, 0.565619, 0.594431, 0.179752, 0.366213, 0.941119, 0.289687, 0.610268, 0.351385, 0.0360112, 0.808765, 0.0219158, 0.726311, 0.821885, 0.505053, 0.393601, 0.802003, 0.486629, 0.0745378, 0.293095, 0.721279, 0.0959323, 0.189982, 0.458306, 0.197874, 0.682533, 0.865253, 0.48575, 0.952798, 0.656869, 0.183761, 0.160079, 0.563724, 0.218579, 0.477994, 0.711483, 0.98865, 0.772193, 0.098602, 0.804028, 0.0478445, 0.0140912, 0.358531, 0.12685, 0.386476, 0.415033, 0.60857, 0.212079, 0.584264, 0.31878, 0.692012, 0.55804, 0.758266, 0.568496, 0.298379, 0.914161, 0.631827, 0.703277, 0.412725, 0.780629, 0.831801, 0.531885, 0.770318, 0.744818, 0.870042, 0.966012, 0.872791, 0.0439762, 0.635564, 0.83381, 0.563948, 0.450552, 0.471029, 0.480525, 0.167484, 0.838567, 0.464977, 0.550288, 0.64475, 0.587563, 0.661885, 0.743214, 0.736044, 0.990045, 0.640839, 0.696208, 0.448824, 0.633417, 0.875786, 0.265168, 0.97023, 0.96676, 0.885306, 0.0864834, 0.575181, 0.175961, 0.085534, 0.0852377, 0.190319, 0.915382, 0.352276, 0.102539, 0.356451, 0.158537, 0.838439, 0.610443, 0.425603, 0.935899, 0.986523, 0.246365, 0.462532, 0.0853146, 0.692733, 0.372798, 0.34892, 0.565835, 0.327587, 0.577304, 0.326466, 0.20986, 0.77368, 0.17855, 0.577021, 0.891355, 0.0532481, 0.154128, 0.283284, 0.205098, 0.945791, 0.663245, 0.506918, 0.330698, 0.328996, 0.30982, 0.106281, 0.809025, 0.497852, 0.548712, 0.548405, 0.64941, 0.814545, 0.448213, 0.699037, 0.416372, 0.193473, 0.983609, 0.034367, 0.850911, 0.418659, 0.0743627, 0.318902, 0.715434, 0.251744, 0.0577854, 0.197381, 0.727104, 0.136417, 0.211629, 0.779168, 0.954151, 0.497768, 0.886959, 0.955738, 0.500886, 0.958923, 0.822482, 0.44359, 0.69157, 0.12148, 0.786116, 0.781982, 0.0694519, 0.264029, 0.907978, 0.7453, 0.0927259, 0.381476, 0.867787, 0.230869, 0.294446, 0.737185, 0.502735, 0.806292, 0.270766, 0.437318, 0.798591, 0.262998, 0.163999, 0.43825, 0.511611, 0.0485982, 0.567658, 0.479723, 0.767343, 0.938719, 0.382932, 0.916514, 0.268318, 0.179268}

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

## See also

[ChartStyle](../../other-advanced/ChartStyle/), [BarSpacing](../../other-advanced/BarSpacing/)

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
