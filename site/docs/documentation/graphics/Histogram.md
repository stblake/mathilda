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
Out[2]= {0.163459, 0.399881, 0.180261, 0.488293, 0.364611, 0.546171, 0.799625, 0.567222, 0.275357, 0.518935, 0.854059, 0.891557, 0.609034, 0.893852, 0.945741, 0.203367, 0.0153462, 0.552036, 0.82574, 0.120352, 0.543773, 0.718595, 0.818073, 0.647427, 0.157353, 0.336319, 0.72078, 0.0387531, 0.871289, 0.949323, 0.813373, 0.963234, 0.918982, 0.73235, 0.0332513, 0.363044, 0.0295465, 0.449457, 0.908902, 0.660356, 0.848165, 0.119556, 0.174597, 0.93679, 0.681477, 0.92874, 0.578528, 0.332455, 0.0199278, 0.879728, 0.404211, 0.657884, 0.118411, 0.232114, 0.144167, 0.789001, 0.459248, 0.640935, 0.772508, 0.586038, 0.969974, 0.191942, 0.458281, 0.75289, 0.52573, 0.868656, 0.00643903, 0.825028, 0.374633, 0.538904, 0.516429, 0.636781, 0.143862, 0.252373, 0.471298, 0.514287, 0.32073, 0.46356, 0.661759, 0.176205, 0.908985, 0.390064, 0.831806, 0.141807, 0.922481, 0.104998, 0.551512, 0.42463, 0.25566, 0.490701, 0.867931, 0.991217, 0.164728, 0.365557, 0.0904831, 0.501906, 0.0848454, 0.702384, 0.0306247, 0.31768, 0.434679, 0.235537, 0.752202, 0.0640635, 0.826504, 0.326468, 0.71664, 0.244142, 0.689698, 0.54554, 0.405133, 0.385448, 0.99899, 0.836612, 0.246385, 0.851748, 0.840174, 0.817124, 0.0592616, 0.234919, 0.198897, 0.479217, 0.483222, 0.509122, 0.995372, 0.114941, 0.409288, 0.558008, 0.828863, 0.979411, 0.537475, 0.0478138, 0.858764, 0.550252, 0.594273, 0.954652, 0.20372, 0.850857, 0.789125, 0.169924, 0.991752, 0.962901, 0.373489, 0.957726, 0.383316, 0.896055, 0.0842293, 0.271931, 0.295567, 0.935604, 0.454057, 0.674897, 0.322378, 0.843365, 0.163633, 0.635681, 0.785748, 0.366094, 0.131396, 0.69922, 0.369056, 0.685819, 0.25458, 0.823105, 0.0108974, 0.582242, 0.607249, 0.465932, 0.798227, 0.0578118, 0.36602, 0.600256, 0.742257, 0.571475, 0.589142, 0.898007, 0.736182, 0.928651, 0.679331, 0.759978, 0.618506, 0.995947, 0.296229, 0.338429, 0.926193, 0.0183857, 0.0732505, 0.997971, 0.0812398, 0.100212, 0.26097, 0.0796988, 0.936997, 0.448482, 0.30116, 0.967142, 0.492933, 0.222438, 0.891999, 0.727065}

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
