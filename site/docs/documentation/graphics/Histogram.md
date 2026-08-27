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
Out[2]= {0.244163, 0.713314, 0.224976, 0.0639541, 0.909797, 0.0701394, 0.127827, 0.805999, 0.124084, 0.408183, 0.163501, 0.781861, 0.847288, 0.848115, 0.101316, 0.5252, 0.916502, 0.0678374, 0.571809, 0.995846, 0.83978, 0.935903, 0.771458, 0.973352, 0.464493, 0.781702, 0.579408, 0.00276135, 0.983198, 0.369131, 0.491101, 0.396882, 0.727233, 0.842599, 0.0898623, 0.247393, 0.482593, 0.904608, 0.450595, 0.0858642, 0.795125, 0.313891, 0.0706362, 0.729018, 0.864282, 0.195639, 0.626901, 0.0612359, 0.117154, 0.257349, 0.261416, 0.556631, 0.880692, 0.734868, 0.692397, 0.0940734, 0.924581, 0.846958, 0.905834, 0.420419, 0.304604, 0.0797307, 0.229237, 0.321182, 0.759437, 0.338299, 0.130389, 0.283755, 0.594734, 0.213802, 0.615958, 0.835749, 0.186752, 0.691226, 0.874303, 0.950438, 0.489498, 0.0296101, 0.258212, 0.977498, 0.697707, 0.520392, 0.610039, 0.133308, 0.520513, 0.933402, 0.195766, 0.994272, 0.109608, 0.377394, 0.825253, 0.724692, 0.727846, 0.420205, 0.142554, 0.900651, 0.797009, 0.00253511, 0.919124, 0.380204, 0.464259, 0.104319, 0.183781, 0.0817113, 0.584081, 0.865813, 0.0683096, 0.330856, 0.329937, 0.975004, 0.972962, 0.0985559, 0.566805, 0.921521, 0.25581, 0.19233, 0.178051, 0.0479284, 0.847013, 0.646028, 0.302099, 0.199475, 0.517436, 0.964877, 0.556736, 0.723791, 0.861801, 0.230068, 0.494914, 0.596513, 0.649441, 0.557038, 0.578154, 0.26336, 0.628363, 0.025725, 0.608854, 0.895514, 0.495079, 0.745635, 0.528465, 0.270966, 0.146849, 0.298602, 0.667878, 0.111662, 0.368146, 0.39874, 0.908115, 0.504391, 0.89932, 0.833793, 0.365466, 0.29501, 0.725555, 0.216676, 0.7082, 0.694024, 0.775463, 0.292231, 0.699224, 0.0969682, 0.259102, 0.150256, 0.949835, 0.862009, 0.0956932, 0.352045, 0.737099, 0.329791, 0.817544, 0.61213, 0.538152, 0.885492, 0.456391, 0.766894, 0.298947, 0.848666, 0.604077, 0.392653, 0.76608, 0.197935, 0.636469, 0.179744, 0.655004, 0.764894, 0.404163, 0.656621, 0.77183, 0.789267, 0.349202, 0.653627, 0.493304, 0.285075, 0.852015, 0.108628, 0.272836, 0.57385, 0.248562, 0.0718005}

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
