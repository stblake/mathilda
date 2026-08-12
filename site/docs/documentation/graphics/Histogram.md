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
Out[2]= {0.526289, 0.472, 0.544357, 0.658365, 0.290012, 0.725488, 0.224113, 0.151975, 0.703448, 0.599117, 0.773641, 0.297808, 0.929722, 0.807791, 0.729235, 0.732263, 0.03861, 0.100051, 0.31998, 0.732719, 0.724715, 0.0863307, 0.569067, 0.503178, 0.269699, 0.927451, 0.867202, 0.728688, 0.541459, 0.370559, 0.924296, 0.103408, 0.786022, 0.576599, 0.986929, 0.510315, 0.983431, 0.183979, 0.591465, 0.0152483, 0.42877, 0.288045, 0.0665673, 0.897174, 0.533838, 0.15016, 0.918833, 0.341523, 0.987495, 0.0588894, 0.459773, 0.970743, 0.40657, 0.416574, 0.924152, 0.799284, 0.305307, 0.911181, 0.464366, 0.403085, 0.356057, 0.828417, 0.433479, 0.917896, 0.358807, 0.788875, 0.636672, 0.166382, 0.804642, 0.6964, 0.711718, 0.00054279, 0.493511, 0.268272, 0.616367, 0.0671583, 0.296618, 0.741174, 0.726503, 0.216634, 0.880104, 0.397038, 0.799341, 0.833518, 0.188418, 0.00884743, 0.42897, 0.331887, 0.327163, 0.283965, 0.397693, 0.724546, 0.538991, 0.580755, 0.689153, 0.170754, 0.869811, 0.412867, 0.0460911, 0.955811, 0.496067, 0.074444, 0.512732, 0.940383, 0.317539, 0.47965, 0.161954, 0.462769, 0.627353, 0.406953, 0.762491, 0.62538, 0.869074, 0.276709, 0.710354, 0.101994, 0.527207, 0.535576, 0.22591, 0.408916, 0.608574, 0.575922, 0.451533, 0.59066, 0.715457, 0.522181, 0.69728, 0.71477, 0.0914863, 0.219534, 0.331362, 0.836247, 0.0939389, 0.37796, 0.0652581, 0.949439, 0.729128, 0.85615, 0.340822, 0.673495, 0.084362, 0.911173, 0.63715, 0.590734, 0.96925, 0.636834, 0.0143595, 0.194061, 0.166235, 0.781234, 0.582432, 0.643007, 0.309424, 0.0767953, 0.0882414, 0.739588, 0.353917, 0.999266, 0.319815, 0.369322, 0.804583, 0.619722, 0.266344, 0.896887, 0.335134, 0.864189, 0.424188, 0.712183, 0.978518, 0.79096, 0.551305, 0.214105, 0.898633, 0.166936, 0.12087, 0.324863, 0.197909, 0.719056, 0.451345, 0.751701, 0.209231, 0.00616041, 0.336623, 0.758938, 0.744696, 0.428423, 0.826336, 0.413575, 0.419952, 0.765312, 0.69726, 0.818987, 0.259742, 0.249779, 0.0138919, 0.318312, 0.469998, 0.9185, 0.781461, 0.723068}

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
