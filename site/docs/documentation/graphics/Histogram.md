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
Out[2]= {0.237231, 0.374445, 0.565766, 0.409079, 0.548907, 0.21016, 0.652042, 0.621096, 0.903705, 0.244033, 0.998326, 0.663092, 0.487048, 0.988259, 0.285837, 0.115686, 0.201192, 0.695576, 0.257594, 0.469508, 0.612651, 0.824651, 0.804327, 0.444043, 0.0334161, 0.955305, 0.258935, 0.20425, 0.134554, 0.699016, 0.468779, 0.391009, 0.961848, 0.519307, 0.254346, 0.137651, 0.801584, 0.660504, 0.864948, 0.586927, 0.917742, 0.936614, 0.544408, 0.527664, 0.29718, 0.802154, 0.754738, 0.059492, 0.150394, 0.628085, 0.583602, 0.700326, 0.284303, 0.337174, 0.63721, 0.832098, 0.411809, 0.304355, 0.843128, 0.92128, 0.817416, 0.937158, 0.253789, 0.618188, 0.79216, 0.530865, 0.365384, 0.851693, 0.530585, 0.615109, 0.423994, 0.933485, 0.445516, 0.974982, 0.32607, 0.700611, 0.756496, 0.499694, 0.750814, 0.458184, 0.573259, 0.0108581, 0.330942, 0.647395, 0.511459, 0.16486, 0.679805, 0.321259, 0.977443, 0.0700589, 0.651091, 0.515003, 0.370104, 0.581186, 0.894053, 0.48973, 0.174798, 0.815803, 0.280629, 0.131743, 0.948266, 0.759397, 0.342507, 0.576124, 0.328183, 0.956571, 0.56431, 0.606028, 0.937722, 0.751981, 0.990932, 0.993808, 0.443058, 0.776485, 0.780467, 0.551943, 0.198564, 0.451093, 0.660177, 0.630649, 0.80991, 0.146543, 0.577419, 0.220213, 0.158129, 0.334406, 0.364532, 0.539079, 0.970147, 0.800423, 0.2774, 0.427133, 0.164928, 0.572882, 0.718588, 0.858568, 0.192241, 0.876377, 0.183892, 0.0452708, 0.479306, 0.929045, 0.0909206, 0.24108, 0.391968, 0.133635, 0.786455, 0.45484, 0.0575281, 0.811539, 0.651936, 0.663715, 0.134136, 0.816239, 0.502966, 0.453494, 0.503355, 0.555362, 0.997315, 0.704096, 0.305336, 0.671408, 0.0413226, 0.196484, 0.906566, 0.086482, 0.689731, 0.911079, 0.711688, 0.623403, 0.110401, 0.529685, 0.0219178, 0.106258, 0.0565692, 0.138183, 0.927996, 0.84653, 0.498723, 0.393964, 0.494825, 0.64314, 0.0553976, 0.862018, 0.605671, 0.372756, 0.425718, 0.979043, 0.0176325, 0.971866, 0.705885, 0.504581, 0.333668, 0.959315, 0.483912, 0.273984, 0.206439, 0.313656, 0.0576264, 0.52576}

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
