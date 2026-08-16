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
Out[2]= {0.554875, 0.218365, 0.657306, 0.64324, 0.55722, 0.43154, 0.53388, 0.100799, 0.514862, 0.764987, 0.584183, 0.495272, 0.0237591, 0.0938345, 0.70125, 0.916055, 0.526967, 0.551743, 0.0914545, 0.52662, 0.43794, 0.982179, 0.901229, 0.404792, 0.281995, 0.669899, 0.516813, 0.319384, 0.00856644, 0.300977, 0.786372, 0.109682, 0.333825, 0.083378, 0.481075, 0.417684, 0.18706, 0.0266203, 0.733357, 0.990577, 0.0365065, 0.120842, 0.0588097, 0.18755, 0.165506, 0.00472379, 0.81303, 0.101897, 0.888827, 0.99712, 0.855612, 0.556146, 0.711846, 0.948466, 0.975534, 0.8093, 0.0795143, 0.622707, 0.034899, 0.834846, 0.547253, 0.228818, 0.178001, 0.498561, 0.517269, 0.188877, 0.0743773, 0.802767, 0.340847, 0.75851, 0.635678, 0.131723, 0.566611, 0.560674, 0.500125, 0.106034, 0.0782597, 0.722787, 0.52095, 0.892778, 0.56952, 0.623875, 0.542694, 0.394523, 0.42208, 0.616541, 0.317404, 0.120341, 0.488656, 0.61699, 0.845778, 0.0879759, 0.155211, 0.862262, 0.361101, 0.765399, 0.707553, 0.0806033, 0.102181, 0.0225218, 0.953231, 0.55852, 0.923704, 0.470853, 0.377672, 0.659791, 0.145015, 0.310377, 0.187206, 0.439086, 0.69927, 0.575896, 0.567092, 0.90355, 0.103255, 0.0706889, 0.780991, 0.151375, 0.919043, 0.584185, 0.713403, 0.833646, 0.689108, 0.905419, 0.532909, 0.816955, 0.529775, 0.147382, 0.027402, 0.2614, 0.899222, 0.437986, 0.682946, 0.26333, 0.807281, 0.130715, 0.226847, 0.489876, 0.701486, 0.188382, 0.95335, 0.590264, 0.315105, 0.982294, 0.860138, 0.650955, 0.641946, 0.408754, 0.670437, 0.376975, 0.854414, 0.926719, 0.672693, 0.15365, 0.474625, 0.738794, 0.548092, 0.488122, 0.632123, 0.505887, 0.184101, 0.624077, 0.908285, 0.521583, 0.194807, 0.314453, 0.590455, 0.74103, 0.782324, 0.0341883, 0.800002, 0.686426, 0.609708, 0.540315, 0.571198, 0.361692, 0.531105, 0.499034, 0.477186, 0.797444, 0.443032, 0.446128, 0.259636, 0.202482, 0.458275, 0.858633, 0.428507, 0.818145, 0.750492, 0.123949, 0.292388, 0.0034257, 0.514753, 0.476178, 0.0504431, 0.238885, 0.995454, 0.424672, 0.883347, 0.402372}

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
