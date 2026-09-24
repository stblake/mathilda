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
Out[2]= {0.69202, 0.690881, 0.709313, 0.150118, 0.187666, 0.945499, 0.276089, 0.567764, 0.680711, 0.797837, 0.779808, 0.28294, 0.610547, 0.331104, 0.453794, 0.651758, 0.930299, 0.807628, 0.016772, 0.0200348, 0.115639, 0.837705, 0.824103, 0.204303, 0.236038, 0.76606, 0.664714, 0.717929, 0.764734, 0.827211, 0.241333, 0.367392, 0.423792, 0.0384548, 0.781493, 0.385961, 0.0805916, 0.211168, 0.669152, 0.860032, 0.176302, 0.159464, 0.784457, 0.939632, 0.781217, 0.947177, 0.831651, 0.450373, 0.0067173, 0.769479, 0.383734, 0.503704, 0.867683, 0.898742, 0.370826, 0.727606, 0.530099, 0.834414, 0.21314, 0.076932, 0.840466, 0.723911, 0.568582, 0.739547, 0.483485, 0.971621, 0.894402, 0.00227273, 0.564126, 0.291537, 0.497164, 0.416576, 0.609903, 0.920561, 0.290441, 0.131721, 0.184904, 0.273455, 0.510234, 0.430705, 0.631596, 0.787577, 0.54864, 0.692714, 0.893727, 0.385517, 0.690127, 0.995913, 0.105539, 0.241375, 0.411758, 0.508411, 0.326449, 0.436279, 0.570323, 0.997599, 0.434812, 0.838711, 0.632932, 0.239633, 0.314202, 0.658179, 0.554992, 0.0675518, 0.0426088, 0.673519, 0.576863, 0.939458, 0.910504, 0.131689, 0.523862, 0.383159, 0.383685, 0.461691, 0.950721, 0.99322, 0.414799, 0.197595, 0.83013, 0.272163, 0.837147, 0.08432, 0.392807, 0.195173, 0.263356, 0.790779, 0.0280253, 0.718866, 0.00586651, 0.773467, 0.201187, 0.669678, 0.18843, 0.429096, 0.18929, 0.0200408, 0.202214, 0.381972, 0.12761, 0.755167, 0.685388, 0.595682, 0.316822, 0.834209, 0.879947, 0.12075, 0.113458, 0.904761, 0.878424, 0.662946, 0.502443, 0.462434, 0.293352, 0.617551, 0.240673, 0.29864, 0.478502, 0.305123, 0.177637, 0.319511, 0.241068, 0.745334, 0.153086, 0.433783, 0.603551, 0.480742, 0.754855, 0.828404, 0.734792, 0.358508, 0.0659264, 0.709785, 0.151552, 0.103116, 0.575402, 0.476557, 0.0299329, 0.0384641, 0.688377, 0.713942, 0.930018, 0.824514, 0.558687, 0.421234, 0.853683, 0.51921, 0.895426, 0.433297, 0.838034, 0.312436, 0.559592, 0.432983, 0.632055, 0.856837, 0.602149, 0.645816, 0.411279, 0.364423, 0.450502, 0.935041}

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
