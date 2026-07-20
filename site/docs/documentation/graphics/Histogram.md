# Histogram

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Histogram[data, opts...]
    Bins the numeric values in data and draws a frequency histogram.
    Bin count defaults to Sturges' rule: ceil(Log2[n]) + 1.
Histogram[data, k, opts...]
    k equal-width bins.
Histogram[data, {step}, opts...]
    Bins of width step.
Histogram[data, {min, max, step}, opts...]
    Explicit range and width.

    Options:
      ChartStyle   color/style list cycling through bins
      BarSpacing   gap fraction of bin width (default 0.2)
      Standard Graphics options pass through.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= data = Table[RandomReal[], {200}]
Out[1]= {0.221688, 0.870888, 0.12876, 0.569992, 0.221897, 0.256753, 0.620016, 0.932122, 0.919757, 0.882481, 0.308026, 0.198228, 0.602547, 0.97653, 0.322893, 0.341688, 0.561818, 0.768395, 0.957481, 0.840476, 0.838587, 0.867618, 0.512536, 0.742414, 0.558439, 0.450289, 0.543026, 0.497358, 0.648099, 0.822109, 0.975153, 0.51641, 0.367635, 0.513985, 0.362038, 0.868808, 0.600727, 0.628272, 0.443499, 0.0892004, 0.331814, 0.542099, 0.331909, 0.641859, 0.609855, 0.391069, 0.0921424, 0.528631, 0.00217859, 0.022109, 0.741695, 0.621449, 0.180581, 0.227325, 0.965069, 0.174037, 0.108688, 0.967555, 0.0571389, 0.743906, 0.949172, 0.662362, 0.952506, 0.535005, 0.853175, 0.542508, 0.636327, 0.858296, 0.678309, 0.903672, 0.425273, 0.51838, 0.57163, 0.198997, 0.928491, 0.398817, 0.0153792, 0.705652, 0.166609, 0.542376, 0.4092, 0.0727793, 0.160092, 0.830068, 0.553478, 0.233212, 0.90503, 0.824535, 0.650188, 0.441774, 0.244004, 0.128008, 0.191017, 0.618797, 0.420244, 0.780165, 0.721903, 0.0797589, 0.378712, 0.957557, 0.521102, 0.979407, 0.464369, 0.163306, 0.63417, 0.287051, 0.061491, 0.883572, 0.776861, 0.530102, 0.183941, 0.461169, 0.0987336, 0.0907663, 0.177086, 0.802317, 0.0352784, 0.437978, 0.87334, 0.478206, 0.0474372, 0.639135, 0.206772, 0.969936, 0.91784, 0.152896, 0.556423, 0.672064, 0.992947, 0.356429, 0.357605, 0.123614, 0.380665, 0.0896029, 0.52867, 0.82454, 0.912847, 0.0379769, 0.202658, 0.424609, 0.381477, 0.486232, 0.819332, 0.104765, 0.838176, 0.0752534, 0.637977, 0.931431, 0.409902, 0.467848, 0.594259, 0.257392, 0.772658, 0.0768824, 0.285639, 0.707626, 0.742571, 0.245426, 0.00998824, 0.0654921, 0.0956331, 0.193023, 0.280781, 0.421955, 0.440124, 0.820011, 0.750137, 0.506523, 0.411173, 0.396663, 0.922134, 0.325052, 0.434298, 0.485811, 0.28888, 0.0266817, 0.848664, 0.514804, 0.167223, 0.0714496, 0.0448969, 0.124686, 0.959857, 0.400977, 0.0826051, 0.0666523, 0.0912369, 0.717309, 0.688934, 0.685154, 0.759059, 0.784748, 0.861009, 0.18395, 0.693747, 0.563032, 0.784622, 0.935085, 0.735978, 0.30624}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
