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
Out[1]= {0.935057, 0.634021, 0.644893, 0.115526, 0.856989, 0.425172, 0.561837, 0.992007, 0.473169, 0.900769, 0.505084, 0.340315, 0.493009, 0.824974, 0.868976, 0.595398, 0.775846, 0.0612878, 0.566164, 0.848359, 0.388734, 0.942848, 0.273556, 0.870697, 0.420343, 0.285931, 0.0184189, 0.779514, 0.566537, 0.348367, 0.0254825, 0.756417, 0.866588, 0.663862, 0.569521, 0.225607, 0.75991, 0.148887, 0.298208, 0.518048, 0.635251, 0.956637, 0.499086, 0.459972, 0.0400711, 0.790538, 0.620301, 0.287987, 0.327684, 0.401936, 0.363889, 0.992179, 0.767683, 0.0272633, 0.661659, 0.683265, 0.761935, 0.892214, 0.644562, 0.952911, 0.742732, 0.251056, 0.751228, 0.971595, 0.370634, 0.723145, 0.224784, 0.883396, 0.201982, 0.463945, 0.180656, 0.701665, 0.169295, 0.115859, 0.822238, 0.86132, 0.412563, 0.510327, 0.607433, 0.315704, 0.784243, 0.627649, 0.769239, 0.471218, 0.199051, 0.325602, 0.85973, 0.00745504, 0.102252, 0.150039, 0.852308, 0.982301, 0.011341, 0.314836, 0.613418, 0.205827, 0.769762, 0.996514, 0.591707, 0.190998, 0.32215, 0.15986, 0.627845, 0.301585, 0.85982, 0.558267, 0.599738, 0.355437, 0.065484, 0.993708, 0.347646, 0.163085, 0.744591, 0.87664, 0.249989, 0.511029, 0.513068, 0.916792, 0.345983, 0.213686, 0.809609, 0.297334, 0.833815, 0.987116, 0.686048, 0.719181, 0.555869, 0.714682, 0.200227, 0.391249, 0.384519, 0.62865, 0.171948, 0.56278, 0.875841, 0.923775, 0.584486, 0.847483, 0.755987, 0.468937, 0.412367, 0.384588, 0.557592, 0.648, 0.357371, 0.530778, 0.517535, 0.678552, 0.626559, 0.901912, 0.149731, 0.509581, 0.64509, 0.308181, 0.109504, 0.413558, 0.00773887, 0.564168, 0.9661, 0.692203, 0.815624, 0.757361, 0.237458, 0.226942, 0.0623514, 0.528454, 0.444427, 0.783814, 0.896304, 0.572705, 0.971553, 0.0797881, 0.603702, 0.2863, 0.534893, 0.39327, 0.995832, 0.636654, 0.87645, 0.604095, 0.29418, 0.697998, 0.590837, 0.884401, 0.857162, 0.647563, 0.915632, 0.585124, 0.775605, 0.255002, 0.87915, 0.391699, 0.399368, 0.676775, 0.506986, 0.430726, 0.452697, 0.964788, 0.818754, 0.726982}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
