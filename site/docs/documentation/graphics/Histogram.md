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
Out[1]= {0.899575, 0.948054, 0.628896, 0.589763, 0.631753, 0.310847, 0.973546, 0.307908, 0.0260113, 0.970786, 0.251562, 0.587218, 0.859885, 0.0528063, 0.640595, 0.133445, 0.93762, 0.274082, 0.942336, 0.436332, 0.636595, 0.67246, 0.0770144, 0.825534, 0.139753, 0.735146, 0.260116, 0.304597, 0.898637, 0.391556, 0.590283, 0.934026, 0.741103, 0.345415, 0.0526863, 0.968573, 0.399691, 0.817443, 0.309447, 0.807634, 0.38809, 0.355142, 0.636883, 0.0191117, 0.361295, 0.992368, 0.137802, 0.0445462, 0.944068, 0.387398, 0.0597762, 0.367451, 0.575104, 0.492564, 0.836797, 0.624523, 0.208622, 0.529702, 0.983735, 0.658633, 0.272116, 0.120605, 0.803297, 0.550474, 0.845999, 0.695447, 0.538163, 0.188346, 0.625114, 0.939472, 0.9025, 0.659236, 0.155097, 0.044435, 0.0762484, 0.126599, 0.28902, 0.564247, 0.512105, 0.467639, 0.334157, 0.853519, 0.27176, 0.562462, 0.0281943, 0.844045, 0.952635, 0.96774, 0.818638, 0.475356, 0.0513118, 0.295332, 0.374502, 0.877734, 0.329813, 0.412349, 0.450765, 0.0841497, 0.0856299, 0.13325, 0.749835, 0.688318, 0.472783, 0.652101, 0.0279092, 0.932337, 0.765428, 0.0308056, 0.78999, 0.291938, 0.79753, 0.291866, 0.00729574, 0.821765, 0.184492, 0.938578, 0.721652, 0.0841076, 0.217988, 0.363826, 0.929263, 0.337202, 0.953892, 0.97948, 0.200339, 0.29495, 0.964053, 0.0493211, 0.822055, 0.209801, 0.575922, 0.812189, 0.312831, 0.821608, 0.549197, 0.687417, 0.35369, 0.960378, 0.292244, 0.765492, 0.0809799, 0.0575503, 0.317469, 0.353929, 0.448358, 0.337742, 0.347479, 0.917663, 0.540326, 0.608017, 0.462837, 0.528178, 0.967624, 0.257446, 0.539685, 0.0948376, 0.619453, 0.971002, 0.405475, 0.71745, 0.735122, 0.163672, 0.912248, 0.412413, 0.0804191, 0.497789, 0.669806, 0.186231, 0.550706, 0.870948, 0.468727, 0.0712072, 0.796365, 0.00600198, 0.779442, 0.202698, 0.405602, 0.141533, 0.986462, 0.232624, 0.677284, 0.220347, 0.69941, 0.030514, 0.973683, 0.475246, 0.0994546, 0.0309971, 0.804714, 0.336263, 0.392022, 0.706637, 0.147014, 0.585064, 0.42953, 0.243266, 0.856627, 0.142376, 0.132674, 0.404484}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
