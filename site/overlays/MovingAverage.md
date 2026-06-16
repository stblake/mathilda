### Worked examples

```mathematica
In[1]:= MovingAverage[{1, 2, 3, 4, 5}, 2]
Out[1]= {3/2, 5/2, 7/2, 9/2}
```

```mathematica
In[1]:= MovingAverage[{1, 2, 3, 4, 5, 6}, 3]
Out[1]= {2, 3, 4, 5}
```

```mathematica
In[1]:= MovingAverage[Table[k^2, {k, 1, 6}], 3]
Out[1]= {14/3, 29/3, 50/3, 77/3}
```

```mathematica
In[1]:= MovingAverage[{a, b, c, d}, {1, 2, 1}]
Out[1]= {1/4 a + 1/2 b + 1/4 c, 1/4 b + 1/2 c + 1/4 d}
```

### Notes

`MovingAverage[list, r]` slides a window of `r` consecutive elements across the
list, returning their averages; the output has length `Length[list] - r + 1`.
Averages are exact rationals, so smoothing the first six squares with a width-3
window gives `{14/3, 29/3, 50/3, 77/3}` rather than decimals. The list-of-weights
form `MovingAverage[list, {w1, ..., wr}]` performs a weighted moving average with
effective weights `wi / Sum[wj]`; with symbolic data and weights `{1, 2, 1}` it
produces the exact binomial smoothing kernel `a/4 + b/2 + c/4`, the discrete
analogue of a triangular filter. The call is left unevaluated when `r < 1` or
`r > Length[list]`.
