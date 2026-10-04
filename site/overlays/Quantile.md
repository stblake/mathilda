### Worked examples

```mathematica
In[1]:= Quantile[{1, 2, 3, 4}, 1/2]  (* the default Type-1 quantile selects an element, not an average *)
Out[1]= 2
```

```mathematica
In[1]:= Median[{1, 2, 3, 4}]  (* while Median averages the two central values, so the two deliberately differ *)
Out[1]= 5/2
```

```mathematica
In[1]:= Quantile[{1, 2, 3, 4}, {1/4, 1/2, 3/4}]  (* a list of probabilities gives a list of quantiles *)
Out[1]= {1, 2, 3}
```

```mathematica
In[1]:= Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]  (* the Quartiles parameterization interpolates — here it is the median *)
Out[1]= 5/2
```

```mathematica
In[1]:= Quantile[{5., 1., 4., 2., 3.}, 0.9]  (* data need not be pre-sorted; a machine-real sample gives a machine-real quantile *)
Out[1]= 5.0
```

```mathematica
In[1]:= Quantile[{1., 2., 3., 4., 5.}, {0.25, 0.5, 0.75}]  (* the list form on real data *)
Out[1]= {2.0, 3.0, 4.0}
```

### Notes

The default parameters are `{{0, 0}, {1, 0}}` — Hyndman–Fan **Type 1**, the
left-continuous inverse CDF — so `Quantile[{1, 2, 3, 4}, 1/2]` selects the order
statistic `x₍⌈nq⌉₎ = 2`, whereas `Median` averages the two central values to get
`5/2`. The two are meant to differ on even-length data. `Quartiles[data]` is
`Quantile[data, {1/4, 1/2, 3/4}, {{1/2, 0}, {0, 1}}]`, which is why the fourth
example reproduces the median.

With `h = a + (n + b)q` and weight `w = c + d (h − ⌊h⌋)`, the result is the
edge-clamped interpolation between `x₍⌊h⌋₎` and `x₍⌈h⌉₎`; for `w ∈ [0, 1]` it is
evaluated as the convex combination `(1 − w) x₍⌊h⌋₎ + w x₍⌈h⌉₎`, which (unlike the
algebraically equal difference form) does not overflow on `Real` neighbours that
straddle zero near the machine range. `Quantile` requires real numeric data and
`q ∈ [0, 1]` (`Quantile::q100` otherwise), computes per column for a matrix, and
materialises a visible `NDArray` to the exact `List` path.
