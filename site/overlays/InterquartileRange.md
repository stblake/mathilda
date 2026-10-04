### Worked examples

```mathematica
In[1]:= InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]  (* the exact interquartile range of an integer sample *)
Out[1]= 4
```

```mathematica
In[1]:= InterquartileRange[{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.}]  (* a machine-real sample gives a machine-real range *)
Out[1]= 5.0
```

```mathematica
In[1]:= InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]  (* on a matrix the IQR is computed per column *)
Out[1]= {6, 6, 6}
```

### Notes

`InterquartileRange[data]` is `q̂₃ − q̂₁`, the spread between the upper and lower
quartiles, computed by composition over `Quartiles` (and so using the same
`{{1/2, 0}, {0, 1}}` parameterization). It is a robust scale estimator — insensitive
to the data beyond the quartiles — and so a good companion to `MedianDeviation`.

The matrix case is handled before the vector case on purpose: a three-column
matrix's per-column quartiles are a 3-list of triples, which a naive shape test
would confuse with a vector's `{q₁, q₂, q₃}`. The head requires numeric data
(`InterquartileRange::rectn` otherwise) and materialises a visible `NDArray` to the
exact `List` path.
