### Worked examples

```mathematica
In[1]:= MedianDeviation[{1, 2, 3, 4}]  (* the exact median absolute deviation (MAD) of an integer sample *)
Out[1]= 1
```

```mathematica
In[1]:= MedianDeviation[{1., 1., 2., 2., 4., 6., 9.}]  (* a machine-real sample gives a machine-real MAD *)
Out[1]= 1.0
```

```mathematica
In[1]:= MedianDeviation[{{1, 2}, {3, 4}, {5, 9}}]  (* on a matrix the MAD is computed per column *)
Out[1]= {2, 2}
```

### Notes

`MedianDeviation[data] = Median[Abs[data − Median[data]]]`, the median absolute
deviation (MAD) about the median — a robust scale estimator whose breakdown point
is far higher than that of the standard deviation. Composed through the evaluator,
so exact input gives exact output.

As with `MeanDeviation`, non-finite entries (`Infinity`, `ComplexInfinity`,
`Indeterminate`) are rejected rather than left half-evaluated, non-numeric data
triggers `MedianDeviation::rectn`, matrix input is reduced per column, and a
visible `NDArray` is materialised to the exact `List` path.
