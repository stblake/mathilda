### Worked examples

```mathematica
In[1]:= MeanDeviation[{1, 2, 3, 4}]  (* the exact mean absolute deviation of an integer sample *)
Out[1]= 1
```

```mathematica
In[1]:= MeanDeviation[{1/2, 3/2}]  (* exact rationals stay exact *)
Out[1]= 1/2
```

```mathematica
In[1]:= MeanDeviation[{2., 4., 4., 4., 5., 5., 7., 9.}]  (* a machine-real sample gives a machine-real deviation *)
Out[1]= 1.5
```

```mathematica
In[1]:= MeanDeviation[{{1, 2}, {3, 4}, {5, 9}}]  (* on a matrix the deviation is computed per column *)
Out[1]= {4/3, 8/3}
```

### Notes

`MeanDeviation[data] = Mean[Abs[data − Mean[data]]]`, the mean absolute deviation
about the mean. It is composed through the evaluator, so exact input gives exact
output. Non-finite entries (`Infinity`, `ComplexInfinity`, `Indeterminate`) are
rejected even though they are `NumericQ` and free of the imaginary unit — they
would otherwise leave the composed `Mean` half-evaluated — and non-numeric data
triggers `MeanDeviation::rectn`.

Matrix input is reduced per column; a visible `NDArray` is materialised to the
exact `List` path.
