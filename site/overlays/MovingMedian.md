### Worked examples

```mathematica
In[1]:= MovingMedian[{1, 3, 2, 8, 5, 4, 9}, 3]
Out[1]= {2, 3, 5, 5, 5}
```

```mathematica
In[1]:= MovingMedian[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 5]
Out[1]= {3, 4, 5, 6, 7, 8}
```

```mathematica
In[1]:= MovingMedian[Table[Mod[n^2, 11], {n, 1, 12}], 4]
Out[1]= {9/2, 9/2, 4, 4, 4, 9/2, 9/2, 5/2, 1}
```

### Notes

`MovingMedian[list, r]` returns a list of length `Length[list] - r + 1`: the
median of each window of `r` consecutive elements as the window slides across
`list`. Even-width windows average the two central order statistics, so results
may be exact rationals (e.g. `9/2`) rather than integers. The robust median
makes it less sensitive to outliers than `MovingAverage`.
