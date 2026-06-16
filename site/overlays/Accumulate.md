### Worked examples

```mathematica
In[1]:= Accumulate[{1, 2, 3, 4, 5}]
Out[1]= {1, 3, 6, 10, 15}
```

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}
```

```mathematica
In[1]:= Accumulate[Table[1/k, {k, 1, 5}]]
Out[1]= {1, 3/2, 11/6, 25/12, 137/60}
```

```mathematica
In[1]:= Accumulate[{{1, 2}, {3, 4}, {5, 6}}]
Out[1]= {{1, 2}, {4, 6}, {9, 12}}
```

### Notes

`Accumulate[list]` gives the running (prefix) sums, equivalent to
`FoldList[Plus, list]`, and the result has the same length as the input. It
works on exact rationals — `Accumulate[Table[1/k, {k, 1, 5}]]` returns the
partial sums of the harmonic series — as well as symbolic terms. For a matrix
(a list of rows), accumulation is performed row-wise, so column totals build up
down the matrix. The head of the input is preserved.
