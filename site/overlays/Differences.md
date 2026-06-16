### Worked examples

```mathematica
In[1]:= Differences[{1, 4, 9, 16, 25}]
Out[1]= {3, 5, 7, 9}
```

```mathematica
In[1]:= Differences[{1, 8, 27, 64, 125, 216}, 3]
Out[1]= {6, 6, 6}
```

```mathematica
In[1]:= Differences[{a, b, c, d}, 1, 2]
Out[1]= {-a + c, -b + d}
```

```mathematica
In[1]:= Differences[{{1, 2, 3}, {4, 6, 8}, {9, 12, 15}}]
Out[1]= {{3, 4, 5}, {5, 6, 7}}
```

```mathematica
In[1]:= FoldList[Plus, 1, Differences[{1, 4, 9, 16, 25}]]
Out[1]= {1, 4, 9, 16, 25}
```

### Notes

`Differences[list]` returns the successive first differences. Iterating `n` times
(`Differences[list, n]`) is the discrete analogue of the `n`-th derivative: the
third differences of the cubes `k^3` are the constant `{6, 6, 6}` (since
`Δ^3 k^3 = 3!`). The step form `Differences[list, n, s]` differences elements `s`
apart. For a nested list, `Differences[m]` differences successive rows
element-wise. `FoldList[Plus, x, Differences[list]]` reconstructs the original
list from its first element and its differences, exhibiting `Differences` as the
inverse of the partial-sum (`Accumulate`) operation.
