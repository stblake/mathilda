### Worked examples

```mathematica
In[1]:= DiagonalMatrix[{1, 2, 3}]
Out[1]= {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}
```

```mathematica
In[1]:= DiagonalMatrix[{a, b}, 1]
Out[1]= {{0, a, 0}, {0, 0, b}, {0, 0, 0}}
```

```mathematica
In[1]:= DiagonalMatrix[{x, y, z}, -1, 4]
Out[1]= {{0, 0, 0, 0}, {x, 0, 0, 0}, {0, y, 0, 0}, {0, 0, z, 0}}
```

```mathematica
In[1]:= DiagonalMatrix[{1, 1, 1, 1}, 2]
Out[1]= {{0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}
```

### Notes

`DiagonalMatrix[list]` places `list` on the main diagonal of an otherwise zero
square matrix. The two-argument form `DiagonalMatrix[list, k]` shifts the band to
the `k`-th diagonal — positive `k` lies above the main diagonal (a superdiagonal),
negative `k` below it (a subdiagonal) — and the matrix grows just large enough to
hold that band, so `DiagonalMatrix[{a, b}, 1]` is `3 x 3`. The three-argument form
`DiagonalMatrix[list, k, n]` pads with zeros to force an explicit `n x n` size, as
in the `2 x` Jordan-style superdiagonal block of ones above.
