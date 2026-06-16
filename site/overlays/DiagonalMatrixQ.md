### Worked examples

```mathematica
In[1]:= DiagonalMatrixQ[{{1, 0}, {0, 2}}]
Out[1]= True

In[2]:= DiagonalMatrixQ[{{1, 2}, {0, 3}}]
Out[2]= False
```

```mathematica
In[1]:= DiagonalMatrixQ[DiagonalMatrix[{a, b, c}]]
Out[1]= True
```

```mathematica
In[1]:= DiagonalMatrixQ[{{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}, 1]
Out[1]= True
```

```mathematica
In[1]:= DiagonalMatrixQ[{{0.0, 1.0*10^-15}, {0, 0.0}}, Tolerance -> 10^-10]
Out[1]= True
```

### Notes

Without a `Tolerance` option the test is structural: only literal numeric zeros off the main diagonal count as zero. Returns `False` on non-matrix, ragged, or higher-rank inputs.
