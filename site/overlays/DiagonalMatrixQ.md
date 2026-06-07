### Worked examples

```mathematica
In[1]:= DiagonalMatrixQ[{{1, 0}, {0, 2}}]
Out[1]= True

In[2]:= DiagonalMatrixQ[{{1, 2}, {0, 3}}]
Out[2]= False
```

### Notes

Without a `Tolerance` option the test is structural: only literal numeric zeros off the main diagonal count as zero. Returns `False` on non-matrix, ragged, or higher-rank inputs.
