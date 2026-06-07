### Worked examples

```mathematica
In[1]:= UpperTriangularMatrixQ[{{1, 2}, {0, 3}}]
Out[1]= True

In[2]:= UpperTriangularMatrixQ[{{1, 0}, {2, 3}}]
Out[2]= False
```

### Notes

A matrix is upper triangular when every entry below the main diagonal is zero. The two-argument form `UpperTriangularMatrixQ[m, k]` tests against the k-th diagonal, and a `Tolerance` option relaxes the zero test for approximate matrices.
