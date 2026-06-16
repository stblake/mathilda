### Worked examples

```mathematica
In[1]:= UpperTriangularMatrixQ[{{1, 2}, {0, 3}}]
Out[1]= True

In[2]:= UpperTriangularMatrixQ[{{1, 0}, {2, 3}}]
Out[2]= False
```

The test is purely structural, so it works on a symbolically generated triangular matrix just as well as on numeric data:

```mathematica
In[1]:= UpperTriangularMatrixQ[Table[If[j >= i, a[i, j], 0], {i, 4}, {j, 4}]]
Out[1]= True
```

The two-argument form shifts the reference diagonal. A strictly upper-triangular matrix (zero main diagonal) passes the `k = 1` test:

```mathematica
In[1]:= UpperTriangularMatrixQ[{{0, 1, 2}, {0, 0, 3}, {0, 0, 0}}, 1]
Out[1]= True
```

For approximate matrices, a `Tolerance` lets tiny round-off entries below the diagonal count as zero:

```mathematica
In[1]:= UpperTriangularMatrixQ[{{1.0, 2}, {1*10^-15, 3}}, Tolerance -> 10^-10]
Out[1]= True
```

### Notes

A matrix is upper triangular when every entry below the main diagonal is zero. The two-argument form `UpperTriangularMatrixQ[m, k]` tests against the k-th diagonal (positive `k` for superdiagonals, negative for subdiagonals), and a `Tolerance` option relaxes the zero test for approximate matrices. Rectangular matrices are supported.
