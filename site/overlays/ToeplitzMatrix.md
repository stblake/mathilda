### Worked examples

```mathematica
In[1]:= ToeplitzMatrix[3]
Out[1]= {{1, 2, 3}, {2, 1, 2}, {3, 2, 1}}
```

A non-symmetric Toeplitz matrix is built from an independent first column and
first row; both lists share their first entry:

```mathematica
In[1]:= ToeplitzMatrix[{a, b, c}, {a, x, y}]
Out[1]= {{a, x, y}, {b, a, x}, {c, b, a}}
```

The tridiagonal `(2, -1)` Toeplitz matrix is the standard second-difference
operator; its determinant equals `n + 1`, here `n = 5`:

```mathematica
In[1]:= Det[ToeplitzMatrix[{2, -1, 0, 0, 0}, {2, -1, 0, 0, 0}]]
Out[1]= 6
```

Its eigenvalues are the exact `2 - 2 cos(k pi/(n+1))` of the discrete Laplacian:

```mathematica
In[1]:= Eigenvalues[ToeplitzMatrix[{2, -1, 0}, {2, -1, 0}]]
Out[1]= {1/2 (4 + 2 Sqrt[2]), 2, 1/2 (4 - 2 Sqrt[2])}
```
