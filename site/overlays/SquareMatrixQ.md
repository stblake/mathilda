### Worked examples

```mathematica
In[1]:= SquareMatrixQ[{{1, 2}, {3, 4}}]
Out[1]= True
```

A rectangular matrix is rejected:

```mathematica
In[1]:= SquareMatrixQ[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= False
```

The test is purely structural, so symbolic entries are fine:

```mathematica
In[1]:= SquareMatrixQ[{{a, b}, {c, d}}]
Out[1]= True
```

It composes with matrix constructors, e.g. an order-5 identity:

```mathematica
In[1]:= SquareMatrixQ[IdentityMatrix[5]]
Out[1]= True
```

Vectors, ragged lists, and other non-square shapes return `False`:

```mathematica
In[1]:= SquareMatrixQ[{1, 2, 3}]
Out[1]= False

In[2]:= SquareMatrixQ[{{1, 2}, {3}}]
Out[2]= False
```

### Notes

`SquareMatrixQ[m]` is `True` exactly when `m` is a rank-2 tensor whose two
dimensions are equal, i.e. `Dimensions[m] == {n, n}`. It returns `False` on
non-list, ragged, rectangular, empty, or higher-rank inputs, and it works for
symbolic as well as numerical entries.
