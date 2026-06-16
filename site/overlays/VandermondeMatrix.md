### Worked examples

```mathematica
In[1]:= VandermondeMatrix[{1, 2, 3, 4}]
Out[1]= {{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}, {1, 4, 16, 64}}
```

Symbolic nodes give the classic generic Vandermonde matrix, with successive powers along each row:

```mathematica
In[1]:= VandermondeMatrix[{a, b, c}]
Out[1]= {{1, a, a^2}, {1, b, b^2}, {1, c, c^2}}
```

Its determinant factors into the product of all pairwise node differences — the celebrated closed form `Product[xj - xi, {i < j}]`:

```mathematica
In[1]:= Factor[Det[VandermondeMatrix[{a, b, c}]]]
Out[1]= (-a + b) (-a + c) (-b + c)
```

The identity scales to higher dimensions, here giving all six pairwise differences for four nodes:

```mathematica
In[1]:= Factor[Det[VandermondeMatrix[{a, b, c, d}]]]
Out[1]= (-a + b) (-a + c) (-b + c) (-a + d) (-b + d) (-c + d)
```

### Notes

`VandermondeMatrix[{x1, ..., xn}]` gives the `n x n` matrix with entry `(i, j)` equal to `xi^(j-1)`; the two-argument form `VandermondeMatrix[nodes, k]` produces an `n x k` rectangular block. Nodes need not be numeric or distinct. The determinant vanishes exactly when two nodes coincide, which is why the Vandermonde system is invertible precisely for distinct interpolation points.
