### Worked examples

```mathematica
In[1]:= RowReduce[{{1, 2}, {3, 4}}]
Out[1]= {{1, 0}, {0, 1}}
```

```mathematica
In[1]:= RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}, {0, 0, 0}}
```

A 4x4 Vandermonde matrix is nonsingular, so it reduces to the identity:

```mathematica
In[1]:= RowReduce[{{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}, {1, 4, 16, 64}}]
Out[1]= {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
```

Reducing an augmented matrix `[A | b]` reads off the unique solution of `A x = b`
in the last column:

```mathematica
In[1]:= RowReduce[{{2, 1, -1, 8}, {-3, -1, 2, -11}, {-2, 1, 2, -3}}]
Out[1]= {{1, 0, 0, 2}, {0, 1, 0, 3}, {0, 0, 1, -1}}
```

A symbolic 2x2 matrix reduces to the identity, certifying generic invertibility:

```mathematica
In[1]:= RowReduce[{{a, b}, {c, d}}]
Out[1]= {{1, 0}, {0, 1}}
```

### Notes

`RowReduce[m]` returns the reduced row-echelon form of `m`. The default method is
a Bareiss-like fraction-free Gauss-Jordan elimination (`"DivisionFreeRowReduction"`),
which avoids intermediate-coefficient swell; pass `Method -> "OneStepRowReduction"`
for classical division-per-pivot, or `Method -> "CofactorExpansion"` to certify
invertibility via the determinant. The pivot columns of the result identify a basis
for the row space, and a zero bottom row (as in the 3x3 magic-style example) signals
linear dependence.
