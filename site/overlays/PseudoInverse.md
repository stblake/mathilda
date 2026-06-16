### Worked examples

For a non-singular square matrix the pseudoinverse coincides with the ordinary inverse:

```mathematica
In[1]:= PseudoInverse[{{1, 2}, {3, 4}}] == Inverse[{{1, 2}, {3, 4}}]
Out[1]= True
```

On a wide (full-row-rank) integer matrix the right inverse is computed exactly in rational arithmetic:

```mathematica
In[1]:= PseudoInverse[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{-17/18, 4/9}, {-1/9, 1/9}, {13/18, -2/9}}
```

The Moore-Penrose defining identity `A . A^+ . A == A` holds exactly:

```mathematica
In[1]:= A = {{1, 2, 3}, {4, 5, 6}}; A . PseudoInverse[A] . A
Out[1]= {{1, 2, 3}, {4, 5, 6}}
```

It handles rank-deficient inputs gracefully; the all-ones matrix has rank 1 and a rank-1 pseudoinverse:

```mathematica
In[1]:= PseudoInverse[{{1, 1}, {1, 1}}]
Out[1]= {{1/4, 1/4}, {1/4, 1/4}}
```

Inexact input is rationalised, solved exactly, and returned at the input precision:

```mathematica
In[1]:= PseudoInverse[{{1., 2.}, {3., 4.}, {5., 6.}}]
Out[1]= {{-1.33333, -0.333333, 0.666667}, {1.08333, 0.333333, -0.416667}}
```

### Notes

`PseudoInverse[m]` computes the Moore-Penrose pseudoinverse via a full-rank
decomposition `m = B . C`, so it never forms an ill-conditioned normal-equation
solve. For exact (integer / rational / complex) input the result is exact, and
for non-singular square `m` it reduces to `Inverse[m]`. Rank-deficient and
rectangular matrices are handled without error; the defining Penrose relations
(here `A . A^+ . A == A`) hold identically. Inexact Real / MPFR input is
rationalised at its working precision, solved in exact arithmetic, and
numericalised back, giving the inexact-in / inexact-out contract.
`Tolerance -> t` drops singular values below `t` times the largest.
