### Worked examples

```mathematica
In[1]:= HilbertMatrix[3]
Out[1]= {{1, 1/2, 1/3}, {1/2, 1/3, 1/4}, {1/3, 1/4, 1/5}}
```

Rectangular Hilbert matrices are available too:

```mathematica
In[1]:= HilbertMatrix[{2, 4}]
Out[1]= {{1, 1/2, 1/3, 1/4}, {1/2, 1/3, 1/4, 1/5}}
```

The Hilbert matrix is the textbook ill-conditioned matrix, yet its exact
rational inverse is integer-valued:

```mathematica
In[1]:= Inverse[HilbertMatrix[3]]
Out[1]= {{9, -36, 30}, {-36, 192, -180}, {30, -180, 180}}
```

Its determinant is a tiny but exact rational (these are reciprocals of the
Hilbert determinants), illustrating the near-singularity:

```mathematica
In[1]:= Det[HilbertMatrix[5]]
Out[1]= 1/266716800000
```

Because entries are kept exact, eigenvalues come out in closed form:

```mathematica
In[1]:= Eigenvalues[HilbertMatrix[2]]
Out[1]= {1/24 (16 + 4 Sqrt[13]), 1/24 (16 - 4 Sqrt[13])}
```

### Notes

`HilbertMatrix[n]` has entries `1/(i + j - 1)` and is the canonical ill-conditioned test matrix. Entries are exact `Rational`s, so `Det`, `Inverse`, and `Eigenvalues` return exact answers; request `WorkingPrecision -> MachinePrecision` (or a digit count) for an inexact version.
