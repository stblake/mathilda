### Worked examples

```mathematica
In[1]:= LatticeReduce[{{1, 1, 1}, {-1, 0, 2}, {3, 5, 6}}]
Out[1]= {{0, 1, 0}, {1, 0, 1}, {-2, 0, 1}}
```

LLL produces a short, nearly-orthogonal basis for the same lattice. Reduction preserves the lattice determinant exactly:

```mathematica
In[1]:= Det[LatticeReduce[{{201, 37}, {1648, 297}}]]
Out[1]= -1279

In[2]:= Det[{{201, 37}, {1648, 297}}]
Out[2]= -1279
```

Integer-relation detection: appending a scaled approximation column to the identity makes a short vector reveal a relation `61*pi - 183*e + 189*phi ~ 0` (the last coordinate is tiny):

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 31415927}, {0, 1, 0, 27182818}, {0, 0, 1, 16180340}}]
Out[1]= {{61, -183, 189, 113}, {-198, 108, 203, -182}, {-235, 146, 211, 323}}
```

Reduction is exact, so rational and Gaussian-rational lattices are handled without rounding:

```mathematica
In[1]:= LatticeReduce[{{1/2, 1}, {1, 1/3}}]
Out[1]= {{1/2, -2/3}, {1, 1/3}}
```

### Notes

`LatticeReduce[m]` returns an LLL-reduced basis for the lattice spanned by the
rows of `m`. Entries may be integers, Gaussian integers, rationals, or Gaussian
rationals; arithmetic is exact (GMP rationals), so results are correct for both
machine-size and arbitrary-precision input. The reduction preserves the lattice,
its determinant, and every linear relation among the rows. The input rows must be
linearly independent.
