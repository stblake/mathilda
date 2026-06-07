---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on square-free factorization and radical simplification."
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", on integer square roots."
---
### Worked examples

```mathematica
In[1]:= Sqrt[50]
Out[1]= 5 Sqrt[2]
```

```mathematica
In[1]:= Sqrt[18]
Out[1]= 3 Sqrt[2]
```

```mathematica
In[1]:= Sqrt[-9]
Out[1]= 3*I
```

```mathematica
In[1]:= Sqrt[1/4]
Out[1]= 1/2
```

### Notes

`Sqrt[n]` is `Power[n, 1/2]`, so it inherits the perfect-square extraction logic:
the largest square factor is pulled out of the radical, reducing `Sqrt[50]` to
`5 Sqrt[2]` and `Sqrt[18]` to `3 Sqrt[2]`. Perfect squares collapse fully, and a
perfect-square rational like `1/4` yields the exact rational `1/2`. Negative
arguments produce a pure imaginary result, with `Sqrt[-9]` giving `3*I`. The
surd is kept in symbolic form when no square factor remains, e.g. `Sqrt[2]`.
