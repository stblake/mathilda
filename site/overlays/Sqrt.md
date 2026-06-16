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
In[1]:= Sqrt[-9]
Out[1]= 3*I
```

Like terms combine and products of surds simplify back to integers:

```mathematica
In[1]:= Sqrt[2] + Sqrt[8]
Out[1]= 3 Sqrt[2]

In[2]:= Sqrt[12]*Sqrt[27]
Out[2]= 18
```

Nested radicals denest under `FullSimplify`:

```mathematica
In[1]:= Sqrt[3 + 2 Sqrt[2]] // FullSimplify
Out[1]= 1 + Sqrt[2]
```

Principal value of a complex root, both symbolically and to 40 digits:

```mathematica
In[1]:= Sqrt[I]
Out[1]= (1 + I)/Sqrt[2]

In[2]:= N[Sqrt[I], 40]
Out[2]= 0.70710678118654752440084436210484903928487 + 0.70710678118654752440084436210484903928487*I
```

The Puiseux series of `Sqrt[1 + x]` is delivered exactly:

```mathematica
In[1]:= Series[Sqrt[1 + x], {x, 0, 5}]
Out[1]= 1 + 1/2 x - 1/8 x^2 + 1/16 x^3 - 5/128 x^4 + 7/256 x^5 + O[x]^6
```

### Notes

`Sqrt[n]` is `Power[n, 1/2]`, so it inherits the perfect-square extraction logic:
the largest square factor is pulled out of the radical, reducing `Sqrt[50]` to
`5 Sqrt[2]` and `Sqrt[18]` to `3 Sqrt[2]`. Perfect squares collapse fully, and a
perfect-square rational like `1/4` yields the exact rational `1/2`. Negative
arguments produce a pure imaginary result, with `Sqrt[-9]` giving `3*I`. The
surd is kept in symbolic form when no square factor remains, e.g. `Sqrt[2]`.
Denesting of compound radicals like `Sqrt[3 + 2 Sqrt[2]]` is not automatic but is
recovered by `FullSimplify`. The branch cut runs along the negative real axis,
so `N[Sqrt[I], 40]` returns the upper-half-plane principal value.
