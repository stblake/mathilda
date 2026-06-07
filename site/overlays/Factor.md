---
status: Stable
references:
  - "B. M. Trager, \"Algebraic factoring and rational function integration\", SYMSAC 1976 — the norm / sqfr_norm / alg_factor approach used for the Extension path."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 8 (polynomial factorization)."
---
### Worked examples

```mathematica
In[1]:= Factor[x^4 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x^2)
```

```mathematica
In[1]:= Factor[6 x^2 + 7 x + 2]
Out[1]= (1 + 2 x) (2 + 3 x)
```

```mathematica
In[1]:= Factor[x^2 + 1, Extension -> I]
Out[1]= (-I + x) (I + x)
```

```mathematica
In[1]:= Factor[x^2 - 2, Extension -> Sqrt[2]]
Out[1]= (Sqrt[2] + x) (-Sqrt[2] + x)
```

### Notes

Over the integers, `Factor` returns irreducible factors and keeps an
integer content split out front; `x^2 + 1` stays irreducible because it has
no rational roots. Supplying `Extension -> I`, `Extension -> Sqrt[c]`, or
`Extension -> c^(1/n)` factors over the corresponding algebraic number field
via Trager's algorithm, so `x^2 + 1` splits over `Q(I)` and `x^2 - 2` over
`Q(Sqrt[2])`. Factors are printed in canonical term order, which places the
constant term before the leading term (`-1 + x` rather than `x - 1`).
