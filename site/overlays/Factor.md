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

```mathematica
In[1]:= Factor[x^10 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x + x^2 + x^3 + x^4) (1 - x + x^2 - x^3 + x^4)
```

```mathematica
In[1]:= Factor[x^4 + 1, Extension -> Sqrt[2]]
Out[1]= (1 - Sqrt[2] x + x^2) (1 + Sqrt[2] x + x^2)
```

```mathematica
In[1]:= Factor[x^4 - 5 x^2 + 6, Extension -> {Sqrt[2], Sqrt[3]}]
Out[1]= (Sqrt[2] + x) (Sqrt[3] + x) (-Sqrt[2] + x) (-Sqrt[3] + x)
```

### Notes

Over the integers, `Factor` returns irreducible factors and keeps an
integer content split out front; `x^2 + 1` stays irreducible because it has
no rational roots. Supplying `Extension -> I`, `Extension -> Sqrt[c]`, or
`Extension -> c^(1/n)` factors over the corresponding algebraic number field
via Trager's algorithm, so `x^2 + 1` splits over `Q(I)` and `x^2 - 2` over
`Q(Sqrt[2])`. Factors are printed in canonical term order, which places the
constant term before the leading term (`-1 + x` rather than `x - 1`).

The cyclotomic factorisation of `x^10 - 1` recovers the degree-1, degree-4
cyclotomic polynomials Φ₁, Φ₂, Φ₅, Φ₁₀ over the integers. `x^4 + 1` (the 8th
cyclotomic, irreducible over Q) splits into two real quadratics over
`Q(Sqrt[2])`. Passing a list to `Extension` factors over the compositum: over
`Q(Sqrt[2], Sqrt[3])` the biquadratic `x^4 - 5 x^2 + 6` splits completely into
four linear factors, the tower being reduced to a single primitive element via
Trager's primitive-element algorithm.
