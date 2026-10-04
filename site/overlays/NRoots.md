### Worked examples

```mathematica
In[1]:= NRoots[x^2 - 2 == 0, x]  (* the two real square roots of 2 *)
Out[1]= x == -1.41421 || x == 1.41421
```

```mathematica
In[1]:= NRoots[x^2 + 1 == 0, x]  (* a complex-conjugate pair *)
Out[1]= x == 0.0 - 1.0*I || x == 0.0 + 1.0*I
```

```mathematica
In[1]:= NRoots[x^3 - 1 == 0, x]  (* the three cube roots of unity *)
Out[1]= x == 1.0 || x == -0.5 - 0.866025*I || x == -0.5 + 0.866025*I
```

```mathematica
In[1]:= NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 20]  (* a digit count selects arbitrary precision *)
Out[1]= x == -1.4142135623730950488 || x == 1.4142135623730950488
```

### Notes

`NRoots[lhs == rhs, x]` returns a disjunction `x == r1 || x == r2 || ...` of
numerical roots of a polynomial equation. A root of multiplicity `k` appears as
`k` identical equations, and a single root yields a bare equation. Real and
complex coefficients are handled.

`Method -> Automatic` uses the **Aberth-Ehrlich** simultaneous iteration (all
roots refined at once from Bini's convex-hull initial placement); an Automatic
machine request first tries LAPACK companion-matrix eigenvalues (`numpy.roots`
exactly) and falls back to Aberth. `Method -> "CompanionMatrix"` forces the
companion eigenvalue route (real QR directly, complex via a real `2n x 2n`
embedding), and `Method -> "JenkinsTraub"` uses the three-stage shifted-deflation
algorithm (CPOLY, ACM TOMS 419). An exact-integer polynomial is first squarefree-
decomposed (Yun), so high multiplicities like `(x^2 - 2)^30` stay well
conditioned.

`PrecisionGoal` selects precision: `Automatic` (or `Infinity`) gives machine
`double`/`Complex`; a positive digit count runs the whole solve in MPFR complex
arithmetic. `AccuracyGoal` (default `MachinePrecision`) drives only the post-solve
residual check — a root whose Newton correction exceeds the goal triggers an
`NRoots::accgl` warning. `NRoots` requires the MPFR build.
