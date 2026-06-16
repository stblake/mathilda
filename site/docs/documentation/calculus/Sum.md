# Sum

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Sum[f, {i, imax}] gives the sum of f for i from 1 to imax. Sum[f, {i, imin, imax}], Sum[f, {i, imin, imax, di}] and Sum[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple iterators give nested sums. Sum[f, i] gives the indefinite sum (antidifference). Symbolic and infinite sums are evaluated in closed form via Method -> "Polynomial" | "Geometric" | "Gosper".
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sum[i^2, {i, 1, 100}]
Out[1]= 338350

In[2]:= Sum[i^2, {i, 1, n}]
Out[2]= 1/6 n (1 + n) (1 + 2 n)

In[3]:= Sum[f[i, j], {i, 1, 3}, {j, 1, i}]
Out[3]= f[1, 1] + f[2, 1] + f[2, 2] + f[3, 1] + f[3, 2] + f[3, 3]
```

```mathematica
In[1]:= Sum[i^3, {i, 1, n}]
Out[1]= 1/4 n^2 (1 + n)^2

In[2]:= Sum[i^2, i]
Out[2]= 1/6 i (-1 + i) (-1 + 2 i)
```

```mathematica
In[1]:= Sum[a^i, i]
Out[1]= a^i/(-1 + a)

In[2]:= Sum[q1^i q2^i, i]
Out[2]= (q1 q2)^i/(-1 + q1 q2)
```

```mathematica
In[1]:= Sum[k k!, k]
Out[1]= Factorial[k]
```

## Implementation notes

**Algorithm.** `Sum` is `HoldAll` so the iterator variable and bounds are not
prematurely evaluated. The dispatcher `builtin_sum` (src/sum/sum.c) strips
trailing options, rewrites multiple iterators `Sum[f, s1, ..., sk]` into nested
single-spec sums (outer bounds may depend on inner variables), and then:

- **Finite explicit expansion** — when a range resolves to a finite integer
  span, or the spec iterates an explicit list, it binds the variable and folds
  the evaluated terms with `Plus` (runaway-guarded by `SUM_MAX_FINITE_TERMS`).
- **Closed-form cascade** — for symbolic bounds, `Infinity`, or the indefinite
  form `Sum[f, i]`, it runs the context-qualified sub-algorithms in order:
  `Sum\`Polynomial` (Faulhaber-style polynomial/telescoping sums,
  src/sum/sum_polynomial.c), `Sum\`Geometric` (geometric/exponential terms,
  src/sum/sum_geometric.c), then `Sum\`Gosper` (src/sum/sum_gosper.c). Each
  stage returns the closed form or comes back unevaluated to fall through; if
  all fall through, `Sum[...]` is returned held. A `Method -> "..."` option
  selects a single strict stage.

The key stage is **Gosper's algorithm** for indefinite hypergeometric
summation. Given a term `t(i)` whose ratio `r(i) = t(i+1)/t(i)` is rational
(`Simplify` reduces factorial ratios, `Together` yields polynomials `a, b`),
it computes the **Gosper–Petkovšek normal form** `r = (a/b)·(c(i+1)/c(i))` with
`gcd(a(i), b(i+h)) = 1` for all integers `h ≥ 0`, found via the dispersion set
(`GOSPER_DISPERSION_MAX`-capped) and GCD peeling. It then solves
`a(i)·x(i+1) − b(i−1)·x(i) = c(i)` for a polynomial `x` by undetermined
coefficients (`SolveAlways`); no solution proves `t` is not Gosper-summable.
The antidifference is `F(i) = (b(i−1)/c(i))·x(i)·t(i)`; the definite finite sum
is `F(imax+1) − F(imin)`. The output stays elementary (`R(i)·t(i)`).

**Data structures.** `Expr*` trees throughout; the Gosper stage builds on
Mathilda's polynomial builtins (`Expand`, `PolynomialGCD`, `PolynomialQuotient`,
degree queries) and `SolveAlways`.

**Complexity / limits.** Finite expansion is linear in the number of terms.
Gosper's procedure is a complete decision procedure for hypergeometric
antidifferences, but only for hypergeometric terms; non-hypergeometric or
non-summable inputs fall through to the held form. No creative-telescoping
(Zeilberger) for parametric/definite hypergeometric sums.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Petkovšek, Wilf & Zeilberger, "A=B" (A K Peters, 1996).
- Graham, Knuth & Patashnik, "Concrete Mathematics", 2nd ed. (Addison-Wesley, 1994), ch. 2 & 6.
- R. W. Gosper, "Decision procedure for indefinite hypergeometric summation", Proc. Natl. Acad. Sci. USA 75 (1978).
- M. Petkovšek, H. S. Wilf, D. Zeilberger, *A = B* (A K Peters, 1996).
- Source: [`src/sum/sum_gosper.c`](https://github.com/stblake/mathilda/blob/main/src/sum/sum_gosper.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sum[k, {k, 1, 10}]
Out[1]= 55
```

```mathematica
In[1]:= Sum[k^2, {k, 1, n}]
Out[1]= 1/6 n (1 + n) (1 + 2 n)
```

Faulhaber's formula extends to high powers, here the fifth power:

```mathematica
In[1]:= Sum[k^5, {k, 1, n}]
Out[1]= 1/12 n^2 (1 + n)^2 (-1 + 2 n + 2 n^2)
```

A finite geometric sum is returned in closed form in the parameter `r`:

```mathematica
In[1]:= Sum[r^k, {k, 0, n}]
Out[1]= -1/(-1 + r) + r^(1 + n)/(-1 + r)
```

Gosper's algorithm handles the arithmetic–geometric summand `k x^k`:

```mathematica
In[1]:= Sum[k x^k, {k, 1, n}]
Out[1]= x/(1 - 2 x + x^2) + (x^(1 + n) (-1 - n - x + (1 + n) x))/(1 - 2 x + x^2)
```

Convergent infinite geometric and exponential series close in symbolic form:

```mathematica
In[1]:= Sum[1/2^k, {k, 0, Infinity}]
Out[1]= 2

In[2]:= Sum[x^k/k!, {k, 0, Infinity}]
Out[2]= E^x
```

### Notes

`Sum` evaluates numeric ranges directly and closes symbolic finite ranges in closed form through the polynomial, geometric, and Gosper (`Method`) routines — so `Sum[k^2, {k, 1, n}]` returns Faulhaber's polynomial and `Sum[k x^k, {k, 1, n}]` is summed by the Gosper backend over a symbolic upper bound. Some infinite sums are recognised: geometric series such as `Sum[1/2^k, {k, 0, Infinity}]` give `2`, and the exponential generating function `Sum[x^k/k!, {k, 0, Infinity}]` returns `E^x`. Zeta-type series such as `Sum[1/k^2, {k, 1, Infinity}]` are **not** evaluated and stay symbolic. `Sum` is `HoldAll`, so the iterator variable is not evaluated before the range is set up.
