# DifferenceDelta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DifferenceDelta[f, i] gives the forward difference (f /. i -> i+1) - f, the discrete analogue of D. It is the left inverse of indefinite Sum.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= DifferenceDelta[i^2, i]
Out[1]= 1 + 2 i

In[2]:= DifferenceDelta[Sum[k k!, k], k]
Out[2]= -Factorial[k] + Factorial[1 + k]
```

### Applications (4)

```mathematica
In[3]:= DifferenceDelta[n^2, n]
Out[3]= 1 + 2 n

In[4]:= DifferenceDelta[f[n], n]
Out[4]= -f[n] + f[1 + n]

In[5]:= DifferenceDelta[n^3, n]
Out[5]= 1 + 3 n + 3 n^2

In[6]:= DifferenceDelta[Binomial[n, k], n]
Out[6]= -Binomial[n, k] + Binomial[1 + n, k]
```

## Algorithm

sum_gosper.c -- Sum`Gosper: Gosper's indefinite hypergeometric summation, plus the DifferenceDelta forward-difference operator.

Given a hypergeometric term t(i) (one whose term ratio t(i+1)/t(i) is a rational function of i), Gosper's algorithm finds a hypergeometric

```text
antidifference F with F(i+1)-F(i) = t(i), or proves none exists.  The output
```

has the shape F = R(i) t(i) with R rational, so no new special functions are needed.

```text
  1. r(i) = t(i+1)/t(i); require it rational (Simplify reduces factorial
     ratios, then Together gives num/den polynomials a, b).
  2. Gosper-Petkovsek normal form r = (a/b)(c(i+1)/c(i)) with
     gcd(a(i), b(i+h)) = 1 for all integers h >= 0, via the dispersion set
     (h with deg gcd(a(i), b(i+h)) > 0) and gcd peeling.
  3. Solve a(i) x(i+1) - b(i-1) x(i) = c(i) for a polynomial x by undetermined
     coefficients (SolveAlways).  No solution => t is not Gosper-summable.
  4. Antidifference F(i) = (b(i-1)/c(i)) x(i) t(i).

  Sum`Gosper[f, i]              -> F(i)                 (indefinite)
  Sum`Gosper[f, i, imin, imax]  -> F(imax+1) - F(imin)  (definite, finite)
```

## Implementation notes

**Algorithm.** `DifferenceDelta[f, i]` computes the forward difference
`(f /. i -> i+1) - f`, the discrete analogue of `D` and the left inverse of
indefinite `Sum`. `builtin_differencedelta` (src/sum/sum_gosper.c) requires the
second argument to be a symbol, substitutes `i -> i+1` into `f` (`shift_var`,
implemented via `ReplaceAll`), subtracts the original `f`, and returns the
`Expand`-ed result. Returns NULL (unevaluated) unless the variable is a symbol.

**Data structures.** Plain `Expr*` tree manipulation built on the existing
`ReplaceAll`, subtraction, and `Expand` builtins (`shift_var`, `sum_sub`,
`sum_eval`). No closed-form machinery — it is a thin structural operator that
lives alongside Gosper's summation because the two are inverse operations.

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [Sum](../../calculus/Sum/)

- Source: [`src/sum/sum_gosper.c`](https://github.com/stblake/mathilda/blob/main/src/sum/sum_gosper.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_sum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum.c)

## Notes & additional examples

### Notes

`DifferenceDelta[f, n]` is the forward difference operator `Δ f = (f /. n -> n+1) - f`,
the discrete analogue of the derivative `D`. On `n^2` it returns `2 n + 1`, the
discrete counterpart of `2 n`; applied to `n^3` it gives `3 n^2 + 3 n + 1`. For an
unknown function head it expands literally to `f[n+1] - f[n]`. The fourth example
is Pascal's rule in disguise: `Binomial[n+1, k] - Binomial[n, k] = Binomial[n, k-1]`.
`DifferenceDelta` is the left inverse of indefinite `Sum`, mirroring the way `D`
inverts the indefinite integral.
