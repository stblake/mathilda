# PolynomialQuotient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialQuotient[p, q, x] gives the quotient of p and q, treated as polynomials in x, with any remainder dropped.
Option: Extension -> alpha (default None) divides over Q(alpha) using the Q(alpha)[x] long-division substrate; Sqrt[c], c^(1/n), and I are recognised forms for alpha.
Extension -> None and Extension -> Automatic are accepted and currently behave as the default (no extension).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialQuotient[x^4+2x+1, x^2+1, x]
Out[1]= -1 + x^2

In[2]:= PolynomialQuotient[x^2+2x+1, x^3, x]
Out[2]= 0

In[3]:= PolynomialQuotient[x^2+x+1, 2x+1, x]
Out[3]= 1/4 + 1/2 x

In[4]:= PolynomialQuotient[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[4]= Sqrt[2] + x

In[5]:= PolynomialQuotient[x^3 - 2, x - 2^(1/3), x, Extension -> 2^(1/3)]
Out[5]= 2^(2/3) + 2^(1/3) x + x^2

In[6]:= PolynomialQuotient[x^2 + 1, x - I, x, Extension -> I]
Out[6]= I + x
```

## Implementation notes

**Algorithm.** `builtin_polynomialquotient` strips an optional `Extension -> α`/`Automatic`
(routing to `polynomialdivrem_with_extension` over `Q(α)` when applicable), then calls the
shared `poly_div_rem(p, q, x, &Q, &R)` and returns the expanded quotient `Q`, discarding `R`.
`poly_div_rem` is textbook **field long division** in `x`: it expands both operands, repeatedly
forms the next quotient term `(lc(R)/lc(q)) · x^(deg R − deg q)`, subtracts `term·q` from the
running remainder `R`, and stops when `deg R < deg q`. A fast path detects exact integer/bigint
leading-coefficient divisions (via `mpz_tdiv_qr`) so the subtraction step never needlessly
introduces rationals; otherwise quotient coefficients are formed symbolically. Coefficients are
extracted with `get_coeff_expanded` against the already-expanded divisor.

**Data structures.** `Expr` polynomial trees in expanded form; quotient and remainder are
returned through `out_Q`/`out_R` pointers.

- `Protected`.
- Default path uses polynomial long division over the field of rational functions in the coefficients.
- Option `Extension -> alpha` (default `None`) lifts $p, q$ into $\mathbb{Q}(\alpha)[x]$ and runs the Q($\alpha$)-aware long division (`qaupoly_divrem`). Recognised forms for $\alpha$: `Sqrt[c]`, `c^(1/n)`, and `I`. `Extension -> None` and `Extension -> Automatic` are accepted and currently behave as the default. The extension path requires univariate input (a single live polynomial variable other than the alpha generator); multivariate inputs fall through to the standard path.
- Threading Extension here keeps the polynomial arithmetic in the Q($\alpha$)[x] substrate and avoids the multivariate Q[$\alpha$, x] subresultant-PRS path that is exponentially slow on Sqrt[$\alpha$]-laden coefficients.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed.), Ch. 2 (division with remainder over a field).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), Ch. 2.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PolynomialQuotient[x^2 - 1, x - 1, x]
Out[1]= 1 + x
```

```mathematica
In[1]:= PolynomialQuotient[x^3 + 2 x^2 + x + 1, x + 1, x]
Out[1]= x + x^2
```

```mathematica
In[1]:= PolynomialQuotient[x^2 + 1, x, x]
Out[1]= x
```

```mathematica
In[1]:= PolynomialQuotient[x^4 - 2, x^2 - Sqrt[2], x, Extension -> Sqrt[2]]
Out[1]= Sqrt[2] + x^2
```

### Notes

`PolynomialQuotient[p, q, x]` performs long division of `p` by `q` in the
variable `x` and returns only the quotient, discarding the remainder. So
`(x^2 + 1) / x` yields `x` (dropping the `1` remainder) and `x^2 - 1` divides
exactly by `x - 1` to give `1 + x`. Pair it with `PolynomialRemainder` to
recover the full relation `p = q*quotient + remainder`. The `Extension -> alpha`
option divides over `Q(alpha)`; the default `None`/`Automatic` divides over
the rationals.
