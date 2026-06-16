# PolynomialRemainder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialRemainder[p, q, x] gives the remainder from dividing p by q, treated as polynomials in x.
Option: Extension -> alpha (default None) computes the remainder over Q(alpha); see PolynomialQuotient for the recognised alpha forms.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialRemainder[x^4+2x+1, x^2+1, x]
Out[1]= 2 + 2 x

In[2]:= PolynomialRemainder[x^3, a x+b, x]
Out[2]= -b^3/a^3

In[3]:= PolynomialRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[3]= 0

In[4]:= PolynomialRemainder[x^2 - 3, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[4]= -1
```

## Implementation notes

**Algorithm.** `builtin_polynomialremainder` is the companion to `PolynomialQuotient`: after
stripping an optional `Extension -> α`/`Automatic` option (routing to
`polynomialdivrem_with_extension` over `Q(α)`), it calls the shared
`poly_div_rem(p, q, x, &Q, &R)` and returns the remainder `R`, discarding the quotient. The
division is field long division in `x` (see `PolynomialQuotient`): subtract
`(lc(R)/lc(q))·x^(deg R − deg q)·q` from the running remainder until `deg R < deg q`, with a
fast exact-integer-division path for pure integer/bigint leading coefficients. The
`PolynomialQuotientRemainder` builtin exposes both outputs `{Q, R}` from a single call.

**Data structures.** Expanded `Expr` polynomial trees; coefficients extracted via
`get_coeff_expanded`.

- `Protected`.
- The degree of the result in $x$ is guaranteed to be smaller than the degree of $q$.
- If the dividend is a multiple of the divisor, then the remainder is zero.
- Option `Extension -> alpha`: see `PolynomialQuotient` for the recognised alpha forms and the fall-through rules.

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
In[1]:= PolynomialRemainder[x^2 - 1, x - 1, x]
Out[1]= 0
```

```mathematica
In[1]:= PolynomialRemainder[x^3 + 2 x^2 + x + 1, x + 1, x]
Out[1]= 1
```

```mathematica
In[1]:= PolynomialRemainder[x^2 + 1, x, x]
Out[1]= 1
```

```mathematica
In[1]:= PolynomialRemainder[x^5 + x + 1, x^2 + 1, x]
Out[1]= 1 + 2 x
```

### Notes

`PolynomialRemainder[p, q, x]` returns the remainder left after dividing `p`
by `q` in `x`; its degree is always strictly less than that of `q`. A zero
remainder, as for `x^2 - 1` divided by `x - 1`, certifies that `q` divides `p`
exactly. Evaluating `p` at a root of a linear divisor recovers the remainder:
`x^3 + 2x^2 + x + 1` at `x = -1` is `1`, matching the output. The companion
`PolynomialQuotient` returns the quotient part; together they satisfy
`p = q*quotient + remainder`. The `Extension -> alpha` option computes the
remainder over `Q(alpha)`.
