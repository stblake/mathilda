# PolynomialLCM

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialLCM[poly1, poly2, ...] gives the least common multiple of the polynomials.
Option Extension -> alpha computes the LCM over Q(alpha) via
lcm(a, b) = a*b / PolynomialGCD[a, b, Extension -> alpha].
Default Extension -> None computes over the rationals.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialLCM[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]
Out[1]= (2 + x) (3 + x) (4 + x) (1 + x)^2

In[2]:= PolynomialLCM[x^4-4, x^4+4 x^2+4]
Out[2]= (-2 + x^2) (4 + 4 x^2 + x^4)

In[3]:= PolynomialLCM[x - Sqrt[2], x + Sqrt[2], Extension -> Sqrt[2]]
Out[3]= -2 + x^2
```

## Implementation notes

- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Preserves explicit factored forms where possible.
- **Option `Extension -> alpha`** computes the LCM over `Q(alpha)` via `lcm(a, b) = a*b / PolynomialGCD[a, b, Extension -> alpha]`, returning the monic, expanded form. Same scope and fallback as `PolynomialGCD`'s extension option.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
