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

**Algorithm.** `builtin_polynomiallcm` mirrors `PolynomialGCD`: it strips an optional
`Extension -> α`/`Automatic` (routing to `polynomiallcm_with_extension` /
`qa_polynomiallcm_with_tower*` over `Q(α)`), force-rationalises and re-numericalises inexact
coefficients, and otherwise pre-decomposes each input with `decompose_to_bp` to handle numeric
content and shared non-numeric factors. The polynomial parts are combined pairwise via the
identity `lcm(a, b) = a·b / gcd(a, b)`, where the GCD comes from `poly_gcd_internal` (the
recursive multivariate subresultant PRS — see `PolynomialGCD`) and the exact division is done
by `exact_poly_div`/`Cancel`. Multiple arguments fold left-to-right, accumulating the running
LCM. When the exact quotient by the GCD can't be certified the code falls back to the plain
product `a·b`.

**Data structures.** `Expr` trees throughout; `BPList` for the base/power decomposition;
`QAUPoly`/`QATower` on the algebraic-extension path.

- `Protected`, `Listable`.
- Handles univariate and multivariate polynomials.
- Treats algebraic numbers (like `I`) as independent variables or constants seamlessly during complex arithmetic evaluations.
- Preserves explicit factored forms where possible.
- **Option `Extension -> alpha`** computes the LCM over `Q(alpha)` via `lcm(a, b) = a*b / PolynomialGCD[a, b, Extension -> alpha]`, returning the monic, expanded form. Same scope and fallback as `PolynomialGCD`'s extension option.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. E. Collins, "Subresultants and Reduced Polynomial Remainder Sequences", JACM 14(1), 1967.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PolynomialLCM[x^2 - 1, x^2 + 2 x + 1]
Out[1]= (-1 + x) (1 + 2 x + x^2)
```

```mathematica
In[1]:= PolynomialLCM[x^6 - 1, x^4 - 1]
Out[1]= (-1 + x^6) (1 + x^2)
```

```mathematica
In[1]:= PolynomialLCM[x^2 - 2, x^2 - Sqrt[2], Extension -> Sqrt[2]]
Out[1]= 2 Sqrt[2] - 2 x^2 - Sqrt[2] x^2 + x^4
```

### Notes

`PolynomialLCM` returns the least common multiple of its polynomial arguments,
computed as `a b / gcd(a, b)`. For `x^2 - 1 = (x-1)(x+1)` and
`(x+1)^2` the shared factor `x+1` appears only once in the LCM, so the result
factors as `(x-1)(x+1)^2`. For `x^6-1` and `x^4-1` the common part is
`x^2-1`, leaving `(x^6-1)(x^2+1)`. The `Extension -> alpha` option computes
the LCM over `Q(alpha)`; with `alpha = Sqrt[2]` the inputs `x^2-2` and
`x^2-Sqrt[2]` are coprime there, so the LCM is their full product
`(x^2-2)(x^2-Sqrt[2])`.
