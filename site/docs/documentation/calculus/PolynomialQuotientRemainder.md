# PolynomialQuotientRemainder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PolynomialQuotientRemainder[p, q, x] returns {Quotient, Remainder}
such that p == Quotient*q + Remainder, with deg(Remainder) < deg(q)
in x. Single-pass companion to PolynomialQuotient/PolynomialRemainder.
Accepts an optional Extension -> alpha rule (default None) to perform
the division over Q(alpha)[x] rather than the rational coefficient field.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PolynomialQuotientRemainder[x^3 + x + 1, x^2 + 1, x]
Out[1]= {x, 1}

In[2]:= PolynomialQuotientRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[2]= {Sqrt[2] + x, 0}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
