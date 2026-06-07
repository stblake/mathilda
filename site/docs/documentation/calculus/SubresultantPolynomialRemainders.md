# SubresultantPolynomialRemainders

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SubresultantPolynomialRemainders[a, b, x] gives the polynomial-remainder
chain {a, b, R_2, R_3, ...} obtained by iterating pseudo-remainder over
K(coeffs)[x] until a constant or zero remainder is reached. Used by the
Lazard-Rioboo-Trager rational integration pipeline; the chain is correct
modulo content scaling, which downstream consumers strip with primitive[].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^4 + 1, 2 x^3, x]
Out[1]= {1 + x^4, 2 x^3, 2}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
