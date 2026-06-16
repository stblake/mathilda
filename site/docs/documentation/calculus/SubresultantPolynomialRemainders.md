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

**Algorithm.** `SubresultantPolynomialRemainders[a, b, x]` returns the
polynomial remainder chain `{a, b, R_2, R_3, ...}` in `K(coeffs)[x]`.
`builtin_subresultantpolynomialremainders` (src/poly/poly.c) expands both
inputs, swaps them if needed so `deg a ≥ deg b` (the documented chain
orientation), seeds the chain with `{a, b}`, and then repeatedly appends
`R = pseudo_rem(prev, cur, x)` — the **pseudo-remainder** — until the remainder
is zero or constant (degree ≤ 0). The chain array grows by doubling.

Note (per the source comment) this is a **pseudo-remainder PRS, not the
Lazard-scaled subresultant chain**: it does not divide out the subresultant
content scaling factors. This is deliberate — its sole consumer, the
Lazard–Rioboo–Trager log part of rational integration (src/calculus/intrat.c),
only uses each chain element's degree in `x` and its primitive part in the
auxiliary variable `t`, and both are invariant under content scaling, so the
simpler pseudo-remainder chain is a correct, cheaper substrate.

**Data structures.** A dynamically grown `Expr**` array of polynomial trees,
returned wrapped in a `List`. Each step relies on `pseudo_rem`,
`get_degree_poly`, `is_zero_poly`, and `expr_expand`.

**Complexity / limits.** Length of the chain is `O(deg b)` reductions; because
content is not removed, intermediate coefficients can swell relative to a true
subresultant PRS, but element degrees and `t`-primitive parts (all the consumer
needs) are exact. Requires the variable to be a symbol.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — polynomial remainder sequences / subresultants.
- W. S. Brown, J. F. Traub, "On Euclid's Algorithm and the Theory of Subresultants", JACM 18 (1971).
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

The pseudo-remainder chain starting from the two input polynomials; here it
terminates as soon as a divisor is reached:

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^4 - 1, x^2 - 1, x]
Out[1]= {-1 + x^4, -1 + x^2}
```

For coprime inputs the chain runs all the way down to a nonzero constant:

```mathematica
In[1]:= SubresultantPolynomialRemainders[x^3 - 2 x + 5, x^2 - 3, x]
Out[1]= {5 - 2 x + x^3, -3 + x^2, 5 + x, 22}
```

### Notes

`SubresultantPolynomialRemainders[a, b, x]` returns the remainder chain
`{a, b, R_2, R_3, ...}` obtained by iterating the pseudo-remainder over
`K(coeffs)[x]` until a constant or zero remainder is reached. The final nonzero
entry is the resultant up to content (a constant precisely when `a` and `b` are
coprime). The chain is correct modulo content scaling, which downstream
consumers strip with `primitive[]`; it is the workhorse of the Lazard–Rioboo–
Trager rational-integration pipeline.
