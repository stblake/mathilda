# PolynomialQuotientRemainder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolynomialQuotientRemainder[p, q, x] returns {Quotient, Remainder}`**

<details>
<summary>Notes</summary>

such that p == Quotient\*q + Remainder, with deg(Remainder) \< deg(q) in x. Single-pass companion to PolynomialQuotient/PolynomialRemainder. Accepts an optional Extension -\> alpha rule (default None) to perform the division over Q(alpha)\[x\] rather than the rational coefficient field.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= PolynomialQuotientRemainder[x^3 + x + 1, x^2 + 1, x]
Out[1]= {x, 1}
```

### Options (1)

```mathematica
In[2]:= PolynomialQuotientRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]
Out[2]= {Sqrt[2] + x, 0}
```

### Applications (5)

```mathematica
In[1]:= PolynomialQuotientRemainder[x^2 - 1, x - 1, x]
Out[1]= {1 + x, 0}
```

```mathematica
In[1]:= PolynomialQuotientRemainder[x^5 + x + 1, x^2 + 1, x]
Out[1]= {-x + x^3, 1 + 2 x}
```

```mathematica
In[1]:= {q, r} = PolynomialQuotientRemainder[x^5 + x + 1, x^2 + 1, x];
In[2]:= Expand[q (x^2 + 1) + r]
Out[2]= 1 + x + x^5
```

```mathematica
In[1]:= PolynomialQuotientRemainder[x^4 - 2, x^2 - Sqrt[2], x, Extension -> Sqrt[2]]
Out[1]= {Sqrt[2] + x^2, 0}
```

## Implementation notes

**Algorithm.** `PolynomialQuotientRemainder[p, q, x]` returns the pair
`{quotient, remainder}` from Euclidean division of `p` by `q` in `x`, computing
both halves in one shot. `builtin_polynomialquotientremainder`
(src/poly/poly.c) calls the shared `poly_div_rem` helper, which expands both
operands, reads the leading coefficient of `q`, and runs the classical
long-division loop: at each step it forms the quotient term
`lc(R)/lc(q) · x^(deg R − deg q)`, subtracts `term·q` from the running remainder,
and repeats until `deg R < deg q`. Constant divisors are short-cut to
`{p/q, 0}`. The quotient is `Expand`-ed before being returned in a `List`.

An optional `Extension -> α` (or `Extension -> Automatic`, which autodetects the
algebraic generators of `p` and `q` via `extension_autodetect_args`) re-runs the
division inside an algebraic number field tower (`QATower`,
`polynomialdivrem_with_extension`), falling back to the plain path on lift
failure.

**Data structures.** `Expr*` polynomial trees; coefficients are exact
(`EXPR_INTEGER`/`EXPR_BIGINT`/`Rational`), with an integer/bigint fast path
(`mpz_tdiv_qr`) that avoids the `Together`/`Cancel` denominator unification when
coefficients stay integral. Extension arithmetic uses the `QATower` algebraic
number representation.

**Complexity / limits.** `O(deg p · deg q)` coefficient operations for the
machine path. The divisor must be nonzero (returns NULL otherwise) and the
variable a symbol. Multivariate inputs are handled coefficient-wise in `x`.

**Attributes:** `Protected`.

## See also

[PolynomialQuotient](../../algebra/PolynomialQuotient/), [PolynomialRemainder](../../algebra/PolynomialRemainder/)

## References

- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra* (Kluwer, 1992) — Ch. 2, polynomial division.
- Source: [`src/poly/poly.c`](https://github.com/stblake/mathilda/blob/main/src/poly/poly.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_extension_auto_builtins.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_auto_builtins.c)
- Tests: [`tests/test_intrat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intrat.c)

## Notes & additional examples

### Notes

`PolynomialQuotientRemainder[p, q, x]` performs a single long division and
returns `{quotient, remainder}` together, satisfying `p == quotient*q +
remainder` with `deg(remainder) < deg(q)`. The third example reconstructs the
dividend `x^5 + x + 1` from the returned pair, verifying the division identity.
The `Extension -> alpha` option carries out the division over `Q(alpha)[x]`;
over `Q(Sqrt[2])` the polynomial `x^4 - 2 = (x^2 - Sqrt[2])(x^2 + Sqrt[2])`
divides exactly, giving a zero remainder. This is the combined form of
`PolynomialQuotient` and `PolynomialRemainder`, computed in one pass.
