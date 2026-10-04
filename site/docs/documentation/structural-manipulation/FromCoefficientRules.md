# FromCoefficientRules

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromCoefficientRules[{expvec -> coeff, ...}, {x1, x2, ...}] reconstructs the polynomial`**

**`Sum[coeff * x1^e1 * x2^e2 * ..., over rules]. The exponent vectors must match the number of`**

<details>
<summary>Notes</summary>

variables. Inverse of CoefficientRules.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FromCoefficientRules[{{2, 0} -> a, {1, 1} -> b, {0, 2} -> c}, {x, y}]
Out[1]= a x^2 + b x y + c y^2
```

### Applications (3)

Rebuild the polynomial

```mathematica
In[2]:= FromCoefficientRules[{{2, 0} -> 1, {1, 1} -> 2, {0, 2} -> 1}, {x, y}]
Out[2]= x^2 + 2 x y + y^2
```

A single variable

```mathematica
In[3]:= FromCoefficientRules[{{2} -> 3, {0} -> 1}, x]
Out[3]= 1 + 3 x^2
```

Round-trips CoefficientRules

```mathematica
In[4]:= FromCoefficientRules[CoefficientRules[1 + x^3 + 7 x y^2, {x, y}], {x, y}]
Out[4]= 1 + x^3 + 7 x y^2
```

## Implementation notes

**Algorithm.** `builtin_fromcoefficientrules` reconstructs a polynomial from the
sparse form `CoefficientRules` produces — it is the inverse of that head.
`FromCoefficientRules[{expvec -> coeff, ...}, {x1, ..., xk}]` reads the variable
list (a `List` or a single bare variable, giving `k`), then walks each rule: the
left side must be a length-`k` `List` of `Integer` exponents and the right side
is the coefficient. For each rule it assembles the monomial
`coeff * x1^e1 * ... * xk^ek` with `internal_power`/`internal_times`, skipping a
variable whose exponent is `0` and dropping the `Power` wrapper when the exponent
is `1`. The monomials are summed with `internal_plus`; an empty rule list
reconstructs the `Integer` `0`.

**Data structures.** No intermediate monomial table — each rule is turned
straight into a `Times` node in a growable `Expr**` `terms` buffer, then folded
into a single `Plus`. A malformed rule (wrong head, a non-`List` or wrong-length
exponent vector, or a non-`Integer` exponent) frees the partial terms and leaves
the call unevaluated (`NULL`).

**Complexity / limits.** Linear in the number of rules times `k`; the real cost
is the evaluator's canonicalisation of the assembled `Plus` (collecting like
terms). A symbolic structural head — exact `Expr` arithmetic, so no
packed/NDArray or `Compile[]` path. `Protected`.

- `Protected`.
- Builds `Sum[coeff * var1^e1 * var2^e2 * ..., over rules]`; the exponent vectors must
  match the number of variables.
- Inverse of `CoefficientRules`: `FromCoefficientRules[CoefficientRules[poly, vars], vars]`
  is `Expand[poly]`.

**Attributes:** `Protected`.

## References

**See also:** [CoefficientRules](../../structural-manipulation/CoefficientRules/)

- Source: [`src/poly/monomials.c`](https://github.com/stblake/mathilda/blob/main/src/poly/monomials.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)

## Notes & additional examples

### Notes

`FromCoefficientRules[{expvec -> coeff, ...}, {x1, ..., xk}]` reconstructs the
polynomial `Sum[coeff x1^e1 ... xk^ek]` — the inverse of `CoefficientRules`. Each
exponent vector must have exactly `k` integer components (matching the variable
list); the variable list may also be a single bare variable. An exponent of `0`
drops the variable and an exponent of `1` drops the `Power` wrapper, and an empty
rule list reconstructs `0`. A malformed rule (wrong length or a non-integer
exponent) leaves the call unevaluated.
