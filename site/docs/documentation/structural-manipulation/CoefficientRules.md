# CoefficientRules

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CoefficientRules[poly, {x1, x2, ...}] gives {expvec -> coeff, ...} for the monomials of poly.`**

**`CoefficientRules[poly] uses Variables[poly]; CoefficientRules[poly, vars, order] sorts by order`**

<details>
<summary>Notes</summary>

(same settings as MonomialList). Modulus -\> m reduces coefficients modulo m. Works whether or not poly is expanded. FromCoefficientRules is the inverse.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= CoefficientRules[(x + y)^3]
Out[1]= {{3, 0} -> 1, {2, 1} -> 3, {1, 2} -> 3, {0, 3} -> 1}

In[2]:= CoefficientRules[a x y^2 + b x^2 z, {x, y, z}, "DegreeReverseLexicographic"]
Out[2]= {{1, 2, 0} -> a, {2, 0, 1} -> b}
```

### Options (1)

```mathematica
In[3]:= CoefficientRules[(x + 1)^5, x, Modulus -> 2]
Out[3]= {{5} -> 1, {4} -> 1, {1} -> 1, {0} -> 1}
```

### Applications (3)

Exponent vector -> coefficient

```mathematica
In[4]:= CoefficientRules[x^2 + 2 x y + y^2, {x, y}]
Out[4]= {{2, 0} -> 1, {1, 1} -> 2, {0, 2} -> 1}
```

A single variable

```mathematica
In[5]:= CoefficientRules[3 x^2 + 1, x]
Out[5]= {{2} -> 3, {0} -> 1}
```

Order by total degree

```mathematica
In[6]:= CoefficientRules[1 + x y + x^3, {x, y}, "DegreeLexicographic"]
Out[6]= {{3, 0} -> 1, {1, 1} -> 1, {0, 0} -> 1}
```

## Implementation notes

**Algorithm.** `builtin_coefficientrules` is the `{exponent-vector -> coefficient}`
rendering of the shared monomial core (see [MonomialList](MonomialList.md)): it
calls the same `parse_and_build` to resolve the variables, reduce by an optional
`Modulus -> m`, `Expand` and split into terms, decompose each term into an
integer exponent vector and coefficient, merge like monomials, and sort by the
requested monomial order (same named orders and explicit weight matrices as
`MonomialList`, default `"Lexicographic"`). It then emits one
`Rule[{e1, ..., ek}, coeff]` per surviving monomial, wrapped in a `List`. Because
the core expands first, it works whether or not `poly` is already expanded.
`FromCoefficientRules` is the exact inverse.

**Data structures.** The same `Mono { int* exps; Expr* coeff; int64_t* key; }`
array as `MonomialList`; the only difference is the output shape — here the
exponent vector becomes a `List` of `Integer`s and is paired with the (moved,
not copied) coefficient under a `Rule`. When FLINT is available and there is no
modulus the monomials are read directly off the packed `fmpq_mpoly`/field
polynomial, bypassing the generic per-term walk.

**Complexity / limits.** Dominated by the expansion/FLINT read-off and the
`O(n log n)` monomial sort. A symbolic structural head — exact exponent vectors
and coefficient trees, so no packed/NDArray or `Compile[]` path. `Protected`.

- `Protected`.
- Returns `{expvec -> coeff, ...}`, one rule per monomial; the exponent vector lists
  the powers of `vars` in order.
- `CoefficientRules[poly]` is equivalent to `CoefficientRules[poly, Variables[poly]]`.
- Same `order` settings and `Modulus -> m` option as `MonomialList`; `vars` may be
  `All`. Works whether or not `poly` is expanded.
- `FromCoefficientRules` is the inverse.
- Note: the no-variable form uses `Variables[poly]`, which Mathilda returns in
  canonical (sorted) order — so e.g. `CoefficientRules[y + x z]` uses the variable
  order `{x, y, z}`.
- Polynomials whose coefficients live in one number field `Q(θ)`
  (`AlgebraicNumber[θ, {..}]`, all sharing one θ — the ParallelMixedTower assembly's
  representation) take a native FLINT read-off (θ→a fresh variable, group over the
  field, read each coefficient back mod θ's minimal polynomial), byte-identical to
  the generic path but without the per-term evaluator work. `MonomialList` shares it.

**Attributes:** `Protected`.

## References

**See also:** [MonomialList](../../structural-manipulation/MonomialList/), [FromCoefficientRules](../../structural-manipulation/FromCoefficientRules/)

- D. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed. (Springer, 2015), ch. 2 — monomial orderings.
- Source: [`src/poly/monomials.c`](https://github.com/stblake/mathilda/blob/main/src/poly/monomials.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)

## Notes & additional examples

### Notes

`CoefficientRules[poly, {x1, ..., xk}]` gives a sparse `{expvec -> coeff, ...}`
view of `poly`: each rule's left side is the length-`k` integer exponent vector
of a monomial and its right side the coefficient. The variables default to
`Variables[poly]`, and an optional third argument sets the monomial order — the
six named orders (`"Lexicographic"` default, `"DegreeLexicographic"`,
`"DegreeReverseLexicographic"`, and `"Negative"` variants) or an explicit weight
matrix — with `Modulus -> m` reducing coefficients modulo `m`. It works whether
or not `poly` is expanded, and `FromCoefficientRules` is the exact inverse.
