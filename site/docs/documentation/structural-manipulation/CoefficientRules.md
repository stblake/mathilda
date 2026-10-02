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

## Examples (3)

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

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)
