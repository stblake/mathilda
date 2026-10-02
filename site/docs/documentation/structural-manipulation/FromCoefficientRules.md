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

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= FromCoefficientRules[{{2, 0} -> a, {1, 1} -> b, {0, 2} -> c}, {x, y}]
Out[1]= a x^2 + b x y + c y^2
```

## Implementation notes

- `Protected`.
- Builds `Sum[coeff * var1^e1 * var2^e2 * ..., over rules]`; the exponent vectors must
  match the number of variables.
- Inverse of `CoefficientRules`: `FromCoefficientRules[CoefficientRules[poly, vars], vars]`
  is `Expand[poly]`.

**Attributes:** `Protected`.

## References

**See also:** [CoefficientRules](../../structural-manipulation/CoefficientRules/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_coefficient_rules.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coefficient_rules.c)
