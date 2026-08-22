# GreaterEqual

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x \>= y or GreaterEqual\[x, y\] yields True if x is greater than or equal to y on numeric inputs, False if strictly less, otherwise unevaluated. Per IEEE 754, an Indeterminate argument gives False.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 3 >= 2
Out[1]= True

In[2]:= 2 >= 2
Out[2]= True

In[3]:= GreaterEqual[2, 5]
Out[3]= False
```

## Implementation notes

`builtin_greaterequal` delegates to `evaluate_inequality(res, 1, 0)`: it walks adjacent argument pairs through `compare_numeric` (exact GMP for integer-like operands, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals) and accepts a pair when the sign is in `{1, 0}` (greater-than-or-equal). All pairs satisfying ≥ → `True`; any pair comparing `-1` → `False`; any non-comparable symbolic pair → NULL (left unevaluated). Implements chained `a >= b >= c`. Same implementation as `Less`/`Greater`/`LessEqual` with a different accepted-sign set.

**Attributes:** `Protected`.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)
- Tests: [`tests/test_rat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_rat.c)
- Tests: [`tests/test_solve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_solve.c)

## Notes & additional examples

### Notes

`>=` is the operator form of `GreaterEqual`. Unlike `Greater`, it yields `True` when the two numeric arguments are equal.
