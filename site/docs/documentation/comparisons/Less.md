# Less

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x \< y or Less\[x, y\] yields True if x is strictly less than y on numeric inputs, False if strictly greater or equal, otherwise unevaluated. Chained forms (x \< y \< z) become Inequality, decided pairwise. Per IEEE 754, an Indeterminate argument gives False.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 2 < 3
Out[1]= True

In[2]:= 5 < 5
Out[2]= False

In[3]:= Less[1, 10]
Out[3]= True
```

## Implementation notes

`builtin_less` delegates to the shared `evaluate_inequality(res, -1, -1)`, which walks adjacent argument pairs and calls `compare_numeric` on each. `compare_numeric` returns -1/0/+1 using exact GMP comparison for integer-like operands, exact cross-multiplied `long double` comparison for rationals, and a `2^-46` relative-tolerance comparison for inexact reals. `Less` accepts a pair only when the sign is in `{-1, -1}` (strictly less). If every pair is strictly increasing → `True`; any pair that compares `0` or `+1` → `False`; any pair that is not numerically comparable → NULL (the whole chain stays unevaluated). This implements the chained semantics of `a < b < c`. `Greater`/`LessEqual`/`GreaterEqual` are the same function with different accepted sign sets.

**Attributes:** `Protected`.

## References

**See also:** [LessEqual](../../comparisons/LessEqual/), [Greater](../../comparisons/Greater/), [GreaterEqual](../../comparisons/GreaterEqual/)

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_mateigen_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_direct.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Notes

`<` is the operator form of `Less`. A chained form such as `a < b < c` parses to `Inequality` and is decided pairwise.
