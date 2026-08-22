# Greater

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x \> y or Greater\[x, y\] yields True if x is strictly greater than y on numeric inputs, False if strictly less or equal, otherwise unevaluated. Chained forms become Inequality. Per IEEE 754, an Indeterminate argument gives False.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= 3 > 2
Out[1]= True

In[2]:= 5 > 5
Out[2]= False

In[3]:= Greater[5, 3]
Out[3]= True

In[4]:= x > 2
Out[4]= x > 2
```

## Implementation notes

`builtin_greater` delegates to the shared `evaluate_inequality(res, 1, 1)`, which walks adjacent argument pairs and calls `compare_numeric` (exact GMP for integer-like, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals; returns -1/0/+1). A pair is accepted only when the sign is `+1` (strictly greater). All pairs strictly decreasing → `True`; any `0`/`-1` pair → `False`; any non-comparable (symbolic) pair → NULL (chain left unevaluated). This is the chained `a > b > c` semantics. Shares its implementation with `Less`/`LessEqual`/`GreaterEqual` via different accepted-sign arguments.

**Attributes:** `Protected`.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)
- Tests: [`tests/test_readwrite.c`](https://github.com/stblake/mathilda/blob/main/tests/test_readwrite.c)

## Notes & additional examples

### Notes

`>` is the operator form of `Greater`. On purely numeric arguments it decides to `True` or `False`; if an argument is symbolic and the relation cannot be settled, the expression is returned unevaluated.
