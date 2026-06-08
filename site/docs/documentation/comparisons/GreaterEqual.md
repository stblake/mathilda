# GreaterEqual

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x >= y or GreaterEqual[x, y]
    yields True if x is greater than or equal to y on numeric inputs,
    False if strictly less, otherwise unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_greaterequal` delegates to `evaluate_inequality(res, 1, 0)`: it walks adjacent argument pairs through `compare_numeric` (exact GMP for integer-like operands, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals) and accepts a pair when the sign is in `{1, 0}` (greater-than-or-equal). All pairs satisfying ≥ → `True`; any pair comparing `-1` → `False`; any non-comparable symbolic pair → NULL (left unevaluated). Implements chained `a >= b >= c`. Same implementation as `Less`/`Greater`/`LessEqual` with a different accepted-sign set.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 3 >= 2
Out[1]= True

In[2]:= 2 >= 2
Out[2]= True

In[3]:= GreaterEqual[2, 5]
Out[3]= False
```

### Notes

`>=` is the operator form of `GreaterEqual`. Unlike `Greater`, it yields `True` when the two numeric arguments are equal.
