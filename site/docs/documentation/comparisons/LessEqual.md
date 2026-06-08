# LessEqual

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x <= y or LessEqual[x, y]
    yields True if x is less than or equal to y on numeric inputs,
    False if strictly greater, otherwise unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_lessequal` delegates to `evaluate_inequality(res, -1, 0)`: it walks adjacent argument pairs through `compare_numeric` (exact GMP for integer-like operands, exact cross-multiplied `long double` for rationals, `2^-46` relative tolerance for inexact reals) and accepts a pair when the sign is in `{-1, 0}` (less-than-or-equal). All pairs satisfying ≤ → `True`; any pair comparing `+1` → `False`; any non-comparable symbolic pair → NULL (left unevaluated). Implements chained `a <= b <= c`. Same implementation as `Less`/`Greater`/`GreaterEqual` with a different accepted-sign set.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2 <= 2
Out[1]= True

In[2]:= 3 <= 2
Out[2]= False

In[3]:= LessEqual[7, 7]
Out[3]= True
```

### Notes

`<=` is the operator form of `LessEqual`. It returns `True` when the arguments are equal, and `False` only when the left side is strictly greater.
