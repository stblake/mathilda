# Less

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x < y or Less[x, y]
    yields True if x is strictly less than y on numeric inputs, False
    if strictly greater or equal, otherwise unevaluated.
Chained forms (x < y < z) become Inequality, decided pairwise.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_less` delegates to the shared `evaluate_inequality(res, -1, -1)`, which walks adjacent argument pairs and calls `compare_numeric` on each. `compare_numeric` returns -1/0/+1 using exact GMP comparison for integer-like operands, exact cross-multiplied `long double` comparison for rationals, and a `2^-46` relative-tolerance comparison for inexact reals. `Less` accepts a pair only when the sign is in `{-1, -1}` (strictly less). If every pair is strictly increasing → `True`; any pair that compares `0` or `+1` → `False`; any pair that is not numerically comparable → NULL (the whole chain stays unevaluated). This implements the chained semantics of `a < b < c`. `Greater`/`LessEqual`/`GreaterEqual` are the same function with different accepted sign sets.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2 < 3
Out[1]= True

In[2]:= 5 < 5
Out[2]= False

In[3]:= Less[1, 10]
Out[3]= True
```

### Notes

`<` is the operator form of `Less`. A chained form such as `a < b < c` parses to `Inequality` and is decided pairwise.
