# Unequal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs != rhs or Unequal[lhs, rhs]
    is the negation of Equal: True if lhs and rhs can be decided unequal,
    False if they can be decided equal, otherwise unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_unequal` is the value-level (not structural) negation of `Equal`. For every argument pair it first checks `expr_eq`; on failure it tries `compare_numeric` to decide equality/inequality numerically. If any pair is found equal it returns `False` immediately. It returns `True` only when *every* pair is provably unequal — either decided by `compare_numeric` or, for non-comparable values, when both sides are distinct raw data (`is_raw_data`). If some pair is neither equal nor provably unequal (e.g. symbolic), it returns `NULL` so the call stays unevaluated. Fewer than two arguments returns `True`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
