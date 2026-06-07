# Slot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
# or Slot[n] represents the n-th argument of a pure function.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Slot[n]` (`#`, `#n`) is an inert marker — `builtin_slot` always returns `NULL` so the node stays unevaluated on its own. Substitution happens only when a pure `Function` is applied: `substitute_slots` walks the function body and replaces each `Slot[n]` with a copy of the `n`-th argument (when `1 ≤ n ≤ arg_count`), stopping recursion at any nested `Function` so inner pure functions are not captured by the outer slots. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/purefunc.c`](https://github.com/stblake/mathilda/blob/main/src/purefunc.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
