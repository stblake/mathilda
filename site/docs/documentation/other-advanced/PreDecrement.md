# PreDecrement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PreDecrement[x] or --x
    decreases the value of x by 1, returning the new value of x.
    --x is equivalent to x = x - 1.

PreDecrement has attribute HoldFirst. In PreDecrement[x], x can be a
symbol or a Part expression referring to an existing value.
If x has no assigned value, PreDecrement::rvalue is emitted and the
expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_predecrement` (`--x`) calls the shared `increment_core(lhs, 1, negate=true, pre=true, "PreDecrement")`. It writes back `Plus[old, Times[-1, 1]]` via an evaluated `Set` and, because `pre=true`, returns the new (decremented) value. Requires an existing OwnValue on the target, else emits `PreDecrement::rvalue` and returns `NULL`. Carries `ATTR_HOLDFIRST | ATTR_PROTECTED`.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
