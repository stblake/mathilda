# SubtractFrom

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SubtractFrom[x, dx] or x -= dx
    subtracts dx from x and returns the new value of x.
    x -= dx is equivalent to x = x - dx.

SubtractFrom has attribute HoldFirst. The first argument x can be a
symbol or a Part expression referring to an existing value; dx may be a
number, a symbolic expression, or a list. If x has no assigned value,
SubtractFrom::rvalue is emitted and the expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_subtractfrom` (`x -= dx`) delegates to the shared `increment_core(lhs, dx, negate=true, pre=true, "SubtractFrom")`. That helper requires `lhs` to be an lvalue with an existing OwnValue (else it emits `SubtractFrom::rvalue` and returns `NULL`), reads the old value via `evaluate`, builds and evaluates `Plus[old, Times[-1, dx]]`, writes the result back through an evaluated `Set[lhs, new]` (preserving lvalue shape such as `Part[...]` via Set's `HoldFirst`), and returns the new value. Has `ATTR_HOLDFIRST` so the target symbol is not evaluated before mutation.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
