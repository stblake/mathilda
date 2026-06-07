# Decrement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Decrement[x] or x--
    decreases the value of x by 1, returning the old value of x.

Decrement has attribute HoldFirst. In Decrement[x], x can be a symbol
or a Part expression referring to an existing value (e.g. list[[2]]--).
If x has no assigned value, Decrement::rvalue is emitted and the
expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_decrement` (`src/core.c`) implements `x--` via the shared `increment_core` helper with a delta of `1`, negate=true, pre=false. `increment_core` requires the target be a symbol with an existing OwnValue (else `Decrement::rvalue`), evaluates the old value, forms and evaluates `Plus[old, Times[-1, 1]]`, and writes it back through an evaluated `Set`. Because pre=false it returns the *old* value (post-decrement). `Decrement` is `ATTR_HOLDFIRST`. The pre-form `--x` is the separate `builtin_predecrement` (pre=true).

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
