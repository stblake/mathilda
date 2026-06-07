# Hold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Hold[expr]
    maintains expr in an unevaluated form.
Hold has attribute HoldAll: its arguments are not evaluated.
Evaluate[expr] inside Hold overrides the hold and evaluates expr once.
Sequence expressions inside Hold are flattened; use HoldComplete to prevent this.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Hold` is not a C builtin. It is registered in `attr.c`'s attribute table as `{"Hold", ATTR_HOLDALL | ATTR_PROTECTED}`, i.e. the symbol `Hold` simply carries the `ATTR_HOLDALL` attribute. When the evaluator (`eval.c`) processes `Hold[args...]` it reads that attribute before evaluating arguments and therefore leaves every argument unevaluated, returning the `Hold[...]` wrapper as-is. There is no per-head logic — holding falls entirely out of the generic attribute-driven argument-evaluation step. `HoldComplete`, `HoldPattern`, and `Unevaluated` are registered the same way with `ATTR_HOLDALLCOMPLETE`/`ATTR_HOLDALL`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
