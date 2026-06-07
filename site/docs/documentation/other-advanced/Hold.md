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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= Hold[1+1, 2+2]
Out[3]= Hold[1 + 1, 2 + 2]
```

### Notes

`Hold` has attribute `HoldAll`, so it keeps every argument unevaluated; wrap an argument in `Evaluate[...]` to force a single evaluation, and use `ReleaseHold` to strip the `Hold` and evaluate the contents.
