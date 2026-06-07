# HoldForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HoldForm[expr] prints as the expression expr, with expr maintained in an unevaluated form.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HoldForm[1 + 1]
Out[1]= 1 + 1
```

## Implementation notes

`HoldForm` has no C handler; it is purely an evaluation/display marker. It is given `ATTR_HOLDALL | ATTR_PROTECTED` in `core_init` (`src/core.c`) so its argument stays unevaluated, and the printer renders `HoldForm[expr]` as just `expr` (the wrapper is invisible). `ReleaseHold` strips it.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
