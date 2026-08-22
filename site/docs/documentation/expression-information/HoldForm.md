# HoldForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HoldForm[expr] prints as the expression expr, with expr maintained in an unevaluated form.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= HoldForm[1 + 1]
Out[1]= 1 + 1
```

## Implementation notes

`HoldForm` has no C handler; it is purely an evaluation/display marker. It is given `ATTR_HOLDALL | ATTR_PROTECTED` in `core_init` (`src/core.c`) so its argument stays unevaluated, and the printer renders `HoldForm[expr]` as just `expr` (the wrapper is invisible). `ReleaseHold` strips it.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_evaluate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_evaluate.c)
- Tests: [`tests/test_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric.c)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)
- Tests: [`tests/test_releasehold.c`](https://github.com/stblake/mathilda/blob/main/tests/test_releasehold.c)
