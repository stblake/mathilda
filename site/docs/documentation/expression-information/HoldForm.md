# HoldForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HoldForm[expr] prints as the expression expr, with expr maintained in an unevaluated form.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= HoldForm[1 + 1]
Out[1]= 1 + 1
```

### Applications (3)

Held, but the wrapper prints invisibly

```mathematica
In[2]:= HoldForm[1 + 1]
Out[2]= 1 + 1
```

The wrapper is really there

```mathematica
In[3]:= FullForm[HoldForm[1 + 1]]
Out[3]= HoldForm[Plus[1, 1]]
```

ReleaseHold strips it and evaluates

```mathematica
In[4]:= ReleaseHold[HoldForm[1 + 1]]
Out[4]= 2
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

## Notes & additional examples

### Notes

`HoldForm[expr]` holds `expr` unevaluated exactly like `Hold`, but the printer renders it
as just `expr` — the wrapper is invisible in output. It is the tool for displaying an
expression in unevaluated form while keeping it a genuine held expression, as `FullForm`
reveals.

`HoldForm` carries `HoldAll` and is Protected; `ReleaseHold` removes it.
