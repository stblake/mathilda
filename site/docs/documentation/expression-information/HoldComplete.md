# HoldComplete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HoldComplete[expr]
    shields expr completely from evaluation.
HoldComplete has attribute HoldAllComplete: it prevents argument evaluation, Sequence flattening, Unevaluated stripping, and Evaluate from firing.
Substitution (via ReplaceAll, etc.) still happens inside HoldComplete.
HoldComplete is removed by one level of ReleaseHold.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Attributes[HoldComplete]
Out[1]= {HoldAllComplete, Protected}

In[2]:= HoldComplete[1+1, Evaluate[1+2], Sequence[3, 4]]
Out[2]= HoldComplete[1 + 1, Evaluate[1 + 2], Sequence[3, 4]]

In[3]:= HoldComplete[Sequence[a, b]]
Out[3]= HoldComplete[Sequence[a, b]]

In[4]:= HoldComplete[f[1+2]] /. f[x_] :> g[x]
Out[4]= HoldComplete[g[1 + 2]]

In[5]:= ReleaseHold[HoldComplete[Sequence[1, 2]]]
Out[5]= Sequence[1, 2]
```

## Implementation notes

`HoldComplete` has no C handler; the attribute table in `src/attr.c` gives it `ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED`, so the evaluator suppresses all argument evaluation and the upvalue/Sequence/Unevaluated-stripping machinery as well. `ReleaseHold` removes the wrapper.

- Attributes: `{HoldAllComplete, Protected}`.
- `HoldComplete` prevents argument evaluation, `Sequence` flattening inside its own arguments, `Unevaluated` wrapper stripping, and `Evaluate` from firing. `Evaluate` cannot override `HoldAllComplete`.
- Structural substitution (via `ReplaceAll`, `Replace`, `ReplacePart`, etc.) still descends into `HoldComplete` because substitution is not part of evaluation.
- `HoldComplete` is removed by one layer of `ReleaseHold`.
- `HoldComplete` is a milder form of `Unevaluated` at top level: `HoldComplete` always keeps the wrapper, while `Unevaluated` is typically stripped by the enclosing function.

**Attributes:** `HoldAllComplete`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
