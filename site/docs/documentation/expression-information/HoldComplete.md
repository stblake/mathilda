# HoldComplete

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HoldComplete[expr]`**

shields expr completely from evaluation.

<details>
<summary>Notes</summary>

HoldComplete has attribute HoldAllComplete: it prevents argument evaluation, Sequence flattening, Unevaluated stripping, and Evaluate from firing. Substitution (via ReplaceAll, etc.) still happens inside HoldComplete. HoldComplete is removed by one level of ReleaseHold.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

## See also

[Sequence](../../expression-information/Sequence/), [Unevaluated](../../expression-information/Unevaluated/), [Evaluate](../../expression-information/Evaluate/), [HoldAllComplete](../../expression-information/HoldAllComplete/), [ReplaceAll](../../assignment-and-rules/ReplaceAll/), [Replace](../../assignment-and-rules/Replace/), [ReleaseHold](../../expression-information/ReleaseHold/)

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_releasehold.c`](https://github.com/stblake/mathilda/blob/main/tests/test_releasehold.c)
- Tests: [`tests/test_sequence.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sequence.c)
- Tests: [`tests/test_unevaluated.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unevaluated.c)
