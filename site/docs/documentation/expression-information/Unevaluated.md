# Unevaluated

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Unevaluated[expr]
    represents the unevaluated form of expr when it appears as the argument to a function.
f[Unevaluated[expr]] effectively works by temporarily holding that argument, then evaluating f[expr] with the wrapper removed.
The wrapper is NOT removed when the argument is held (e.g. when f has HoldAll, HoldFirst, or HoldRest applies) or when f has attribute HoldAllComplete.
Sequence expressions directly inside Unevaluated are NOT flattened: Length[Unevaluated[Sequence[a, b]]] gives 2.
Unevaluated has attributes {HoldAllComplete, Protected}, so its own argument is itself protected from evaluation, Sequence flattening, and wrapper stripping.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Length[Unevaluated[Plus[5, 6, 7, 8]]]
Out[1]= 4

In[2]:= Length[Unevaluated[Sequence[a, b]]]
Out[2]= 2

In[3]:= Hold[Evaluate[Unevaluated[1+2]]]
Out[3]= Hold[Unevaluated[1 + 2]]

In[4]:= SetAttributes[f, HoldAll]; f[Unevaluated[1+2]]
Out[4]= f[Unevaluated[1 + 2]]

In[5]:= HoldComplete[Unevaluated[1+2]]
Out[5]= HoldComplete[Unevaluated[1 + 2]]

In[6]:= Attributes[Unevaluated]
Out[6]= {HoldAllComplete, Protected}
```

## Implementation notes

`Unevaluated` has no C handler; the attribute table (`src/attr.c`) gives it `ATTR_HOLDALLCOMPLETE | ATTR_PROTECTED`. The evaluator (`src/eval.c`) special-cases it: in a non-held argument position `f[Unevaluated[expr]]` passes `expr` itself (unevaluated) to `f`, stripping the wrapper. In held positions and under `HoldAllComplete` heads the wrapper is preserved.

- Attributes: `{HoldAllComplete, Protected}`.
- `f[Unevaluated[expr]]` passes `expr` to `f` as if `f` temporarily held that single argument; the `Unevaluated` wrapper is then stripped before `f`'s body runs, effectively yielding `f[expr]` with `expr` unevaluated.
- The wrapper is **not** stripped when the enclosing function holds the argument (e.g. `f` has `HoldAll`, or `HoldFirst`/`HoldRest` applies to that position).
- The wrapper is **not** stripped when the enclosing function has `HoldAllComplete`.
- Stripping happens **after** `Sequence` flattening, so a `Sequence` directly inside `Unevaluated` survives into the argument slot (`Length[Unevaluated[Sequence[a, b]]]` gives `2`).
- Nested `Unevaluated` wrappers are stripped one layer per evaluation step.
- As a top-level expression, `Unevaluated[expr]` evaluates to itself (because of `HoldAllComplete`).

**Attributes:** `HoldAllComplete`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
