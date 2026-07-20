# Evaluate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Evaluate[expr]
    causes expr to be evaluated even if it appears as the argument of a function whose attributes specify that it should be held unevaluated.
Evaluate only overrides HoldFirst, HoldRest, and HoldAll attributes when it appears directly as the head of the function argument that would otherwise be held.
Evaluate does not override HoldAllComplete.
Evaluate with other than one argument reduces to Sequence: Evaluate[] gives Sequence[] and Evaluate[a, b] gives Sequence[a, b].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Evaluate[1+1]
Out[1]= 2

In[2]:= Hold[Evaluate[1+1], 2+2]
Out[2]= Hold[2, 2 + 2]

In[3]:= Hold[Evaluate[1+1], Evaluate[2+2], Evaluate[3+3]]
Out[3]= Hold[2, 4, 6]

In[4]:= Hold[f[Evaluate[1+2]]]
Out[4]= Hold[f[Evaluate[1 + 2]]]

In[5]:= Hold[Evaluate[Sin[Pi/6]]]
Out[5]= Hold[1/2]

In[6]:= Hold[Evaluate[2^10]]
Out[6]= Hold[1024]

In[7]:= HoldForm[Evaluate[1+1]]
Out[7]= 2

In[8]:= Hold[Evaluate[Length[{a,b,c}]]]
Out[8]= Hold[3]
```

## Implementation notes

`builtin_evaluate` (`src/core.c`) returns a copy of its single argument. Its real effect happens earlier: the evaluator forces `Evaluate[expr]` arguments to be evaluated even inside a `Hold*` head's held positions, so by the time the builtin runs the argument is already evaluated and it merely unwraps it.

- `Protected`.
- `Evaluate` only overrides `HoldFirst`, `HoldRest`, and `HoldAll` attributes when it appears directly as the head of the function argument that would otherwise be held.
- `Evaluate` does not override `HoldAllComplete`.
- `Evaluate` works only on the first level, directly inside a held function. It does not penetrate into deeper subexpressions.
- Outside of held contexts, `Evaluate` acts as identity.
- `Evaluate` with any number of arguments other than one reduces to `Sequence`: `Evaluate[]` gives `Sequence[]` and `Evaluate[a, b]` gives `Sequence[a, b]`. The resulting `Sequence` then splices wherever it appears, so `Hold[Evaluate[]]` gives `Hold[]` and `Hold[Evaluate[1+1, 2+2]]` gives `Hold[2, 4]`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
