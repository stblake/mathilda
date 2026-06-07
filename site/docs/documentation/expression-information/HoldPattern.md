# HoldPattern

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
HoldPattern[expr]
    is equivalent to expr for pattern matching, but maintains expr in an unevaluated form.
HoldPattern has attributes {HoldAll, Protected}.
The left-hand sides of rules and assignments are normally evaluated before being used for matching; wrap the LHS in HoldPattern to stop that evaluation (e.g. HoldPattern[_+_] -> 0 matches any two-term sum, whereas _+_ -> 0 would match a pattern simplified by Plus before the rule is applied).
HoldPattern is removed by one layer of ReleaseHold.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= HoldPattern[_+_] -> 0
Out[1]= HoldPattern[_ + _] -> 0

In[2]:= a + b /. HoldPattern[_+_] -> 0
Out[2]= 0

In[3]:= Cases[{a -> b, c -> d}, HoldPattern[a -> _]]
Out[3]= {a -> b}

In[4]:= MatchQ[a + b, HoldPattern[_+_]]
Out[4]= True
```

## Implementation notes

- Attributes: `{HoldAll, Protected}`.
- `HoldPattern[p]` is equivalent to `p` in the pattern matcher; the matcher transparently unwraps a single-argument `HoldPattern` before matching.
- Useful on the left-hand side of rules and assignments, because those positions are normally evaluated before being used for matching. Wrapping in `HoldPattern` stops that evaluation and preserves the literal pattern shape.
- `HoldPattern` is removed by one layer of `ReleaseHold`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
