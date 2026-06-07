# If

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
If[cond, t]
    gives t if cond evaluates to True; gives Null otherwise.
If[cond, t, f]
    gives t if cond is True, f if False, and is left unevaluated
    if cond is neither.
If[cond, t, f, u]
    also supplies u as the result when cond is neither True nor False.
If has attribute HoldRest: only the branch chosen by cond is evaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= If[True, x, y]
Out[1]= x

In[2]:= If[a < b, 1, 0, Indeterminate]
Out[2]= Indeterminate
```

## Implementation notes

- `HoldRest`, evaluating only the chosen branch.
- Remains unevaluated if the condition is undetermined and `u` is not provided.
- `If[condition, t]` returns `Null` if `condition` evaluates to `False`.

**Attributes:** `HoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= If[3 > 2, yes, no]
Out[1]= yes
```

```mathematica
In[1]:= If[x > 0, pos, neg]
Out[1]= If[x > 0, pos, neg]
```

```mathematica
In[1]:= If[1 == 2, a, b, c]
Out[1]= b
```

```mathematica
In[1]:= If[PrimeQ[7], prime, composite]
Out[1]= prime
```

### Notes

`If[cond, t, f]` evaluates `t` when `cond` is `True` and `f` when it is `False`.
`If` has the `HoldRest` attribute, so only the selected branch is evaluated — the
other is never run. When the condition is neither `True` nor `False` (e.g. the
symbolic `x > 0`), `If` returns unevaluated rather than guessing. The optional
fourth argument `u` is returned for that indeterminate case in the four-argument
form `If[cond, t, f, u]`.
