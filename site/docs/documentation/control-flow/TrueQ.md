# TrueQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrueQ[expr] yields True if expr is True, and False otherwise.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_trueq` is a one-argument predicate (no Hold attributes, so its argument is already evaluated). It returns the symbol `True` only when the argument is exactly the interned symbol `True` (pointer equality on `SYM_True`), and `False` for everything else — so unlike a bare condition it never stays symbolic. A non-unary call returns `NULL`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TrueQ[True]
Out[1]= True

In[2]:= TrueQ[x]
Out[2]= False

In[3]:= TrueQ[1 < 2]
Out[3]= True
```

### Notes

`TrueQ` yields `True` only when its argument evaluates to the symbol `True`; anything else (including unevaluated symbolic expressions) yields `False`, making it useful as a safe guard in conditions.
