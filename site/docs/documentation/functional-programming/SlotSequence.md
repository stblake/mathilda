# SlotSequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

## or SlotSequence\[n\] represents arguments from the n-th onward.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= (Plus[##] &)[1, 2, 3]
Out[1]= 6

In[2]:= (f[##2] &)[a, b, c]
Out[2]= f[b, c]

In[3]:= FullForm[##]
Out[3]= SlotSequence[1]
```

## Implementation notes

`SlotSequence[n]` (`##`, `##n`) is an inert marker — `builtin_slotsequence` always returns `NULL`. During pure-`Function` application, `substitute_slots` replaces `SlotSequence[n]` with a `Sequence[...]` of copies of arguments `n` through the last (`arg_count − n + 1` of them), which then splices into the surrounding call when the `Sequence` is flattened by the evaluator. Like `Slot`, recursion stops at nested `Function` nodes. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

- Source: [`src/purefunc.c`](https://github.com/stblake/mathilda/blob/main/src/purefunc.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_purefunc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_purefunc.c)

## Notes & additional examples

### Notes

`##` (`SlotSequence[1]`) splices *all* arguments of the enclosing pure function into the surrounding expression; `##n` (Out[2]) splices the arguments from the n-th onward.
