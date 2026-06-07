# SlotSequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
## or SlotSequence[n] represents arguments from the n-th onward.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`SlotSequence[n]` (`##`, `##n`) is an inert marker — `builtin_slotsequence` always returns `NULL`. During pure-`Function` application, `substitute_slots` replaces `SlotSequence[n]` with a `Sequence[...]` of copies of arguments `n` through the last (`arg_count − n + 1` of them), which then splices into the surrounding call when the `Sequence` is flattened by the evaluator. Like `Slot`, recursion stops at nested `Function` nodes. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/purefunc.c`](https://github.com/stblake/mathilda/blob/main/src/purefunc.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= (Plus[##] &)[1, 2, 3]
Out[1]= 6

In[2]:= (f[##2] &)[a, b, c]
Out[2]= f[b, c]

In[3]:= FullForm[##]
Out[3]= SlotSequence[1]
```

### Notes

`##` (`SlotSequence[1]`) splices *all* arguments of the enclosing pure function into the surrounding expression; `##n` (Out[2]) splices the arguments from the n-th onward.
