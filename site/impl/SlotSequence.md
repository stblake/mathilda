---
source: src/purefunc.c
---
`SlotSequence[n]` (`##`, `##n`) is an inert marker — `builtin_slotsequence` always returns `NULL`. During pure-`Function` application, `substitute_slots` replaces `SlotSequence[n]` with a `Sequence[...]` of copies of arguments `n` through the last (`arg_count − n + 1` of them), which then splices into the surrounding call when the `Sequence` is flattened by the evaluator. Like `Slot`, recursion stops at nested `Function` nodes. `ATTR_PROTECTED`.
