---
source: src/attr.c
---
**An attribute symbol, not a function.** `SequenceHold` is one of the evaluation
attributes, carried in the `SymbolDef` bitflags as `ATTR_SEQUENCEHOLD`
(`1 << 13`, `src/attr.h`). It names itself in the static `builtin_attrs[]` table
(`src/attr.c`) as a `Protected` symbol, and `attr_from_symbol` maps the symbol
`SequenceHold` back to the bit, so `SetAttributes[h, SequenceHold]` and
`Attributes[h]` round-trip through it.

**Effect.** The evaluator's `Sequence`-splicing pass (step 2.5, `src/eval.c`) is
guarded by `if (hold_all_complete || (attrs & ATTR_SEQUENCEHOLD))`: when the head
being evaluated carries the bit, `Sequence[...]` objects in its arguments are left
unflattened. `HoldAll` alone does **not** suppress splicing — `Hold[Sequence[a,
b]]` still gives `Hold[a, b]` — so `SequenceHold` is the specific lever that stops
it. `HoldAllComplete` implies it (the same gate). The assignment and rule heads
`Set`, `SetDelayed`, `Rule`, and `RuleDelayed` carry `SequenceHold` in the table,
which is what lets a definition store and later deliver a spliceable `Sequence`.

**Scope.** The bit only changes the splicing decision; it is independent of the
`Hold*` family, which governs argument *evaluation*. Any user head that sets it is
honoured automatically by the same gate — there is no per-head special-casing.
Attributes of the symbol itself: `Protected`.
