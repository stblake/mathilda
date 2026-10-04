---
source: src/assoc_ops.c
---
**Algorithm.** `builtin_missingq` returns `True` exactly when its single argument has
head `Missing` — `head_is(arg, SYM_Missing)` — for any arity (`Missing[]`,
`Missing["reason"]`, `Missing["reason", data]`). A call with other than one argument
emits `MissingQ::argx` and is left unevaluated.

**Data structures.** None beyond the head comparison against the interned
`SYM_Missing`.

**Complexity / limits.** O(1). It is purely a head test, so it says nothing about
*why* the value is missing; pair it with `DeleteMissing` to remove such entries, or
with `Lookup`'s default to replace them.
