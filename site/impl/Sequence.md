---
source: src/eval.c
---
**Algorithm.** `Sequence` has no C builtin — it is a structural marker the
evaluator splices. `flatten_sequences` (`src/eval.c`) scans a function's argument
list for any child whose head is `SYM_Sequence`; if it finds one it rebuilds the
argument array, copying each `Sequence`'s contents in place of the wrapper, and
frees the emptied wrapper. `Sequence[]` therefore contributes zero arguments (it
evaporates) and `Sequence[e]` contributes one (it acts as the identity). The pass
runs at evaluation step 2.5, **before** `Flat` / `Listable` / `Orderless`, so
`f[a, Sequence[b, c], d]` is `f[a, b, c, d]` by the time attributes are consulted.

**Data structures.** The splice allocates one fresh `Expr**` argument array of the
combined length and re-homes the surviving pointers into it (`expr_copy` for the
spliced children, pointer move for the rest), then invalidates the node's cached
hash. `Sequence` is also the object produced by `BlankSequence` /
`BlankNullSequence` bindings and by `SlotSequence` (`##`), so the same splice
machinery resolves `f[a,b,c] /. f[x__] -> x`.

**Complexity / limits.** `O(n)` in the argument count of the enclosing call.
Splicing is suppressed when the enclosing head carries `SequenceHold` or
`HoldAllComplete` (checked at the step-2.5 gate), which is how an assignment or
rule can carry a `Sequence` that splices only at the eventual call site. A bare
top-level `Sequence[...]` with no enclosing function is left intact. Attributes
`Protected`.
