---
source: src/assoc_ops.c
---
**Algorithm.** `Key` is an inert wrapper: `Key[k]` represents the key `k` of an
association and does not evaluate on its own (it is registered `PROTECTED` with
no builtin body in `assoc_ops.c`). It acquires meaning in three places. The
operator form `Key[k][assoc]` is handled in the evaluator (`eval.c`): it looks
`k` up and returns the value or `Missing["KeyAbsent", k]` — the curried
complement of `assoc[Key[k]]`. As a `Part` specification, `assoc[[Key[k]]]`
reaches `k` through the association (`part.c`). And wherever a key spec is read —
`Lookup`, `KeyDrop`, `JoinAcross` — `Key[k]` is unwrapped to `k`.

**Data structures.** `Key[k]` is a one-argument `EXPR_FUNCTION`; the lookups it
drives go through `assoc_lookup_value`, i.e. the association's cached
open-addressing key index.

**Complexity / limits.** The wrapper itself is `O(1)`; the lookup it triggers is
`O(1)` amortised via the key index. `Key` exists so that record pipelines such as
`GroupBy[records, Key["field"]]` and `SortBy[records, Key["field"]]` can name a
key as an ordinary operator; it is distinct from the bare key only where a key is
itself something that would otherwise evaluate.
