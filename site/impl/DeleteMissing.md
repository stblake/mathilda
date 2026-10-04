---
source: src/assoc_ops.c
---
**Algorithm.** `DeleteMissing` drops `Missing[...]` elements. The one-argument
form `DeleteMissing[expr]` (the original `builtin_delete_missing`, in
`patterns.c`) rewrites to `DeleteCases[expr, _Missing]`. The two- and
three-argument forms are handled by `ops_deletemissing` (`assoc_ops.c`):
`DeleteMissing[expr, n]` recurses with `dm_rec` to levels `1..n` (association
values counting one level down), and `DeleteMissing[expr, n, d]` deletes the
level-`1..n` elements that *contain* a `Missing[...]` at depth `d` or less.

**Data structures.** `dm_rec` rebuilds lists and associations bottom-up, dropping
a child when `has_missing_within` finds a qualifying `Missing`; an association
child is re-wrapped with `assoc_entry_with_value` so keys survive. A packed-list
argument is returned unchanged — a machine buffer holds no `Missing`.

**Complexity / limits.** Linear in the tree size visited (bounded by the level
spec `n`). An argument that is not a list, association, or packed list raises
`DeleteMissing::invrp`; a bad level or depth raises `DeleteMissing::arg2` — both
through the message funnel.
