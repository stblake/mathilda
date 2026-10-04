---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
  - "Garcia-Molina, Ullman and Widom, *Database Systems: The Complete Book*, 2nd ed. (Pearson, 2008), §15.4 (hash join)."
source: src/assoc_ops.c
---
**Algorithm.** `builtin_joinacross` is the relational join over two lists of
associations (rows). The join spec is a key `k`, `Key[k]`, a renaming
`k1 -> k2`, or a `List` of these (join on several columns at once). The right
rows are hashed on their join-key tuple; each left row's tuple is then probed,
and every matching right row is merged into the left row. A fourth argument
selects the join type (`"Inner"` default, `"Left"`, `"Right"`, `"Outer"`), and
the option `KeyCollisionFunction -> Left | Right | f` decides how colliding
non-join keys are kept.

**Data structures.** A `KeyIndex`/hash table over the right rows keyed by the
owned join-key tuple (`join_tuple` builds a `List` of the key values, returning
`NULL` when a join key is absent so the row never matches). Output order follows
Mathematica 15: matched pairs left-major, then unmatched left rows, then
unmatched right rows; unmatched rows carry the other side's common keys as
`Missing["Unmatched"]`, and disagreeing key sets are padded to the union with
`Missing["NotAvailable"]`.

**Complexity / limits.** `O(|left| + |right| + |output|)` — the hash index
replaces the quadratic scan of nested-loop matching. Invalid arguments route
their diagnostics through `ops_msg` (so `Quiet[]`/`Check[]` are honoured).
