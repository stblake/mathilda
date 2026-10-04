---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables)."
source: src/assoc.c
---
**Algorithm.** `builtin_merge` makes one hash pass over all entries of all
associations. The first time a key is seen it is registered (first-seen order fixes
the output key order) with a growable value bucket; each later occurrence appends its
value. After the pass, each key's collected `List` of values is wrapped as
`f[{v1, v2, …}]` and the result is one re-canonicalised association. `Merge` also
accepts an association *of* associations, delegating to `Merge[Values[col], f]`.

**Data structures.** A single open-addressing `KeyIndex` over the distinct keys, plus
a per-key growable `Expr**` bucket (`vals`/`vcap`/`vcnt`) holding owned value copies.

**Complexity / limits.** O(N) in the total number of entries for the collection pass,
plus the cost of one `f` application per distinct key. The value bucket doubles on
growth, so building it is amortised O(1) per value.
