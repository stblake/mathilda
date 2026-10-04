---
source: src/funcprog.c
---
**Algorithm.** `builtin_mapindexed` wraps each selected part of `expr` in
`f[part, {position}]`, where `position` is the list of indices locating the part.
It parses an optional level spec (default `{1}`) with `parse_level_spec_strict`
and a trailing `Heads -> True` option, then rebuilds the tree bottom-up in
`mi_at_level`: each child is visited at `level + 1` carrying an `MIPath` frame
that records its index component, and a node at a level inside the spec is wrapped
with its accumulated position (`mi_position` reverses the linked frames into a
`List`). An association's parts are its *values*, positioned by `Key[k]`, and its
head is never mapped. `Rational` and `Complex` are atomic here, as they are for
`Depth` and `Level`.

A visible `NDArray` is atomic to the generic traversal, so the default level `{1}`
iterates its leading axis directly (`mi_ndarray_axis`); any other spec materialises
the array to a nested list first and repacks nothing. An empty level range (e.g.
`{3, 1}`) selects nothing and returns the expression itself, which keeps a packed
array packed.

**Data structures.** Positions are built from a stack-allocated `MIPath` linked
list (`{path, component, length}` frames) threaded through the recursion, so no
index vector is copied per node. The original subtree depth is accumulated into
`*out_depth` as the recursion unwinds — needed only by a negative level bound —
keeping the whole traversal `O(n)` rather than recomputing a depth at each node.

**Complexity / limits.** `O(n)` over the parts at a non-negative spec, with early
pass-through once past the maximum level so the common `MapIndexed[f, list]` never
walks each element's whole subtree. Negative-bound specs keep the full depth-aware
descent.
