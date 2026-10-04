---
source: src/list/join.c
---
**Algorithm.** `builtin_catenate` flattens one level of a single collection: the
elements `ei` of `Catenate[{e1, e2, ...}]` must share a head and their arguments
are concatenated under it (`{{1,2},{3,4}}` → `{1,2,3,4}`). Associations take
part through their *values* (Mathematica 15): `Catenate[{<|a->1|>, <|b->2|>}]` is
`{1, 2}`, mixed list/association parts concatenate their elements and values, and
the outer collection may itself be an association. This differs from `Join`,
which concatenates several arguments.

**Data structures.** The buffer path runs first: a *packed* argument arrives as
one rank-2 row-major array (the pack gate absorbs a list of packed vectors before
`Catenate` is called), and catenating its rows is a reshape handled by
`ndstruct_catenate`; a *visible* `NDArray` list is de-listed and re-evaluated.
Otherwise the result is an ordinary `List` built by copying each part's elements
in order.

**Complexity / limits.** `O(total elements)`. Returns `NULL` (unevaluated) for a
non-list/association argument, mixed heads among the parts, or a part that is not
a list/association when an association is involved.
