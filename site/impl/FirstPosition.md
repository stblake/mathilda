---
source: src/patterns.c
---
**Algorithm.** `builtin_first_position` (`src/patterns.c`) gives the position of
the first subexpression matching a pattern, in depth-first order. Rather than
re-implement traversal it **delegates to `Position` with a first-match cap** and
takes the first result, inheriting `Position`'s levelspec parsing, `Heads`
handling, association value → `Key[...]` remapping, and depth-first ordering
unchanged. A visible `NDArray` first argument is an atom to the matcher, so
`patterns_delist_visible` materialises it before matching.

The handler splits the trailing `Heads -> True|False` option from the positional
arguments `expr`, `pattern`, `default`, `levelspec`. It then builds the delegated
call: for a 2-argument association it uses the 2-argument `Position` form (the
only one that does the value→`Key` remap); otherwise it passes the supplied
levelspec — or the `{0, Infinity}` default, which includes level 0 (the whole
expression, position `{}`) — plus the match count `1` so `Position` stops at the
first hit, forwarding any `Heads` option. If `Position` returns a non-empty list
the first element is the answer; otherwise the held `default` is returned
(evaluated by the fixed-point loop only when returned), falling back to
`Missing["NotFound"]`.

**Attributes & limits.** `FirstPosition` carries `HoldRest | Protected`, which is
how `default` stays unevaluated until it is actually produced. It takes two to
four positional arguments.
