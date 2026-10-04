---
source: src/funcprog.c
---
**Algorithm.** `builtin_scan` is `Map` run for effect: it applies `f` to each
selected part of `expr`, discards every result, and returns `Null`. It parses an
optional level spec (default `{1}`) and a `Heads -> True` option, then walks the
tree in `scan_at_level` depth-first, **leaves before roots** — sub-parts (and,
with `Heads -> True`, the head) are visited before the node itself. Each visit
goes through `scan_apply`, which builds `f[part]`, evaluates it, and throws the
result away unless control flow intervenes: an in-flight `Throw` is propagated to
an enclosing `Catch`, and a `Return[ret]` is classified against the `Scan`
boundary (`eval_classify_return`) and, when consumed there, becomes the call's
return value.

An association at the default level scans its values (mirroring `Map`). A
numeric-closed body at the default level takes the `numloop_scan` fast path — it
still runs the body so a non-finite element hands control back to the interpreter.
A visible `NDArray` at the default level iterates its leading axis directly
(scalar leaf for rank 1, sub-array row otherwise); any other spec materialises to
a nested list first.

**Data structures.** `scan_apply` adopts the part copy into the `f[...]` call,
evaluates, and frees. The traversal returns `NULL` to continue and a non-`NULL`
sentinel (`Throw`/`Return`) to stop and hand back, so no result list is ever
accumulated — the point of `Scan` over `Map`.

**Complexity / limits.** `O(n)` over the parts; a non-negative upper bound stops
the descent once past the maximum level, so `Scan[f, list]` does not walk each
element's whole subtree. `get_depth` is consulted only for a negative bound.
