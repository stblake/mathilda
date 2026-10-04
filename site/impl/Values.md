---
source: src/assoc.c
---
**Algorithm.** The base `builtin_values` (`keys_or_values(res, want_keys=false)`)
walks the entries once and copies each rule's right-hand side into a fresh `List`,
also accepting a single `Rule` or a bare list of rules. The extended forms come from
the `assoc_ops_init` wrapper (`ops_values` → `keys_values` → `kv_extract`):
`Values[assoc, f]` wraps each value as `f[v]`, and both forms thread recursively over
nested lists of associations and lists of rules.

**Data structures.** Plain `Expr` trees; one `List` of `expr_copy`'d values in
insertion order. No hash index — it is a structural read.

**Complexity / limits.** O(n). A non-rule, non-association list element triggers
`Values::invrl` (or leaves the call unevaluated inside an association).
