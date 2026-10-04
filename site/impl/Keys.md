---
source: src/assoc.c
---
**Algorithm.** The base `builtin_keys` (`keys_or_values(res, want_keys=true)`) walks
the entries once and copies each rule's key into a fresh `List`. It also accepts a
single `Rule` (returning its key) or a bare list of rules, for Wolfram parity. The
extended forms live in a wrapper installed by `assoc_ops_init` (`ops_keys` →
`keys_values` → `kv_extract`): `Keys[assoc, f]` wraps each key as `f[k]`, and both
forms thread recursively over nested lists of associations and lists of rules.

**Data structures.** Plain `Expr` trees. One `List` of `expr_copy`'d keys is built;
no hash index is needed because the operation is a straight structural read in
insertion order.

**Complexity / limits.** O(n) in the number of entries. A list element that is not a
rule or association yields the `Keys::invrl` message (or, inside an association,
leaves the call unevaluated).
