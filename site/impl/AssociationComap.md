---
source: src/assoc_ops.c
---
**Algorithm.** `builtin_associationcomap` is the co-map (reverse of
`AssociationMap`): given a `List` of functions `{f1, f2, ...}` and a value `x`,
it builds `<|f1 -> f1[x], f2 -> f2[x], ...|>`. Each function becomes both a key
and the head of an unevaluated application `fi[x]` (built with `mk_call1`), and
the rule array is canonicalised through `assoc_from_rules`.

**Data structures.** A flat rule array handed to `assoc_from_rules`, whose
transient `KeyIndex` de-duplicates the function keys; `x` is copied once per
function into its application.

**Complexity / limits.** `O(k)` for `k` functions (plus the later cost of
evaluating each `fi[x]`). The first argument must be a `List`, otherwise
`AssociationComap::invl` is emitted through the message funnel and the call is
left unevaluated.
