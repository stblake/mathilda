---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `AssociationMap` has two forms. The key form
`AssociationMap[f, {k1, ...}]` (`builtin_associationmap`, `assoc.c`) builds
`<|k1 -> f[k1], ...|>`, each `f[k]` left for the evaluator. The association form
`AssociationMap[f, assoc]` (`ops_associationmap`, `assoc_ops.c`) applies `f` to
each entry *as a rule* `k -> v` and splices the result: `f` may return a rule, a
list of rules, an association, or `Nothing` (which contributes no entry). Both
forms canonicalise the collected rules through `assoc_from_rules`.

**Data structures.** The rule array feeds `assoc_from_rules`, whose transient
`KeyIndex` open-addressing hash set de-duplicates keys (first position, last
value). The association form keeps the raw applications in an unevaluated
`Association` if `f` ever returns something that is not a rule / rule-list /
association.

**Complexity / limits.** `O(n)` applications of `f` plus `O(n)` hashing. The key
form requires a `List` of keys; an invalid application in the association form
raises `AssociationMap::invrlf` through the message funnel.
