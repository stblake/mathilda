---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_association` normalises its arguments — any mix of
`Rule`/`RuleDelayed` nodes, `List`s of such rules, and existing associations
(which are spliced) — into one flat rule array via `collect_entries`, then hands
it to `assoc_from_rules`. That pass de-duplicates keys with the rule that *first
occurrence fixes position, last occurrence fixes value*, keeping each entry's own
head so a `RuleDelayed` stays delayed. When the input is already canonical (all
direct 2-arg rules, no splicing, no duplicate keys) the builtin returns `NULL`,
so the common `<|a -> 1, b -> 2|>` literal costs nothing to re-evaluate.

**Data structures.** De-duplication is driven by a transient `KeyIndex`: a
fixed-capacity open-addressing hash set over borrowed `Expr*` keys, sized once to
the maximum key count (power of two, load factor < 0.5) so it never rehashes and
every probe is branch-predictable. Keys are compared with `expr_eq` and hashed
with `expr_hash`. A separate persistent `AssocIndex` (`assoc_index.c`) is
attached lazily to the surviving node on the first single-key read, not at
construction — the fixed-point evaluator keeps the original literal node and
discards the rebuilt one, so an eagerly-attached index would be thrown away.

**Complexity / limits.** Construction and de-duplication are amortised `O(n)` in
the number of entries. Once the lazy index exists, `Lookup`/`KeyExistsQ` and
`assoc[key]` are `O(1)` amortised; without it a single-key probe falls back to an
`O(n)` linear `assoc_scan`. A malformed argument (anything not a rule / rule-list
/ association) leaves the node unevaluated.
