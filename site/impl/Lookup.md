---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables)."
source: src/assoc.c
---
**Algorithm.** `Lookup` is `HoldAll` so that the default stays lazy: `builtin_lookup`
evaluates only the association and the key(s), holding the default and evaluating it
once *per absent key*. `lookup_core` then dispatches on the key shape:

1. a single key → `assoc_lookup_value` (one hash probe);
2. `Key[k]` → the one literal key `k`, even when `k` is itself a list;
3. a `List` of keys → one `KeyIndex` is built over the association, then each key is
   an O(1) probe (O(n + m) overall);
4. a `List` of associations → the lookup threads, one `Lookup` per element.

A hit copies the stored value; a miss yields the evaluated default, or
`Missing["KeyAbsent", key]` when no default was given.

**Data structures.** The association's persistent `AssocIndex`, built lazily and
cached on the first single-key read (eager attachment does not survive the
fixed-point evaluation step, so it is deferred to the reader); the compiled
evaluator pre-builds it at its marshalling boundary so no worker thread mutates a
shared node. The multi-key path uses a transient open-addressing `KeyIndex`.

**Complexity / limits.** O(1) amortised per single key, O(n + m) for a list of `m`
keys. When the arguments do not form a valid lookup the call is returned with its
first two arguments evaluated.
