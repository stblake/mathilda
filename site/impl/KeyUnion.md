---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables)."
source: src/assoc.c
---
**Algorithm.** `builtin_keyunion` equalises a list of associations to a common key
set so they can be processed row-wise. First it collects the union of all keys in
first-appearance order with one hash pass over every entry. Then each association is
rebuilt against a hash index of *its own* keys: for each union key present it copies
the stored value, and for one absent it fills `Missing["KeyAbsent", key]`. The result
is the list of equalised associations. The `assoc_ops_init` wrapper (`ops_keyunion`)
adds the rule/rule-list element forms and `KeyUnion[{…}, f]`, which fills an absent
key `k` with `f[k]` instead.

**Data structures.** A union `KeyIndex` over all distinct keys (keys borrowed), and a
fresh per-association `KeyIndex` so each union key is resolved in one probe rather
than by a membership scan.

**Complexity / limits.** O(total entries) for the union, then O(m · u) probes for `m`
associations and `u` union keys — linear rather than the O(keys · entries) of
repeated scans.
