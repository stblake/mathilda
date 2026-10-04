# KeyUnion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyUnion[{assoc1, assoc2, ...}]`**

Gives the list of associations padded to the union of all their keys; a key absent from an association is filled with Missing\["KeyAbsent", key\].

**`KeyUnion[{assoc1, assoc2, ...}, f]`**

Fills each absent key k with f\[k\] instead.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeyUnion[{<|"a" -> 1|>, <|"b" -> 2|>}, 0 &]
Out[1]= {<|"a" -> 1, "b" -> 0|>, <|"a" -> 0, "b" -> 2|>}
```

### Applications (1)

Pad both to keys a, b, c

```mathematica
In[2]:= KeyUnion[{<|a -> 1, b -> 2|>, <|b -> 3, c -> 4|>}]
Out[2]= {<|a -> 1, b -> 2, c -> Missing["KeyAbsent", c]|>, <|a -> Missing["KeyAbsent", a], b -> 3, c -> 4|>}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [Missing](../../data-structures/Missing/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

Every returned association carries the same keys, in first-appearance order across
the inputs, so the list is ready for tabular / row-wise processing. A key missing
from a given association is filled with `Missing["KeyAbsent", key]`; the two-argument
`KeyUnion[{…}, f]` fills it with `f[key]` instead.
