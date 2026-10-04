# KeyComplement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyComplement[{assoc1, assoc2, ...}]`**

Gives the entries of assoc1 whose keys occur in none of the other associations. Elements may also be rules or lists of rules.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeyComplement[{<|"a" -> 1, "b" -> 2, "c" -> 3|>, <|"b" -> 0|>, <|"c" -> 0|>}]
Out[1]= <|"a" -> 1|>
```

### Applications (1)

Keys of the first not in the rest

```mathematica
In[2]:= KeyComplement[{<|a -> 1, b -> 2, c -> 3|>, <|b -> 9|>}]
Out[2]= <|a -> 1, c -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_keycomplement` takes a single `List` of associations
`{a1, a2, ...}` and returns the entries of the first association `a1` whose key
appears in *none* of the later associations — the set-difference of key sets,
carrying `a1`'s values and preserving `a1`'s order. Each candidate key of `a1` is
tested against the other associations with `assoc_lookup_value`, and the entry is
kept only when every other lookup misses.

**Data structures.** The other associations are probed through their cached
`AssocIndex` (open-addressing hash), so each membership test is `O(1)` amortised;
the surviving entries are deep-copied into a fresh `Association`.

**Complexity / limits.** `O(|a1| · m)` lookups for `m` associations, each lookup
`O(1)` amortised. An empty list argument raises `KeyComplement::empt`; a
non-list-of-associations raises `KeyComplement::invar`. Both diagnostics route
through the message funnel (`ops_msg`), so they honour `Quiet[]` and `Check[]`.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`KeyComplement[{a1, a2, ...}]` returns the entries of the first association `a1`
whose key appears in *none* of the later associations — the set difference of key
sets, carrying `a1`'s values and preserving `a1`'s order. It is the association
counterpart of `Complement`, useful for finding the records present in a baseline
but absent from every comparison set. The argument must be a non-empty list of
associations (or rule lists); the sibling `KeyIntersection` keeps the common keys
instead.
