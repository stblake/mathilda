# AssociationThread

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationThread[{k...}, {v...}]  |  AssociationThread[keys -> values]`**

Builds \<|k1 -\> v1, ...|\> from parallel key and value lists.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= AssociationThread[{"a", "b"}, {1, 2}]
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= AssociationThread[{"a", "b"} -> {1, 2}]
Out[2]= <|"a" -> 1, "b" -> 2|>
```

### Applications (2)

```mathematica
In[3]:= AssociationThread[{a, b, c}, {1, 2, 3}]
Out[3]= <|a -> 1, b -> 2, c -> 3|>
```

The keys -> values rule form

```mathematica
In[4]:= AssociationThread[{a, b} -> {1, 2}]
Out[4]= <|a -> 1, b -> 2|>
```

## Implementation notes

**Algorithm.** `builtin_associationthread` accepts either the two-argument form
`AssociationThread[keys, values]` or the single rule form
`AssociationThread[keys -> values]`, both requiring two equal-length `List`s. It
zips them into `Rule[key_i, value_i]` nodes and canonicalises the result through
`assoc_from_rules`, so a repeated key keeps its first position and takes its last
value (the same de-duplication rule as `Association`).

**Data structures.** The temporary rule array is handed to `assoc_from_rules`,
which drives a transient `KeyIndex` open-addressing hash set over the keys to
de-duplicate in one pass; the rules are deep-copied into the result and the
temporaries freed.

**Complexity / limits.** `O(n)` amortised. Returns `NULL` (unevaluated) unless
both arguments are `List`s of the same length, matching Wolfram's shape
requirement.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`AssociationThread[keys, values]` zips two equal-length lists into an
association, pairing `keys[[i]]` with `values[[i]]`. The single-argument rule form
`AssociationThread[keys -> values]` is equivalent. It is the inverse of taking
`Keys` and `Values` apart, and the quickest way to build an association from two
parallel lists. Duplicate keys collapse with the usual *first position, last
value* rule; the two lists must have the same length.
