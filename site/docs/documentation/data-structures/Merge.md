# Merge

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Merge[{assoc1, assoc2, ...}, f]`**

Combines associations, applying f to the list of values collected for each key (e.g. Merge\[{...}, Total\]). The list may also hold rules and lists of rules, each rule contributing one value.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Merge[{<|"a" -> 1|>, <|"a" -> 2, "b" -> 3|>}, Total]
Out[1]= <|"a" -> 3, "b" -> 3|>

In[2]:= Merge[<|"g1" -> <|"a" -> 1|>, "g2" -> <|"a" -> 2, "b" -> 3|>|>, Total]
Out[2]= <|"a" -> 3, "b" -> 3|>

In[3]:= Merge[{"a" -> 1, "b" -> 2, "a" -> 3}, Total]
Out[3]= <|"a" -> 4, "b" -> 2|>
```

### Applications (2)

Sum colliding keys

```mathematica
In[4]:= Merge[{<|a -> 1, b -> 2|>, <|a -> 10, c -> 3|>}, Total]
Out[4]= <|a -> 11, b -> 2, c -> 3|>
```

Keep every value in a list

```mathematica
In[5]:= Merge[{<|a -> 1|>, <|a -> 2|>, <|b -> 3|>}, Identity]
Out[5]= <|a -> {1, 2}, b -> {3}|>
```

## Implementation notes

**Algorithm.** `builtin_merge` makes one hash pass over all entries of all
associations. The first time a key is seen it is registered (first-seen order fixes
the output key order) with a growable value bucket; each later occurrence appends its
value. After the pass, each key's collected `List` of values is wrapped as
`f[{v1, v2, …}]` and the result is one re-canonicalised association. `Merge` also
accepts an association *of* associations, delegating to `Merge[Values[col], f]`.

**Data structures.** A single open-addressing `KeyIndex` over the distinct keys, plus
a per-key growable `Expr**` bucket (`vals`/`vcap`/`vcnt`) holding owned value copies.

**Complexity / limits.** O(N) in the total number of entries for the collection pass,
plus the cost of one `f` application per distinct key. The value bucket doubles on
growth, so building it is amortised O(1) per value.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

For each key, `f` receives the `List` of all values seen under it across the input
associations, in first-seen key order. Common combiners are `Total`, `Max`, `Mean`
and `Identity` (which keeps the raw value lists). The collection runs in one O(N)
hash pass.
