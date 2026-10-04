# GatherBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GatherBy[list, f]`**

Gathers elements with equal f\[element\] into sublists, in first-appearance order: {{group1}, {group2}, ...}.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[1]= {{1, 3, 5}, {2, 4, 6}}

In[2]:= GatherBy[<|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>, EvenQ]
Out[2]= {<|"a" -> 1, "c" -> 3|>, <|"b" -> 2, "d" -> 4|>}
```

### Applications (3)

```mathematica
In[3]:= GatherBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[3]= {{1, 3, 5}, {2, 4, 6}}

In[4]:= GatherBy[Range[9], Mod[#, 3] &]
Out[4]= {{1, 4, 7}, {2, 5, 8}, {3, 6, 9}}
```

Group rows by their first column

```mathematica
In[5]:= GatherBy[{{1, a}, {2, b}, {1, c}}, First]
Out[5]= {{{1, a}, {1, c}}, {{2, b}}}
```

## Implementation notes

**Algorithm.** `builtin_gatherby` delegates to `assoc_gather_core`, the single
grouping engine shared with `Gather`. It evaluates `f[x]` once per element, looks
the key up in a `KeyIndex`, and appends the element to that key's buffer —
exactly like `GroupBy`, but the result is the plain list of groups with the group
keys dropped. Groups appear in first-appearance order and keep input order
within each group. Over an `Association` the entries are gathered by `f[value]`
into sub-associations (keys preserved), returned as a list.

**Data structures.** A `KeyIndex` open-addressing hash set over owned keys, with
per-group doubling `Expr**` buffers; the outer result is a `List` of `List`s (or
of sub-associations).

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Returns
`NULL` for a non-2-argument call or a non-list/association first argument. `f` is
arbitrary, so no packed fast path applies. `Gather[list]` is the identity-key
special case (`assoc_gather_core(list, NULL)`), which skips the per-element
function application entirely.

**Attributes:** `Protected`.

## References

**See also:** [GroupBy](../../data-structures/GroupBy/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`GatherBy[list, f]` gathers elements with equal `f[x]` into sublists, in
first-appearance order — like `GroupBy`, but returning the groups as a plain list
of lists with the group keys dropped. It is the natural tool for grouping records
by a field (`GatherBy[rows, First]`) or partitioning numbers by a residue class.
Over an association the entries are gathered by `f[value]` into sub-associations.
Use `GroupBy` when you want the keys kept, or `Counts`/`CountsBy` when you only
need the sizes.
