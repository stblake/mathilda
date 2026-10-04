# KeySortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeySortBy[assoc, f]`**

Sorts an association by f applied to each key (stable).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]
Out[1]= <|"a" -> 2, "cc" -> 3, "bbb" -> 1|>
```

### Applications (1)

Order keys by length

```mathematica
In[2]:= KeySortBy[<|"bbb" -> 1, "a" -> 2, "cc" -> 3|>, StringLength]
Out[2]= <|"a" -> 2, "cc" -> 3, "bbb" -> 1|>
```

## Implementation notes

**Algorithm.** `builtin_keysortby` computes `f[key]` once per entry, storing the
owned `f`-values in a parallel array, then runs a **stable insertion sort** keyed by
those values (compared with `expr_compare`). Stability is deliberate: entries whose
`f`-values tie keep their association order. The entries are rebuilt into a fresh
`Association` in the sorted order.

**Data structures.** Two parallel `Expr**` arrays — the entry copies and their
`f`-values — reordered together by the insertion sort.

**Complexity / limits.** O(n²) comparisons in the worst case (insertion sort, chosen
for a cheap stable order) with one `f` evaluation per entry, which is well suited to
the modest sizes associations usually reach. `KeySort` sorts by the keys themselves;
`SortBy` is the list/association analogue that sorts by a function of the values.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

Entries are ordered by `f` applied to each key, with ties keeping the association's
original order (the sort is stable). `f` is evaluated once per key. Contrast
`KeySort` (sort by the keys directly) and `SortBy` (sort by a function of the
values).
