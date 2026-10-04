# Gather

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Gather[list]`**

Gathers identical elements of list into sublists, giving {{group1}, {group2}, ...}. Sublists appear in order of the first occurrence of their element, and elements keep their input order within a sublist. Equal elements are collected from anywhere in the list, not only from adjacent runs (unlike Split). Gather\[list\] is equivalent to GatherBy\[list, Identity\].

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
Out[1]= {{1}, {7, 7}, {3, 3}, {2}, {9}}

In[2]:= Gather[{a, b, a}]
Out[2]= {{a, a}, {b}}

In[3]:= Gather[{}]
Out[3]= {}
```

### Applications (3)

```mathematica
In[4]:= Gather[{1, 7, 3, 7, 2, 3, 9}]
Out[4]= {{1}, {7, 7}, {3, 3}, {2}, {9}}

In[5]:= Gather[{a, b, a}]
Out[5]= {{a, a}, {b}}

In[6]:= Gather[{3, 1, 3, 2, 1}]
Out[6]= {{3, 3}, {1, 1}, {2}}
```

## Algorithm

--------------------------------------------------------------------------- gather.c — Gather[list], the identity case of GatherBy.

Gather partitions a list into sublists of structurally identical elements: every element appears in exactly one sublist, two elements share a sublist iff expr_eq holds between them, sublists appear in order of the first occurrence of their element, and within a sublist elements keep their input order. Unlike Split, grouping is not restricted to adjacent runs:

```text
    Gather[{1, 7, 3, 7, 2, 3, 9}]  ->  {{1}, {7, 7}, {3, 3}, {2}, {9}}
    Gather[{a, b, a}]              ->  {{a, a}, {b}}
```

The grouping itself is not reimplemented here. assoc_gather_core (assoc.c) is the same hash-indexed O(n) engine that backs GatherBy; passing a NULL key function selects the identity key, so Gather[l] and GatherBy[l, Identity] agree by construction, and the identity path skips the n Identity[x] applications that spelling it as GatherBy[l, Identity] would evaluate. --------------------------------------------------------------------------

## Implementation notes

**Algorithm.** `builtin_gather` (in `src/list/gather.c`) partitions a list into
sublists of structurally identical elements: two elements share a sublist iff
`expr_eq` holds, sublists appear in order of first occurrence, and input order is
kept within each sublist. Unlike `Split`, grouping is not restricted to adjacent
runs (`Gather[{1, 7, 3, 7, 2, 3, 9}]` → `{{1}, {7, 7}, {3, 3}, {2}, {9}}`). The
grouping is not reimplemented here — it calls `assoc_gather_core` (`assoc.c`)
with a `NULL` key function, which selects the identity key.

**Data structures.** The shared hash-indexed engine's `KeyIndex`
(open-addressing over `Expr*`) plus per-group doubling buffers; the identity path
takes a copy of each element as its own group key rather than evaluating
`Identity[x]`.

**Complexity / limits.** `O(n)` amortised. Passing `NULL` rather than
`f = Identity` makes `Gather[l] === GatherBy[l, Identity]` structural rather than
coincidental and avoids `n` redundant `Identity[x]` evaluations. Requires a
single argument; returns `NULL` otherwise.

**Attributes:** `Protected`.

## References

**See also:** [Split](../../structural-manipulation/Split/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/list/gather.c`](https://github.com/stblake/mathilda/blob/main/src/list/gather.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`Gather[list]` partitions the elements into sublists of identical elements: two
elements share a sublist when they are structurally equal, sublists appear in
order of first occurrence, and input order is kept within each sublist. Unlike
`Split`, the grouping is not limited to adjacent runs — equal elements anywhere in
the list land together. It is the identity case of `GatherBy` (`Gather[l]` is
`GatherBy[l, Identity]`) and is itself the key-less cousin of `Tally`, which
reports counts rather than the grouped elements.
