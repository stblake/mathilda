# SortBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SortBy[list, f]`**

Sorts the elements of list by the canonical order of f applied to each element.

**`SortBy[assoc, f]`**

Sorts an association by f applied to each value.

**`SortBy[list, f, p]`**

Compares the f values with the ordering function p.

**`SortBy[f]`**

Operator form: SortBy\[f\]\[expr\] is SortBy\[expr, f\].

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>

In[5]:= Sort[<|"a" -> 2, "b" -> 3, "c" -> 1|>, Greater]
Out[5]= <|"b" -> 3, "a" -> 2, "c" -> 1|>

In[6]:= SortBy[<|"a" -> {1, 2}, "b" -> {0, 5}, "c" -> {1, 1}|>, First, Greater]
Out[6]= <|"a" -> {1, 2}, "c" -> {1, 1}, "b" -> {0, 5}|>

In[7]:= Ordering[<|"a" -> 2, "b" -> 2, "c" -> 1|>, All, Greater]
Out[7]= {2, 1, 3}

In[8]:= ReverseSort[<|"a" -> 2, "b" -> 2, "c" -> 1|>]
Out[8]= <|"a" -> 2, "b" -> 2, "c" -> 1|>
```

### Applications (3)

Order by the first part

```mathematica
In[9]:= SortBy[{{2, 1}, {1, 2}, {3, 0}}, First]
Out[9]= {{1, 2}, {2, 1}, {3, 0}}
```

Order by magnitude

```mathematica
In[10]:= SortBy[{-3, 1, -2}, Abs]
Out[10]= {1, -2, -3}
```

Key extracted from each record

```mathematica
In[11]:= SortBy[{<|n -> 3|>, <|n -> 1|>, <|n -> 2|>}, #n &]
Out[11]= {<|n -> 1|>, <|n -> 2|>, <|n -> 3|>}
```

## Implementation notes

**Algorithm.** `builtin_sort_by` orders a collection by `f` applied to each element
(by `f[value]` for an association, the keys following). `f` is evaluated once per
element and the resulting key is cached in a `SortByPair`, then `qsort`
(`sortby_pair_cmp`) compares by that key via `expr_compare`, breaking ties by the
element itself and finally by original position — a stable canonical order. Two
variants extend this: a *list* key `{f1, f2, …}` builds the tuple `{f1[e], …}` and
sorts lexicographically; the three-argument `SortBy[coll, f, p]` first puts elements
in canonical order of their subjects, then merge-sorts by the ordering function `p`
on the `f`-values. The one-argument `SortBy[f]` returns the operator form
`Function[SortBy[#, f]]`.

**Data structures.** A `SortByPair{key, payload, pos, multi}` array; the C library
`qsort` for the two-argument form, a merge sort (`p_sort_perm`) for the `p` form.

**Complexity / limits.** O(n log n) comparisons with exactly one `f` evaluation per
element (the key is computed up front, not re-evaluated inside the comparator). The
`{f1, f2, …}` form sorts by `f1`, then `f2`, … as tie-breakers.

- With an ordering function, the sort is Mathematica 15's top-down merge sort
  (the merge keeps the left element unless `p[left, right]` is `False` or
  `-1`), so ties under a strict `p` such as `Greater` land exactly where
  Mathematica puts them. `Ordering[..., p]` and `KeySort[assoc, p]` share it.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Join](../../data-structures/Join/), [Greater](../../comparisons/Greater/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

The key `f[element]` is computed once per element and cached, so an expensive `f`
runs `n` times, not `n log n`. Ties fall back to the canonical order of the elements
and then to original position (a stable order). A list of functions `{f1, f2, …}`
sorts by `f1`, breaking ties by `f2`, and so on; the three-argument form takes an
explicit ordering function.
