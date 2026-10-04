# Intersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Intersection[list]`**

gives the sorted list of distinct elements in list.

**`Intersection[l1, l2, ...]`**

gives the sorted list of elements common to all the li (set intersection). The li must share a head, which need not be List.

**`Intersection[l1, ..., SameTest -> f]`**

uses f\[a, b\] to decide whether elements a and b are the same.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]
Out[1]= {1, 3}

In[2]:= Intersection[f[a, b], f[c, a], f[b, b, a]]
Out[2]= f[a]

In[3]:= Intersection[Divisors[45], Divisors[78]]
Out[3]= {1, 3}
```

### Options (1)

```mathematica
In[4]:= Intersection[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 7.7, 1.3}, SameTest -> (Floor[#1] == Floor[#2] &)]
Out[4]= {1.9, 3.4, 7.6}
```

### Applications (4)

The elements common to all, sorted and distinct

```mathematica
In[5]:= Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]
Out[5]= {1, 3}
```

Two integer ranges overlap on 5..10

```mathematica
In[6]:= Intersection[Range[1, 10], Range[5, 15]]
Out[6]= {5, 6, 7, 8, 9, 10}
```

Common divisors

```mathematica
In[7]:= Intersection[Divisors[60], Divisors[45]]
Out[7]= {1, 3, 5, 15}
```

A custom equality test

```mathematica
In[8]:= Intersection[{1.1, 3.4, 0.5, 7.6, 1.9}, {1.2, 3.3, 7.7}, SameTest -> (Floor[#1] == Floor[#2] &)]
Out[8]= {1.9, 3.4, 7.6}
```

## Implementation notes

**Algorithm.** `builtin_intersection` gives the sorted list of elements common to
all operands (or, for a single argument, its sorted distinct elements), using the
head of the first argument (need not be `List`); every operand must share that
head. It first locates a trailing `SameTest -> f` option. With the default test
it is `O(total)`: deduplicate the first operand into a candidate array through the
file-local chained hash set (`ht_insert`/`ht_find`, keyed by `expr_hash`/
`expr_eq`), then for each later operand build a hash set of its members and keep
only candidates present in it. With `SameTest -> f` it switches to an `O(n^2)`
path that treats `f[a, b] === True` as the equivalence relation and keeps the
canonically-greatest member of each class. The surviving candidates are sorted
with `expr_compare` for the result.

**Data structures.** A `HashTable` of chained `HashNode` buckets for the default
path; a borrowed/`expr_copy`'d `Expr**` candidate array carried through the
intersection passes. The int64 buffer fast path uses plain `int64_t*` arrays.

**Complexity / limits.** `O(total)` with the default test, `O(n^2)` with
`SameTest`. A rank-1 buffer of exact integers (from either the invisible packed
`List` or an explicit `NDArray[...]`) takes the machine fast path
`setop_packed` — a sorted-merge intersection over `int64` words (or direct
range-indexing when the value range is bounded) — so `Intersection` is on
`pack.c`'s `AWARE` and `INT64_OK` lists and keeps whichever representation it was
given; `setop_any_nd` guards it, and reals, a custom `SameTest`, or a non-`List`
head take the general path. There is **no** `Compile[]` lowering
(`CompileDiagnostics` reports `Compiled -> False`). `Flat`, `OneIdentity`,
`Protected`.

- `Flat`, `OneIdentity`, `Protected`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; the empty intersection is `{}`.
- With `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]`
  is `True`; the canonically-greatest member of each class is kept.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/), [List](../../other-advanced/List/)

- Source: [`src/list/setops.c`](https://github.com/stblake/mathilda/blob/main/src/list/setops.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
- Tests: [`tests/test_intersection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intersection.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`Intersection[l1, l2, ...]` gives the sorted list of distinct elements common to
every operand; `Intersection[list]` alone is the sorted distinct elements of one
list. All operands must share a head, which need not be `List` (so
`Intersection[f[a, b], f[c, a]]` is `f[a]`). The default comparison is canonical
structural equality, done in `O(total)` with a hash set; `SameTest -> f` switches
to an `O(n^2)` path keeping the canonically-greatest member of each class. A
rank-1 buffer of exact integers — whether an invisible packed list or an explicit
`NDArray` — takes a machine sorted-merge fast path and keeps its representation;
reals are routed through the general path. `Flat`, `OneIdentity`, `Protected`.
