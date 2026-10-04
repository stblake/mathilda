# Complement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Complement[eall, e1, e2, ...]`**

gives the sorted list of distinct elements in eall that are not in any of the ei (set difference). The arguments must share a head, which need not be List.

**`Complement[eall, e1, ..., SameTest -> f]`**

uses f\[a, b\] to decide whether elements a and b are the same.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Complement[{a, b, c, d, e}, {a, c}, {d}]
Out[1]= {b, e}

In[2]:= Complement[f[a, b, c, d], f[c, a], f[b, b, a]]
Out[2]= f[d]

In[3]:= Complement[{b, e, d, a, b, c, d}, {b, c}]
Out[3]= {a, d, e}
```

### Options (1)

```mathematica
In[4]:= Complement[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 1.3}, SameTest -> (Floor[#1] == Floor[#2] &)]
Out[4]= {0.5, 7.1}
```

### Applications (4)

The first list minus everything in the rest

```mathematica
In[5]:= Complement[{a, b, c, d, e}, {a, c}, {d}]
Out[5]= {b, e}
```

The odd numbers in 1..10

```mathematica
In[6]:= Complement[Range[10], Range[2, 10, 2]]
Out[6]= {1, 3, 5, 7, 9}
```

Any shared head works

```mathematica
In[7]:= Complement[f[a, b, c, d], f[c, a], f[b, b, a]]
Out[7]= f[d]
```

A custom equality test

```mathematica
In[8]:= Complement[{1.1, 3.4, 0.5, 7.6, 1.9}, {1.2, 3.3}, SameTest -> (Floor[#1] == Floor[#2] &)]
Out[8]= {0.5, 7.6}
```

## Implementation notes

**Algorithm.** `builtin_complement` gives the sorted distinct elements of its
first argument that appear in *none* of the later operands (set difference),
using the head of the first argument (need not be `List`); all operands must
share that head. It is the structural twin of `Intersection` with the membership
test inverted — a candidate survives when it is absent from every later operand.
After locating a trailing `SameTest` option (`Automatic` means the default), the
default path is `O(total)`: deduplicate the first operand through the file-local
chained hash set, then drop any candidate that hits a later operand's hash set.
`SameTest -> f` switches to an `O(n^2)` path with `f[a, b] === True` as the
equivalence relation; unlike `Intersection` (which keeps the greatest member of a
class), `Complement` keeps the canonically-*smallest* member, so the first
operand is sorted ascending and the first element seen for each class is the
representative. The result is deduplicated and sorted with `expr_compare`.

**Data structures.** The same `HashTable` of chained `HashNode` buckets
(`expr_hash`/`expr_eq` keys) and an `Expr**` candidate array as `Intersection`;
the buffer fast path uses `int64_t*` arrays. Unlike `Union`/`Intersection`,
`Complement` is order-sensitive in its first argument, so it is *not* `Flat`/
`OneIdentity`.

**Complexity / limits.** `O(total)` default, `O(n^2)` with `SameTest`. A rank-1
buffer of exact integers (invisible packed `List` or explicit `NDArray[...]`)
takes `setop_packed`'s machine fast path — a sorted set-difference over `int64`
words (`acc \ b`), or direct range-indexing when bounded — so `Complement` is on
`pack.c`'s `AWARE` and `INT64_OK` lists and keeps its input's representation;
reals, a custom `SameTest`, or a non-`List` head fall to the general path. There
is **no** `Compile[]` lowering (`CompileDiagnostics` reports `Compiled -> False`).
`Protected`.

- `Protected`. Unlike `Union`/`Intersection`, `Complement` is order-sensitive in
  its first argument, so it is *not* `Flat`/`OneIdentity`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; deduplicated and sorted into standard
  order. If nothing survives the removals the result is `{}`.
- Default option `SameTest -> Automatic` (`Options[Complement]`). With
  `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]` is
  `True`; the canonically-smallest member of each class is kept.

**Attributes:** `Protected`.

## References

**See also:** [Union](../../structural-manipulation/Union/), [Intersection](../../structural-manipulation/Intersection/), [Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/), [List](../../other-advanced/List/)

- Source: [`src/list/setops.c`](https://github.com/stblake/mathilda/blob/main/src/list/setops.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`Complement[eall, e1, e2, ...]` gives the sorted distinct elements of `eall` that
appear in none of the later lists — a set difference. Unlike `Union` and
`Intersection` it is order-sensitive in its first argument, so it is not `Flat`.
All operands must share a head, which need not be `List`. The default comparison
is canonical structural equality, computed in `O(total)` with a hash set;
`SameTest -> f` uses `f[a, b] === True` as the equivalence relation (keeping the
canonically-smallest member of each class), and `SameTest -> Automatic` is the
default. A rank-1 buffer of exact integers (packed list or explicit `NDArray`)
takes a machine sorted-set-difference fast path and keeps its representation.
`Protected`.
