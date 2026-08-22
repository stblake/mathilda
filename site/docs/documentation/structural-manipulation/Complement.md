# Complement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Complement[eall, e1, e2, ...]`**

gives the sorted list of distinct elements in eall that are not in any of the ei (set difference). The arguments must share a head, which need not be List.

**`Complement[eall, e1, ..., SameTest -> f]`**

uses f\[a, b\] to decide whether elements a and b are the same.

## Examples (4)

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

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
