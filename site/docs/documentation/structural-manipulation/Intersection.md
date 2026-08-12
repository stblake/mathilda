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

## Examples (4)

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

## Implementation notes

- `Flat`, `OneIdentity`, `Protected`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; the empty intersection is `{}`.
- With `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]`
  is `True`; the canonically-greatest member of each class is kept.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## See also

[Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/), [List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_intersection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intersection.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
