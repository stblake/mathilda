# Intersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Intersection[list]
    gives the sorted list of distinct elements in list.
Intersection[l1, l2, ...]
    gives the sorted list of elements common to all the li (set
    intersection). The li must share a head, which need not be List.
Intersection[l1, ..., SameTest -> f]
    uses f[a, b] to decide whether elements a and b are the same.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]
Out[1]= {1, 3}

In[2]:= Intersection[f[a, b], f[c, a], f[b, b, a]]
Out[2]= f[a]

In[3]:= Intersection[Divisors[45], Divisors[78]]
Out[3]= {1, 3}
```

## Implementation notes

- `Flat`, `OneIdentity`, `Protected`.
- All expressions must have the same head, which need not be `List`.
- Result has the same head as the inputs; the empty intersection is `{}`.
- With `SameTest -> f`, elements `a`, `b` are treated as equal when `f[a, b]`

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
