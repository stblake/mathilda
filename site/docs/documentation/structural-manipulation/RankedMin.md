# RankedMin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RankedMin[list, n]
    gives the n-th smallest element of list.
RankedMin[list, -n]
    gives the n-th largest element of list.
    RankedMin[list, 1] is Min[list] and RankedMin[list, -1] is Max[list]. Yields a definite result when every element is a real number; +-Infinity are ordered as +-infinity. Has a packed-array fast path and is compilable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RankedMin[{12, 13, 11}, 2]
Out[1]= 12

In[2]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]
Out[2]= 3

In[3]:= RankedMax[{2.5, E, 12, 15, 485}, -2]
Out[3]= E

In[4]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]
Out[4]= Infinity
```

## Implementation notes

- `Protected`.
- `RankedMax[list, k]` is `RankedMin[list, -k]`.
- `RankedMin[list, 1]` is `Min[list]`; `RankedMin[list, -1]` is `Max[list]`.
- Yields a definite result whenever every element is a real number, including

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
