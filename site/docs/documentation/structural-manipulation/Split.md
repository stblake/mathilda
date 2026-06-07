# Split

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Split[list]
    splits list into runs of consecutive identical elements, returning
    a list of these runs.
Split[list, test]
    groups runs of consecutive elements ei, ej for which test[ei, ej]
    yields True.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Split[{a, a, a, b, b, a, a, c}]
Out[1]= {{a, a, a}, {b, b}, {a, a}, {c}}

In[2]:= Split[{1, 2, 3, 4, 3, 2, 1}, Less]
Out[2]= {{1, 2, 3, 4}, {3}, {2}, {1}}
```

## Implementation notes

**Algorithm.** `builtin_split` partitions the list into maximal runs of consecutive equal
elements, each run wrapped in the list's head. It scans left to right comparing each element to
its predecessor: `expr_eq` by default, or an optional two-argument test evaluated as
`test[prev, curr]`. A run boundary is placed wherever the test fails (and at the end), and the
accumulated run is emitted. Output is a list of run sublists.

- `Protected`.
- Uses `SameQ` as the default test.
- Result has the same head as the input.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
