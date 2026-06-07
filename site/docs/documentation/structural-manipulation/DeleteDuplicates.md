# DeleteDuplicates

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DeleteDuplicates[list]
    returns list with duplicate elements removed, keeping the first
    occurrence of each element and preserving the original order.
DeleteDuplicates[list, test]
    treats two elements as duplicates when test[a, b] yields True.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DeleteDuplicates[{a, a, b, a, c, b, a}]
Out[1]= {a, b, c}
```

## Implementation notes

- `Protected`.
- Preserves the order of first occurrences.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
