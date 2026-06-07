# Tally

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Tally[list] counts the number of occurrences of each distinct element in list.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Tally[{a, a, b, a, c, b, a}]
Out[1]= {{a, 4}, {b, 2}, {c, 1}}
```

## Implementation notes

**Algorithm.** `builtin_tally` counts distinct elements, returning `{element, multiplicity}`
pairs in first-occurrence order. With the default sameness test it uses a chained hash table
(`expr_hash` for bucketing, `expr_eq` for equality) for `O(n)` expected counting; with a custom
two-argument test it falls back to an `O(n²)` linear scan, evaluating `test[a, b]` per
comparison. Multiplicities are tracked in a parallel `int64_t` array.

- `Protected`.
- Returns a list of `{element, count}` pairs.
- Elements appear in the order of their first occurrence.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
