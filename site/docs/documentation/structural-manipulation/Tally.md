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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Tally[{a, b, a, c, b, a}]
Out[1]= {{a, 3}, {b, 2}, {c, 1}}
```

```mathematica
In[1]:= Tally[Table[Mod[n^2, 5], {n, 0, 20}]]
Out[1]= {{0, 5}, {1, 8}, {4, 8}}
```

```mathematica
In[1]:= Tally[Table[GCD[n, 12], {n, 1, 12}]]
Out[1]= {{1, 4}, {2, 2}, {3, 2}, {4, 2}, {6, 1}, {12, 1}}
```

### Notes

`Tally[list]` returns `{element, count}` pairs for each distinct element, in the
order of first appearance. It is a compact way to read off the distribution of a
computed sequence — for example, the multiplicities of the quadratic residues
mod 5 (`{0, 1, 4}` appearing 5, 8, and 8 times among `n = 0 .. 20`), or the
divisor structure of `GCD[n, 12]`.
