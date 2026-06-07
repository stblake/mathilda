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

**Algorithm.** `builtin_deleteduplicates` (in `src/list.c`) keeps the first occurrence of each distinct element. The default (no test) path builds a `HashTable` keyed by expression hash/equality and inserts each element only if `ht_find` reports it absent, giving expected O(N) behavior. When a custom equivalence test is supplied as the second argument it falls back to an O(N²) scan, evaluating `test[elem, kept_j]` and treating a `True` result as a duplicate.

**Data structures.** The `HashTable` from `src/list.c`'s hashing utilities (open-addressing over `Expr*` via `expr_hash`/`expr_eq`); a pre-sized `Expr**` collects the survivors before the result `List` (the original head) is built.

- `Protected`.
- Preserves the order of first occurrences.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
