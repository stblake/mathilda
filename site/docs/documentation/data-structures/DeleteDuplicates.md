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
In[1]:= DeleteDuplicates[<|"a" -> 1, "b" -> 1, "c" -> 2, "d" -> 2, "e" -> 3|>]
Out[1]= <|"a" -> 1, "c" -> 2, "e" -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_deleteduplicates` (in `src/list.c`) keeps the first occurrence of each distinct element. The default (no test) path builds a `HashTable` keyed by expression hash/equality and inserts each element only if `ht_find` reports it absent, giving expected O(N) behavior. When a custom equivalence test is supplied as the second argument it falls back to an O(N²) scan, evaluating `test[elem, kept_j]` and treating a `True` result as a duplicate.

**Data structures.** The `HashTable` from `src/list.c`'s hashing utilities (open-addressing over `Expr*` via `expr_hash`/`expr_eq`); a pre-sized `Expr**` collects the survivors before the result `List` (the original head) is built.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= DeleteDuplicates[{1,2,1,3,2,4}]
Out[1]= {1, 2, 3, 4}
```

```mathematica
In[1]:= DeleteDuplicates[{a,b,a,c,b,d,a}]
Out[1]= {a, b, c, d}
```

```mathematica
In[1]:= DeleteDuplicates[{1,2,3,4,5,6}, Mod[#1,3]==Mod[#2,3]&]
Out[1]= {1, 2, 3}
```

```mathematica
In[1]:= DeleteDuplicates[{1,-1,2,-2,3,-3}, Abs[#1]==Abs[#2]&]
Out[1]= {1, 2, 3}
```

### Notes

`DeleteDuplicates[list]` keeps the first occurrence of each distinct element and preserves the original order. With a two-argument equivalence test `DeleteDuplicates[list, test]` two elements are treated as duplicates when `test[a, b]` is `True`, so you can deduplicate by an arbitrary relation rather than literal equality — picking one representative per residue class modulo 3, or one per absolute value, in the examples above.
