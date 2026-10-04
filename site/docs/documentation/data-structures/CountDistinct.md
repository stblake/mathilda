# CountDistinct

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CountDistinct[expr]`**

Gives the number of distinct elements of expr (of its values, for an association). One hash pass, O(n).

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= CountDistinct[{1, 2, 1, 3, 2}]
Out[1]= 3

In[2]:= CountDistinct[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= 2
```

### Applications (3)

```mathematica
In[3]:= CountDistinct[{a, b, a, c, b, a}]
Out[3]= 3

In[4]:= CountDistinct[{1, 1, 2, 3, 3, 3, 4}]
Out[4]= 4
```

Over an association, the distinct values

```mathematica
In[5]:= CountDistinct[<|x -> 1, y -> 1, z -> 2|>]
Out[5]= 2
```

## Implementation notes

**Algorithm.** `builtin_countdistinct` returns the number of distinct elements of
a non-atomic expression (the distinct *values*, for an association). It is the
shared `count_distinct` core with a `NULL` function: every element is added to an
`ExprSet` and the final set size is returned as an integer. One hash pass — no
sorting, no `DeleteDuplicates` round-trip.

**Data structures.** `ExprSet`, a hash set over `Expr*` keyed by `expr_hash` and
compared with `expr_eq`, sized once to the element count. A visible `NDArray`
argument is de-listed up front (`ops_delist_visible`) so the count runs over
ordinary elements.

**Complexity / limits.** `O(n)` amortised. A non-atomic first argument is
required; otherwise `CountDistinct::normal` is emitted (through the message
funnel) and the call is left unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [SameQ](../../comparisons/SameQ/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`CountDistinct[list]` is the number of distinct elements — `Length[Union[list]]`
computed in a single hash pass, without building the sorted set. Over an
association it counts the distinct *values*. It answers the "how many different
things are here?" question directly; use `CountDistinctBy` to count distinct
values of a function of each element.
