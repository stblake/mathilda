# CountDistinctBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CountDistinctBy[expr, f]`**

Gives the number of distinct values of f\[e\] over the elements e of expr (the values, for an association).

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= CountDistinctBy[{1, 2, 3, 4, 5}, EvenQ]
Out[1]= 2

In[2]:= CountDistinctBy[{"apple", "avocado", "banana"}, StringTake[#, 1] &]
Out[2]= 2
```

### Applications (3)

Two classes: even and odd

```mathematica
In[3]:= CountDistinctBy[{1, 2, 3, 4, 5, 6}, EvenQ]
Out[3]= 2
```

|x| takes three distinct values

```mathematica
In[4]:= CountDistinctBy[{-1, 1, -2, 3}, Abs]
Out[4]= 3
```

```mathematica
In[5]:= CountDistinctBy[Range[10], Mod[#, 3] &]
Out[5]= 3
```

## Implementation notes

**Algorithm.** `builtin_countdistinctby` returns the number of distinct values of
`f[element]`. It shares the `count_distinct` core with `CountDistinct`, passing
the second argument as the key function: `f` is evaluated once per element
(`eval_call1`) and the result added to an `ExprSet`; the final set size is the
answer.

**Data structures.** An `ExprSet` hash set over the owned `f`-values
(`expr_hash`/`expr_eq`), plus an array that owns each `f[element]` result for the
duration of the pass. A visible `NDArray` argument is de-listed first
(`ops_delist_visible`).

**Complexity / limits.** `O(n)` evaluations of `f` plus `O(n)` hashing. Requires
a non-atomic first argument (`CountDistinctBy::normal` otherwise, via the message
funnel); `f` is arbitrary, so no packed fast path applies.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`CountDistinctBy[list, f]` counts the distinct values of `f[element]` — the number
of groups `GatherBy[list, f]` would produce, without materialising the groups. It
is a one-pass hash count: `f` is applied once per element and the distinct results
are tallied. Use it to ask how many categories a key function induces.
