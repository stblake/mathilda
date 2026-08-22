# Union

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Union[list]`**

gives the sorted list of distinct elements in list.

**`Union[l1, l2, ...]`**

gives the sorted list of distinct elements appearing in any of the input lists (set union).

<details>
<summary>Notes</summary>

Comparison is by canonical structural equality.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Union[{1, 2, 1, 3, 6, 2, 2}]
Out[1]= {1, 2, 3, 6}

In[2]:= Union[{a, b, a, c}, {d, a, e, b}, {c, a}]
Out[2]= {a, b, c, d, e}
```

### Applications (3)

```mathematica
In[3]:= Union[{3, 1, 2, 1, 3}]
Out[3]= {1, 2, 3}

In[4]:= Union[{1, 2, 3}, {2, 3, 4}, {5}]
Out[4]= {1, 2, 3, 4, 5}

In[5]:= Union[{x, Sin[y], x, 1, Sin[y]}]
Out[5]= {1, x, Sin[y]}
```

## Options & behaviour

`Union`, `Intersection`, `Complement` and `DeleteDuplicates` have a machine
fast path over a rank-1 buffer of exact integers — the domain where a set
operation is most often a graph traversal or a k-mer count. It is reached from
**either** array representation, the invisible packed `List` and an explicit
`NDArray[...]`, and the result keeps whichever it was given. Until 2026-08-01
only the packed form reached it, so the same call on the same values ran 145×
slower when the argument was written as an `NDArray`.

Reals are deliberately excluded: `0.` and `-0.` compare equal and print
differently, so which of two equal elements survived would depend on the
representation. A custom `SameTest`, a non-`List` head, or any other element
type takes the general path and gives the same answer.

## Implementation notes

**Algorithm.** `builtin_union` concatenates the elements of all argument lists (which must
share a common head), sorts the combined `Expr**` array with `qsort` under the canonical
`expr_compare` order, then removes adjacent duplicates — `expr_eq` by default, or an optional
`SameTest -> f` which is evaluated per adjacent pair. The result is the sorted, deduplicated
list. (`DeleteDuplicates` in the same file does the order-preserving variant using a hash table
keyed on `expr_hash`/`expr_eq`.)

- `Flat`, `OneIdentity`, `Protected`, `ReadProtected`.
- All expressions must have the same head.
- Result has the same head as the inputs.

**Attributes:** `Flat`, `OneIdentity`, `Protected`, `ReadProtected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/), [Intersection](../../structural-manipulation/Intersection/), [Complement](../../structural-manipulation/Complement/), [DeleteDuplicates](../../data-structures/DeleteDuplicates/), [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Union[list]` gives the sorted list of distinct elements; `Union[l1, l2, ...]` gives the set union. Comparison is by canonical structural equality, and the output is always sorted with duplicates removed.
