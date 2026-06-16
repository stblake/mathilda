# Commonest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Commonest[list] gives a list of the elements that are the most common in list.
Commonest[list, n] gives a list of the n most common elements in list.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Commonest[{b, a, c, 2, a, b, 1, 2}]
Out[1]= {b, a, 2}

In[2]:= Commonest[{b, a, c, 2, a, b, 1, 2}, 4]
Out[2]= {b, a, c, 2}

In[3]:= Commonest[{b, a, c, 2, a, b, 1, 2}, UpTo[6]]
Out[3]= {b, a, c, 2, 1}

In[4]:= Commonest[{1, 2, 2, 3, 3, 3, 4}]
Out[4]= {3}

In[5]:= Commonest[{a, E, Sin[y], E, a, 7}]
Out[5]= {a, E}
```

## Implementation notes

**Algorithm.** `builtin_commonest` (in `src/list.c`) tallies element multiplicities with a `HashTable` (one pass, O(N)), recording each distinct element's count and first-appearance index. The tallies are packed into `CommonestItem` records and sorted by descending count (ties broken by first occurrence) to find the maximum multiplicity. Without a count argument it returns every element sharing the maximum count; `Commonest[list, n]` / `Commonest[list, UpTo[n]]` returns up to the `n` most frequent.

**Data structures.** Parallel `Expr**` / `int64_t*` arrays for unique elements and their multiplicities, a `HashTable` mapping element → index, and a `CommonestItem` array (`element`, `count`, `first_index`) used for the stable sort.

- `Protected`.
- When several elements occur with equal frequency, `Commonest` picks first the ones that occur first in `list`.
- `Commonest[list, n]` returns the `n` commonest elements in the order they appear in `list`.
- `Commonest[list, UpTo[n]]` returns the `n` commonest elements, or as many as are available.
- A message `Commonest::dstlms` is generated if there are fewer distinct elements than requested by an integer `n`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Commonest[{1, 2, 2, 3, 3, 3, 4}]
Out[1]= {3}
```

```mathematica
In[1]:= Commonest[{a, b, a, c, b, a, d, b}, 2]
Out[1]= {a, b}
```

```mathematica
In[1]:= Commonest[Table[Mod[k^2, 7], {k, 0, 20}]]
Out[1]= {1, 4, 2}
```

### Notes

`Commonest[list]` returns every element tied for the highest frequency; `Commonest[list, n]` returns the `n` most common. The quadratic-residue example above shows the three nonzero residues mod 7 each occurring equally often.
