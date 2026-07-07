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
In[1]:= Tally[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[1]= {{1, 2}, {2, 1}}

In[2]:= Commonest[<|"a" -> 1, "b" -> 1, "c" -> 2|>]
Out[2]= {1}
```

## Implementation notes

**Algorithm.** `builtin_commonest` (in `src/list.c`) tallies element multiplicities with a `HashTable` (one pass, O(N)), recording each distinct element's count and first-appearance index. The tallies are packed into `CommonestItem` records and sorted by descending count (ties broken by first occurrence) to find the maximum multiplicity. Without a count argument it returns every element sharing the maximum count; `Commonest[list, n]` / `Commonest[list, UpTo[n]]` returns up to the `n` most frequent.

**Data structures.** Parallel `Expr**` / `int64_t*` arrays for unique elements and their multiplicities, a `HashTable` mapping element → index, and a `CommonestItem` array (`element`, `count`, `first_index`) used for the stable sort.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

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
