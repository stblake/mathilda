# Catenate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Catenate[{e1, e2, ...}]`**

Concatenates the ei (which must share a head) into one, flattening a single level. Associations contribute their values: Catenate\[{\<|a -\> 1|\>, \<|b -\> 2|\>}\] is {1, 2}; Catenate\[assoc\] catenates the association's values.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Catenate[{<|"a" -> 1|>, <|"b" -> 2|>}]
Out[1]= {1, 2}

In[2]:= Insert[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "d" -> 4, Key["b"]]
Out[2]= <|"a" -> 1, "d" -> 4, "b" -> 2, "c" -> 3|>

In[3]:= Insert[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "a" -> 9, 3]
Out[3]= <|"b" -> 2, "a" -> 9, "c" -> 3|>

In[4]:= Pick[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {True, False, True}]
Out[4]= <|"a" -> 1, "c" -> 3|>

In[5]:= Partition[<|"a" -> 1, "b" -> 2|>, 1]
Out[5]= Partition[<|"a" -> 1, "b" -> 2|>, 1]

In[6]:= Reverse[<|"x" -> {1, 2}, "y" -> {3, 4}|>, 2]
Out[6]= <|"x" -> {2, 1}, "y" -> {4, 3}|>
```

### Applications (4)

```mathematica
In[7]:= Catenate[{{1, 2}, {3, 4}, {5}}]
Out[7]= {1, 2, 3, 4, 5}
```

Ragged lengths are fine

```mathematica
In[8]:= Catenate[{{1}, {2, 3}, {4, 5, 6}}]
Out[8]= {1, 2, 3, 4, 5, 6}
```

Associations take part through their values

```mathematica
In[9]:= Catenate[{<|a -> 1|>, <|b -> 2|>}]
Out[9]= {1, 2}
```

The collection may itself be an association

```mathematica
In[10]:= Catenate[<|x -> {1, 2}, y -> {3}|>]
Out[10]= {1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_catenate` flattens one level of a single collection: the
elements `ei` of `Catenate[{e1, e2, ...}]` must share a head and their arguments
are concatenated under it (`{{1,2},{3,4}}` → `{1,2,3,4}`). Associations take
part through their *values* (Mathematica 15): `Catenate[{<|a->1|>, <|b->2|>}]` is
`{1, 2}`, mixed list/association parts concatenate their elements and values, and
the outer collection may itself be an association. This differs from `Join`,
which concatenates several arguments.

**Data structures.** The buffer path runs first: a *packed* argument arrives as
one rank-2 row-major array (the pack gate absorbs a list of packed vectors before
`Catenate` is called), and catenating its rows is a reshape handled by
`ndstruct_catenate`; a *visible* `NDArray` list is de-listed and re-evaluated.
Otherwise the result is an ordinary `List` built by copying each part's elements
in order.

**Complexity / limits.** `O(total elements)`. Returns `NULL` (unevaluated) for a
non-list/association argument, mixed heads among the parts, or a part that is not
a list/association when an association is involved.

- `Insert` removes an existing entry with the inserted key, so the new position
  wins; a non-rule element returns the association unchanged, and an
  out-of-range position or absent key leaves the call unevaluated
  (Mathematica 15).
- An association used as a `Pick` selector is atomic: the whole expression if
  it matches the pattern, else `Sequence[]`, as in Mathematica 15 (it no longer
  builds malformed `Rule[]` nodes).

**Attributes:** `Protected`.

## References

**See also:** [Insert](../../data-structures/Insert/), [Pick](../../data-structures/Pick/), [Partition](../../data-structures/Partition/)

- Source: [`src/list/join.c`](https://github.com/stblake/mathilda/blob/main/src/list/join.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`Catenate[{e1, e2, ...}]` flattens one level: the parts must share a head and
their elements are concatenated under it, so `Catenate[{{1,2},{3,4}}]` is
`{1,2,3,4}`. It differs from `Join`, which concatenates several *arguments*
rather than one list of parts. Associations join through their values
(`Catenate[{<|a->1|>, <|b->2|>}]` is `{1, 2}`), and the outer collection may
itself be an association. A packed matrix catenates as a reshape of its row-major
buffer, so the operation stays on the buffer.
