# Join

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Join[list1, list2, ...]`**

Concatenates lists or other expressions that share the same head.

**`Join[list1, list2, ..., n]`**

Joins the objects at level n in each of the lists. Handles ragged arrays by concatenating successive elements at level n.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>

In[5]:= Sort[<|"a" -> 2, "b" -> 3, "c" -> 1|>, Greater]
Out[5]= <|"b" -> 3, "a" -> 2, "c" -> 1|>

In[6]:= SortBy[<|"a" -> {1, 2}, "b" -> {0, 5}, "c" -> {1, 1}|>, First, Greater]
Out[6]= <|"a" -> {1, 2}, "c" -> {1, 1}, "b" -> {0, 5}|>

In[7]:= Ordering[<|"a" -> 2, "b" -> 2, "c" -> 1|>, All, Greater]
Out[7]= {2, 1, 3}

In[8]:= ReverseSort[<|"a" -> 2, "b" -> 2, "c" -> 1|>]
Out[8]= <|"a" -> 2, "b" -> 2, "c" -> 1|>
```

### Applications (3)

```mathematica
In[9]:= Join[{1, 2}, {3, 4}, {5}]
Out[9]= {1, 2, 3, 4, 5}
```

Any shared head, not just List

```mathematica
In[10]:= Join[f[a], f[b, c]]
Out[10]= f[a, b, c]
```

Associations merge, later values winning

```mathematica
In[11]:= Join[<|a -> 1|>, <|b -> 2|>]
Out[11]= <|a -> 1, b -> 2|>
```

## Implementation notes

`builtin_join` (in `src/list.c`) concatenates its arguments via the helper `join_at_level`. A trailing integer argument is interpreted as a level specification (default 1): at level 1 the arguments' top-level elements are spliced into a single result sharing the first list's head; deeper levels splice element-wise at the corresponding depth. Returns `NULL` if no lists remain or the level is below 1.

- With an ordering function, the sort is Mathematica 15's top-down merge sort
  (the merge keeps the left element unless `p[left, right]` is `False` or
  `-1`), so ties under a strict `p` such as `Greater` land exactly where
  Mathematica puts them. `Ordering[..., p]` and `KeySort[assoc, p]` share it.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [SortBy](../../data-structures/SortBy/), [Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Greater](../../comparisons/Greater/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

`Join[e1, e2, ...]` concatenates its arguments, which must share a head, into a
single expression under that head — the standard list-concatenation operator, and
it works equally on `f[...]` expressions. Associations are merged key-wise, with
later associations overriding earlier values on shared keys. A trailing integer
argument gives a level specification (default 1) so that `Join[..., n]`
concatenates at depth `n`. Use `Catenate` to flatten one list of parts rather
than several arguments.
