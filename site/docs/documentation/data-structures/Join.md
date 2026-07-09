# Join

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Join[list1, list2, ...]
    Concatenates lists or other expressions that share the same head.
Join[list1, list2, ..., n]
    Joins the objects at level n in each of the lists.
    Handles ragged arrays by concatenating successive elements at level n.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>
```

## Implementation notes

`builtin_join` (in `src/list.c`) concatenates its arguments via the helper `join_at_level`. A trailing integer argument is interpreted as a level specification (default 1): at level 1 the arguments' top-level elements are spliced into a single result sharing the first list's head; deeper levels splice element-wise at the corresponding depth. Returns `NULL` if no lists remain or the level is below 1.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
