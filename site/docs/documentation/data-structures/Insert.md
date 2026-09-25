# Insert

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Insert[expr, elem, n] inserts elem at position n in expr.`**

## Examples (6)

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

## Implementation notes

`builtin_insert` (in `src/part.c`) handles the 3-arg form `Insert[expr, elem, pos]` by delegating to the helper `expr_insert`, which inserts `elem` before position `pos` (negative indices count from the end, and a position may be a list path for nested insertion), deep-copying the surrounding structure and preserving the original head.

- `Insert` removes an existing entry with the inserted key, so the new position
  wins; a non-rule element returns the association unchanged, and an
  out-of-range position or absent key leaves the call unevaluated
  (Mathematica 15).
- An association used as a `Pick` selector is atomic: the whole expression if
  it matches the pattern, else `Sequence[]`, as in Mathematica 15 (it no longer
  builds malformed `Rule[]` nodes).

**Attributes:** none registered.

## References

**See also:** [Catenate](../../data-structures/Catenate/), [Pick](../../data-structures/Pick/), [Partition](../../data-structures/Partition/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_part.c`](https://github.com/stblake/mathilda/blob/main/tests/test_part.c)
