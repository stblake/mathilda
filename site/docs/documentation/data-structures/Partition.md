# Partition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Partition[list, n]`**

partitions list into non-overlapping sublists of length n; trailing elements that do not fill a block are discarded.

**`Partition[list, n, d]`**

uses offset d between successive sublists; d = 1 gives a moving window, d = n gives non-overlapping blocks.

## Examples (11)

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

### Applications (5)

```mathematica
In[7]:= Partition[{a, b, c, d, e, f}, 2]
Out[7]= {{a, b}, {c, d}, {e, f}}

In[8]:= Partition[{1, 2, 3, 4, 5}, 2]
Out[8]= {{1, 2}, {3, 4}}

In[9]:= Partition[{1, 2, 3, 4, 5}, 2, 1]
Out[9]= {{1, 2}, {2, 3}, {3, 4}, {4, 5}}

In[10]:= Map[Total, Partition[Range[12], 4]]
Out[10]= {10, 26, 42}

In[11]:= Map[(#[[2]] - #[[1]] &), Partition[{1, 4, 9, 16, 25}, 2, 1]]
Out[11]= {3, 5, 7, 9}
```

## Implementation notes

**Algorithm.** `builtin_partition` splits a list into sublists of length `n` with offset `d`
(default `d = n`, i.e. non-overlapping blocks), via the recursive `partition_rec`. At each level
it reads the block size `n` and offset `d` for that level (a plain integer applies to level 0,
or a `List` gives a per-level spec), computes the number of full blocks `(len − n)/d + 1`, and
emits each sublist `args[i·d .. i·d + n)` wrapped in the list's head. An `UpTo[n]` size allows a
short final block. It recurses into each element so multi-level specs partition nested arrays.
Trailing partial blocks (when no `UpTo`) are dropped, following the no-padding default.

- `Insert` removes an existing entry with the inserted key, so the new position
  wins; a non-rule element returns the association unchanged, and an
  out-of-range position or absent key leaves the call unevaluated
  (Mathematica 15).
- An association used as a `Pick` selector is atomic: the whole expression if
  it matches the pattern, else `Sequence[]`, as in Mathematica 15 (it no longer
  builds malformed `Rule[]` nodes).

**Attributes:** `Protected`.

## References

**See also:** [Catenate](../../data-structures/Catenate/), [Insert](../../data-structures/Insert/), [Pick](../../data-structures/Pick/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`Partition[list, n]` cuts `list` into consecutive non-overlapping length-`n`
sublists, discarding a trailing remainder that cannot fill a full block.
`Partition[list, n, d]` advances by offset `d` between successive sublists:
`d = n` reproduces the non-overlapping blocks, while smaller `d` produces
overlapping moving windows (`d = 1` slides one element at a time). Pairing
`Partition` with `Map`/`Total` is the idiomatic way to express block reductions
and finite-difference / sliding-window computations.
