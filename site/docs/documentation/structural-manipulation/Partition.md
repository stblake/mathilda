# Partition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Partition[list, n]`**

partitions list into non-overlapping sublists of length n; trailing elements that do not fill a block are discarded.

**`Partition[list, n, d]`**

uses offset d between successive sublists; d = 1 gives a moving window, d = n gives non-overlapping blocks.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Partition[{a, b, c, d, e}, 2]
Out[1]= {{a, b}, {c, d}}

In[2]:= Partition[{a, b, c, d, e}, 2, 1]
Out[2]= {{a, b}, {b, c}, {c, d}, {d, e}}

In[3]:= Partition[{a, b, c, d, e}, UpTo[2]]
Out[3]= {{a, b}, {c, d}, {e}}

In[4]:= Partition[{{1, 2, 3}, {4, 5, 6}}, {2, 2}]
Out[4]= {{{{1, 2}}, {{4, 5}}}}
```

### Applications (5)

```mathematica
In[5]:= Partition[{a, b, c, d, e, f}, 2]
Out[5]= {{a, b}, {c, d}, {e, f}}

In[6]:= Partition[{1, 2, 3, 4, 5}, 2]
Out[6]= {{1, 2}, {3, 4}}

In[7]:= Partition[{1, 2, 3, 4, 5}, 2, 1]
Out[7]= {{1, 2}, {2, 3}, {3, 4}, {4, 5}}

In[8]:= Map[Total, Partition[Range[12], 4]]
Out[8]= {10, 26, 42}

In[9]:= Map[(#[[2]] - #[[1]] &), Partition[{1, 4, 9, 16, 25}, 2, 1]]
Out[9]= {3, 5, 7, 9}
```

## Implementation notes

**Algorithm.** `builtin_partition` splits a list into sublists of length `n` with offset `d`
(default `d = n`, i.e. non-overlapping blocks), via the recursive `partition_rec`. At each level
it reads the block size `n` and offset `d` for that level (a plain integer applies to level 0,
or a `List` gives a per-level spec), computes the number of full blocks `(len − n)/d + 1`, and
emits each sublist `args[i·d .. i·d + n)` wrapped in the list's head. An `UpTo[n]` size allows a
short final block. It recurses into each element so multi-level specs partition nested arrays.
Trailing partial blocks (when no `UpTo`) are dropped, following the no-padding default.

- `Protected`.
- Works on any expression with arguments.
- `Partition[list, n, d]` only includes full sublists of length `n` unless `UpTo` is used.

**Attributes:** `Protected`.

## References

**See also:** [UpTo](../../structural-manipulation/UpTo/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Notes

`Partition[list, n]` cuts `list` into consecutive non-overlapping length-`n`
sublists, discarding a trailing remainder that cannot fill a full block.
`Partition[list, n, d]` advances by offset `d` between successive sublists:
`d = n` reproduces the non-overlapping blocks, while smaller `d` produces
overlapping moving windows (`d = 1` slides one element at a time). Pairing
`Partition` with `Map`/`Total` is the idiomatic way to express block reductions
and finite-difference / sliding-window computations.
