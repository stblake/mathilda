# Partition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Partition[list, n]
    partitions list into non-overlapping sublists of length n; trailing
    elements that do not fill a block are discarded.
Partition[list, n, d]
    uses offset d between successive sublists; d = 1 gives a moving
    window, d = n gives non-overlapping blocks.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Algorithm.** `builtin_partition` splits a list into sublists of length `n` with offset `d`
(default `d = n`, i.e. non-overlapping blocks), via the recursive `partition_rec`. At each level
it reads the block size `n` and offset `d` for that level (a plain integer applies to level 0,
or a `List` gives a per-level spec), computes the number of full blocks `(len − n)/d + 1`, and
emits each sublist `args[i·d .. i·d + n)` wrapped in the list's head. An `UpTo[n]` size allows a
short final block. It recurses into each element so multi-level specs partition nested arrays.
Trailing partial blocks (when no `UpTo`) are dropped, matching Mathematica's no-padding default.

- `Protected`.
- Works on any expression with arguments.
- `Partition[list, n, d]` only includes full sublists of length `n` unless `UpTo` is used.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

Split a list into non-overlapping blocks of a given length (trailing
under-filled elements are dropped):

```mathematica
In[1]:= Partition[{a, b, c, d, e, f}, 2]
Out[1]= {{a, b}, {c, d}, {e, f}}

In[2]:= Partition[{1, 2, 3, 4, 5}, 2]
Out[2]= {{1, 2}, {3, 4}}
```

The offset `d` turns `Partition` into a sliding window — `d = 1` gives every
consecutive overlapping pair:

```mathematica
In[1]:= Partition[{1, 2, 3, 4, 5}, 2, 1]
Out[1]= {{1, 2}, {2, 3}, {3, 4}, {4, 5}}
```

Combined with `Map`, non-overlapping blocks let you reduce a stream in chunks —
the block sums of `1..12` taken four at a time:

```mathematica
In[1]:= Map[Total, Partition[Range[12], 4]]
Out[1]= {10, 26, 42}
```

A length-2 sliding window computes first differences. Applied to the squares
`1, 4, 9, 16, 25` it recovers the classical identity that consecutive squares
differ by successive odd numbers:

```mathematica
In[1]:= Map[(#[[2]] - #[[1]] &), Partition[{1, 4, 9, 16, 25}, 2, 1]]
Out[1]= {3, 5, 7, 9}
```

### Notes

`Partition[list, n]` cuts `list` into consecutive non-overlapping length-`n`
sublists, discarding a trailing remainder that cannot fill a full block.
`Partition[list, n, d]` advances by offset `d` between successive sublists:
`d = n` reproduces the non-overlapping blocks, while smaller `d` produces
overlapping moving windows (`d = 1` slides one element at a time). Pairing
`Partition` with `Map`/`Total` is the idiomatic way to express block reductions
and finite-difference / sliding-window computations.
