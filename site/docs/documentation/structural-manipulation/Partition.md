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
