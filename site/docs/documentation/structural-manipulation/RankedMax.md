# RankedMax

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RankedMax[list, n]`**

gives the n-th largest element of list.

**`RankedMax[list, -n]`**

gives the n-th smallest element of list. RankedMax\[list, n\] is RankedMin\[list, -n\]; RankedMax\[list, 1\] is Max\[list\] and RankedMax\[list, -1\] is Min\[list\]. Yields a definite result when every element is a real number. Has a packed-array fast path and is compilable.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= RankedMin[{12, 13, 11}, 2]
Out[1]= 12

In[2]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]
Out[2]= 3

In[3]:= RankedMax[{2.5, E, 12, 15, 485}, -2]
Out[3]= E

In[4]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]
Out[4]= Infinity
```

### Applications (3)

The largest element is Max

```mathematica
In[5]:= RankedMax[{12, 13, 11}, 1]
Out[5]= 13
```

The 2nd smallest, via a negative index

```mathematica
In[6]:= RankedMax[{2.5, E, 12, 15, 485}, -2]
Out[6]= E
```

Ties and infinities rank by value

```mathematica
In[7]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]
Out[7]= Infinity
```

## Implementation notes

**Algorithm.** `builtin_ranked_max` selects the `n`-th largest element of a list
— `RankedMax[list, n]`, with `RankedMax[list, -n]` the `n`-th smallest, so
`RankedMax[list, n]` is exactly `RankedMin[list, -n]` (`RankedMax[list, 1]` is
`Max[list]`, `RankedMax[list, -1]` is `Min[list]`). It validates an integer,
nonzero `n` and delegates to the shared `ranked_select` core with `is_max =
true`, which negates `n` before computing the 1-based ascending rank `r`. As for
`RankedMin`, one scan picks the exact `expr_compare` path when every element is a
real numeric atom, else the `double`-key path (`±HUGE_VAL` for `±Infinity`, the
value of a real atom, else the machine-precision numericalisation of a symbolic
real); a non-real element makes it decline. `ranked_select_idx` quickselects the
`r`-th position (`O(m)` average, Hoare partition, stable index tiebreak) and the
exact element there is returned.

**Data structures.** The same `ranked_select` machinery as `RankedMin`: a
`size_t* idx` quickselected in place, an optional `double* keys`, and a
`RankedCtx` read by `ranked_cmp`; the result is an `expr_copy` of the selected
element in its exact form.

**Complexity / limits.** `O(m)` average quickselect. An `NDArray` argument takes
the buffer order-statistic path `ndred_ranked_max` (`int64` exactly, reals via an
`O(m)` quickselect) — `RankedMax` is on `pack.c`'s `AWARE` and `INT64_OK` lists.
It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, RankedMax[v, 2]]` reports `Compiled -> True`.
`Protected`.

- `Protected`.
- `RankedMax[list, k]` is `RankedMin[list, -k]`.
- `RankedMin[list, 1]` is `Min[list]`; `RankedMin[list, -1]` is `Max[list]`.
- Yields a definite result whenever every element is a real number, including
  symbolic real constants (`Pi`, `E`, `Sqrt[2]`, `Pi + E`), which order by value;
  `Infinity`/`-Infinity` rank as `±∞`. Returns the element in its exact form.
- Exact for arbitrary-precision integers and rationals; a symbolic non-real
  element (a free symbol or a non-real complex), an empty list, or `|n|` out of
  range leaves the call unevaluated.
- Packed-array fast path (int64 exact, real via O(*n*) quickselect) and a
  `Compile[]` lowering, so `RankedMin[v, k]`/`RankedMax[v, k]` compile and
  auto-compile.

**Attributes:** `Protected`.

## References

**See also:** [RankedMin](../../structural-manipulation/RankedMin/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ranked.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ranked.c)

## Notes & additional examples

### Notes

`RankedMax[list, n]` is the `n`-th largest element, the mirror of `RankedMin`:
`RankedMax[list, n]` equals `RankedMin[list, -n]`, so `RankedMax[list, 1]` is
`Max[list]` and `RankedMax[list, -1]` is `Min[list]`; a negative index
`RankedMax[list, -n]` gives the `n`-th smallest. It returns a definite result
whenever every element is a real number (symbolic real constants order by value,
`±Infinity` as `±∞`) and returns the element in its exact form; a non-real
element, an empty list, or an out-of-range index leaves the call unevaluated.
Selection is an `O(n)` quickselect sharing `RankedMin`'s core, with a
packed-array fast path and a `Compile[]` lowering.
