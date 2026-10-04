# ArrayPad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayPad[array, m]`**

pads array with m elements of 0 on every side of every level.

**`ArrayPad[array, {m, n}]`**

pads with m elements at the start and n at the end of each dimension.

**`ArrayPad[array, {{m1, n1}, {m2, n2}, ...}]`**

pads with mi, ni elements at level i; a negative amount removes elements.

**`ArrayPad[array, amounts, padding]`**

uses the given padding: a constant c, a cyclic list {c1, c2, ...}, or one of "Fixed", "Periodic", "Reflected", "Reversed", "ReversedNegation", "ReflectedDifferences", "ReversedDifferences", "Extrapolated" (which takes the option InterpolationOrder).

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= ArrayPad[{1, 2, 3}, 1]
Out[1]= {0, 1, 2, 3, 0}

In[2]:= ArrayPad[{1, 2, 3}, 2, "Fixed"]
Out[2]= {1, 1, 1, 2, 3, 3, 3}

In[3]:= ArrayPad[{a, b, c}, 3, "Extrapolated"]
Out[3]= {4 a - 3 b, 3 a - 2 b, 2 a - b, a, b, c, -b + 2 c, -2 b + 3 c, -3 b + 4 c}

In[4]:= ArrayPad[Range[10], -2]
Out[4]= {3, 4, 5, 6, 7, 8}
```

### Applications (4)

One layer of zeros around a matrix

```mathematica
In[5]:= ArrayPad[{{1, 2}, {3, 4}}, 1]
Out[5]= {{0, 0, 0, 0}, {0, 1, 2, 0}, {0, 3, 4, 0}, {0, 0, 0, 0}}
```

Two before, none after

```mathematica
In[6]:= ArrayPad[{1, 2, 3, 4}, {2, 0}]
Out[6]= {0, 0, 1, 2, 3, 4}
```

Wrap the array around cyclically

```mathematica
In[7]:= ArrayPad[{1, 2, 3}, 2, "Periodic"]
Out[7]= {2, 3, 1, 2, 3, 1, 2}
```

Mirror the edge values outward

```mathematica
In[8]:= ArrayPad[{1, 2, 3}, 2, "Reflected"]
Out[8]= {3, 2, 1, 2, 3, 2, 1}
```

## Implementation notes

**Algorithm.** `builtin_array_pad` adds (or, for a negative amount, removes)
padding around a nested-`List` array. The amount spec is parsed by
`ap_parse_amounts` into per-level `lo[]`/`hi[]` arrays: an integer `m` pads `m`
on every side of every dimension, `{m, n}` puts `m` before and `n` after on each
dimension, and `{{m1,n1}, ...}` gives per-level amounts. A trailing
`InterpolationOrder` option is stripped first (`options_extract`).

**Two builders.** For a constant or cyclic-list padding (default `0`),
`ap_build_const` constructs each level at rectangular target width
`orig_dim + lo + hi`, copying original elements where the shifted index lands
inside the source and otherwise cyclically indexing the padding block
(`ap_pad_at`). For a named value-dependent scheme — `"Fixed"`, `"Periodic"`,
`"Reflected"`, `"Reversed"`, `"ReversedNegation"`, `"ReflectedDifferences"`,
`"ReversedDifferences"`, `"Extrapolated"` (classified in `pad_schemes.c`) —
`ap_pad_valdep` extends each axis's fiber by `pad_scheme_extend` outer-to-inner,
with an `ArrayPad::mindimsize` guard for difference schemes on an axis shorter
than 2.

**Data structures / limits.** Rank is capped at `AP_MAX_RANK` (64); amounts live
in two fixed `int64_t[64]` arrays. A packed/`NDArray` argument takes the rank-1
`ndstruct_arraypad` buffer fast path, otherwise it is materialised once with
`ndarray_to_nested_list`; `ArrayPad` is on `pack.c`'s `AWARE` list (whole
elements move by `memcpy`). `ATTR_PROTECTED`. A bad amounts spec or non-list
array leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/list/array_pad.c`](https://github.com/stblake/mathilda/blob/main/src/list/array_pad.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_array_pad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_pad.c)

## Notes & additional examples

### Notes

`ArrayPad[array, m]` pads `m` elements on every side of every level;
`ArrayPad[array, {m, n}]` puts `m` before and `n` after; and
`ArrayPad[array, {{m1,n1}, ...}]` gives per-level amounts. A **negative** amount
removes elements from that side instead.

The optional third argument is the padding: a constant (default `0`), a cyclic
list of constants, or a named scheme — `"Fixed"`, `"Periodic"`, `"Reflected"`,
`"Reversed"`, `"ReversedNegation"`, `"ReflectedDifferences"`,
`"ReversedDifferences"`, or `"Extrapolated"` (which honours
`InterpolationOrder`). The difference schemes need an axis of length at least 2
and otherwise emit `ArrayPad::mindimsize`.
