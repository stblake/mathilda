---
source: src/list/array_pad.c
---
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
