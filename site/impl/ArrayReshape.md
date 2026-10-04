---
source: src/list/array_reshape.c
---
**Algorithm.** `builtin_array_reshape` arranges the flattened elements of a list
into a rectangular `dims` array. `ar_parse_dims` reads the dims spec (a
non-negative integer, or a `List` of them) and computes the element total with an
overflow guard. `ar_collect` then fully flattens the input, descending `List`
heads only, into a growing buffer of borrowed leaf pointers. The output is forced
to exactly `total` owned leaves — copied/truncated when the input has enough, or
extended by `pad_scheme_extend` (default fill `0`, or a given padding scheme)
when it has too few — and `ar_build` folds that flat sequence into the nested
rectangular shape in row-major order. So up to the shared length,
`Flatten[ArrayReshape[list, dims]] == Flatten[list]`.

**Data structures / limits.** Rank is capped at `AR_MAX_RANK` (64); the leaf
buffer starts at 16 and doubles. A packed/`NDArray` first argument takes the
`ndstruct_arrayreshape` buffer fast path (a dims-header `memcpy`, no
per-element work); `ArrayReshape` is on `pack.c`'s `AWARE` list.

**Complexity / limits.** O(total) element copies plus the flatten pass.
`ATTR_PROTECTED`. A bad dims spec (empty list, negative or non-integer entry,
rank overflow, or an int64-overflowing total) or a non-list, non-NDArray first
argument leaves the call unevaluated.
