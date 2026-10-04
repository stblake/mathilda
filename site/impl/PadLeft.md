---
source: src/list/pad.c
---
**Algorithm.** `PadLeft` and `PadRight` are exact mirrors and share one recursive
engine `pr_build`, selected by a `pad_left` flag. `pad_dispatch` parses the
length spec (arg 2): an `Integer` `n` (a single level), a `List {n1, ..., nk}`
(target length `ni` at level `i`, building a full nested array), or
`Automatic`/omitted (pad the ragged input to full rectangular — `pr_scan_dims`
records the maximum width seen at each level). Padding (arg 3) defaults to the
`Integer` `0`, may be a single element, or a `List` tiled cyclically; the margin
(arg 4) leaves that many padding elements on the far side (negative truncates).
`pr_build` lays out each level into `N = |n|` slots, placing the original row
left-aligned (for right padding) or right-aligned (for left padding) and filling
the rest from the padding block via `pr_pad_at`, which indexes into a `List`
padding with floored modulo so the tiling phase is continuous. A negative length
pads on the opposite side; the head of `list` need not be `List` and is
preserved.

**Data structures.** A heap `int64_t* dimv` of per-level target lengths, an
`int64_t* coords` path accumulated down the recursion (used to index the cyclic
padding), and the per-level source row read straight off the input node's args.
`pr_is_atomic` stops the dimension scan and the builder from descending into
`Rational`/`Complex` nodes.

**Complexity / limits.** `O(total output elements)`. A rank-1 visible/packed
`NDArray` with a dtype-compatible fill (or pure truncation) takes the native
buffer path `ndstruct_pad` — `PadLeft`/`PadRight` are on `pack.c`'s `AWARE` and
`INT64_OK` lists; otherwise the array is materialised once and the `List` engine
runs (padding an exact `0` into a `float64` buffer yields a mixed list no uniform
buffer can hold, so it deliberately does not repack). There is **no** `Compile[]`
lowering (`CompileDiagnostics` reports `Compiled -> False` for this head).
`Protected`.
