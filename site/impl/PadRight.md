---
source: src/list/pad.c
---
**Algorithm.** `PadRight` is the exact mirror of `PadLeft`, sharing the recursive
`pr_build` engine with `pad_left == 0`. `pad_dispatch` parses the same specs: an
`Integer` length `n`, a `List {n1, ..., nk}` of per-level lengths (building a full
nested array, with a nested padding block tiled), or `Automatic`/omitted (pad a
ragged array to full rectangular, dimensions found by `pr_scan_dims`). The
default padding is the `Integer` `0`; a single element repeats, a `List` repeats
cyclically (via `pr_pad_at`'s floored-modulo indexing); the margin argument
leaves padding on the *left* (`PadRight`) and a negative margin truncates leading
elements. For a non-negative length `PadRight` places the original row
left-aligned and pads on the right; a negative length pads on the left. The head
of `list` need not be `List`.

**Data structures.** Identical to `PadLeft`: a heap `int64_t* dimv` of per-level
lengths, an `int64_t* coords` path threaded down the recursion to drive the
cyclic padding, and per-level source rows read directly from the input args;
`pr_is_atomic` keeps the scan and builder out of `Rational`/`Complex` nodes.

**Complexity / limits.** `O(total output elements)`. A rank-1 visible/packed
`NDArray` with a compatible fill (or pure truncation) takes the native buffer
path `ndstruct_pad` — `PadRight` is on `pack.c`'s `AWARE` and `INT64_OK` lists;
otherwise the array is materialised once and the `List` engine runs (it does not
repack, because padding an exact `0` into a `float64` buffer gives a mixed list).
There is **no** `Compile[]` lowering (`CompileDiagnostics` reports
`Compiled -> False`). `Protected`.
