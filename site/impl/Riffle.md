---
source: src/list/riffle.c
---
**Algorithm.** `builtin_riffle` interleaves separators into the gaps of a list.
`Riffle[list, x]` places `x` in every gap; `Riffle[list, {x1, ..., xk}]` consumes
the `xi` in order and cycles back to `x1`, filling gaps left to right. A list of
`n` elements has exactly `n - 1` gaps, so the output has `2n - 1` slots: the gap
following element `i` takes separator index `i mod k`, and separators past the
last gap are simply never indexed.

**Edge invariants.** `n <= 1` (no gaps) or an empty separator list copies the
input through unchanged — checked *before* the `2n - 1` sizing, since with
`n == 0` that expression underflows `size_t`. The head of the first argument is
preserved rather than forced to `List`, so `Riffle[f[a, b], x]` gives
`f[a, x, b]` and `Riffle[{}, 0]` is `{}` with no special case.

**Data structures / limits.** One pass, O(n) element copies, a single
exactly-sized allocation (the output length is known up front). A packed/`NDArray`
first argument takes the `ndstruct_riffle` buffer fast path, falling back to
`ndstruct_delist_repack`; `Riffle` is on `pack.c`'s `AWARE` list. An atom first
argument stays unevaluated. No special attributes beyond `ATTR_PROTECTED`.
