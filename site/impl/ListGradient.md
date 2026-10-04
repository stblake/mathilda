---
references:
  - "B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. **51** (1988) 699-706."
source: src/list/list_gradient.c
---
**Algorithm.** `builtin_list_gradient` is a from-first-principles port of
`numpy.gradient`: the numerical gradient of a sampled rank-*k* array by finite
differences, one same-shape array per differentiated axis (a rank-1 `f` gives
one vector, a rank-*k* `f` gives `{g_1, ..., g_k}`). The default reproduces
numpy exactly — a second-order central difference in the interior and a
first-order one-sided difference at the two extreme endpoints (numpy
`edge_order = 1`). One kernel underlies *every* Method/order/window/grid
combination: Fornberg's finite-difference-weight recurrence, which returns the
first-derivative weights `w[]` of `f'(x0) ~= sum_k w[k] f(z[k])` on arbitrarily
spaced nodes — so a coordinate vector (a non-uniform grid) needs no special
case. `lg_stencil` picks the node window per output index: symmetric in the
interior, shifting inward near an edge, and for the Centered method dropping the
two extreme endpoints to a one-sided `(m-1)`-point stencil; Forward/Backward
stay order `p` everywhere. Options are `Method -> "Centered"|"Forward"|
"Backward"`, `DifferenceOrder -> p` (default 2), `WindowLength -> m`
(default `p+1`), and `Axis -> All | a | {a1, ...}`.

**Data structures.** Two instantiations of the one recurrence.
`fd_weights_double` drives the machine-buffer fast path — a `float64`/`float32`
NDArray, packed or visible; its interior run is a fused tap-outer `axpy` over the
contiguous block (special-cased for the common 2- and 3-point windows) that the
compiler turns into SIMD FMA, with uniform-spacing weight vectors memoised by
`(count, offset)`. A complex machine buffer is carried through the same real
weights over its interleaved `(re, im)` doubles. `fd_weights_expr` builds `Expr`
weights the evaluator reduces — exact Rationals for an integer/rational grid,
symbolic weights for a symbolic grid — so integer and symbolic inputs answer
exactly. The buffer path declines (returns `NULL`) for anything it cannot
represent in a float buffer (int64/bool dtype, symbolic spacing), and
`ndstruct_delist_repack` re-enters the exact/symbolic List path with a
materialised List. On a packed multi-axis call the stacked rank-(*k*+1) output is
built in place to skip the transparency gate's copy; a visible NDArray keeps the
List-of-arrays shape.

**Complexity / limits.** `O(N * m)` for `N` elements and window `m`.
`ListGradient` is on `pack.c`'s `AWARE` list but deliberately **not**
`INT64_OK`: the central difference of exact integers is a Rational
(`(f[i+1]-f[i-1])/2`), which no int64 buffer holds, so an int64 argument
materialises and the exact List path answers. The rank-1 real-array,
default-centered form lowers in `Compile[]` (`compile_ndtables.c`) and
auto-compiles. A differentiated axis of length `< 2` is undefined and the call
declines.
