---
source: src/vectors.c
---
**Algorithm.** `builtin_unit_vector` builds the length-`n` vector with a `1` in
position `k` and `0` elsewhere. `UnitVector[k]` is the 2-D unit vector
(`UnitVector[2, k]`); `UnitVector[n, k]` the general form. It counts the leading
required (non-rule) arguments, reports a trailing non-option via
`UnitVector::nonopt` and a zero-argument call via `UnitVector::argt`, and
requires `1 <= k <= n` with both positive machine integers — out-of-range or
non-integer arguments leave the call unevaluated (symbolic arguments flow
through).

**Component precision.** The `WorkingPrecision` option selects the component
representation, mirroring `HilbertMatrix`: `Infinity` (default) gives exact
integer components, `MachinePrecision` gives machine `Real`s, and a digit count
`d` above machine precision gives `d`-digit MPFR reals (degrading to machine
precision when `USE_MPFR=0`, with one `UnitVector::wprec` warning). Last valid
setting wins; an unparseable value is ignored.

**Buffer fast path.** For the exact-integer and machine-real modes every
component shares one dtype, so `ndbuild_open` opens a packed `NDArray`
(`NDT_INT64` / `NDT_FLOAT64`), writes zeros, and sets the single `1` — making
`UnitVector[10^6, 1]` a buffer fill rather than 10⁶ boxed nodes for one nonzero
(130 ms → the NumPy-grade path). MPFR components have no buffer form and keep the
`List`. `ndbuild_open` declining (packing off, small `n`) falls through to the
boxed `List`. `ATTR_PROTECTED`.
