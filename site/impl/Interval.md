---
source: src/interval.c
references:
  - "R. E. Moore, R. B. Kearfott and M. J. Cloud, *Introduction to Interval Analysis* (SIAM, 2009)."
---
**Algorithm.** `Interval[{min, max}, ...]` represents a set of real values — the union of
the given inclusive ranges. `builtin_interval` (`src/interval.c`) canonicalises: it orders
each endpoint pair, sorts the ranges by their lower endpoint, and merges overlapping or
touching ranges, so `Interval[{1,3},{2,5}]` collapses to `Interval[{1,5}]`. It is a plain
`EXPR_FUNCTION` with head `SYM_Interval` (the `Complex`/`Rational` pattern), recognised by
`is_interval`, and marked `Protected` only — deliberately **not** `NumericFunction` or
`Listable`, because an interval is a *set*, not a number threaded over.

Arithmetic and elementary functions thread through intervals to produce **rigorous
enclosures**. The endpoint arithmetic goes through the ordinary evaluator (`Plus`, `Times`,
`Power`, `Sin`, …), so exact endpoints stay exact or symbolic; only an inexact (Real/MPFR)
result is nudged one ULP **outward** by `iv_widen` (`nextafter`, or `mpfr_nextbelow`/
`_nextabove`), which keeps the enclosure sound without any MPFI dependency. The hooks live
at the real dispatch sites: `plus.c`/`times.c`/`power.c` for arithmetic, and
`interval_apply_function` / the general derivative-sign certifier (`interval_thread_call`,
wired into `eval.c` after a builtin returns `NULL`) for monotone and certified-monotone
unary functions. Non-monotone heads such as `Sin`/`Cos` pin an extremum to exact `±1` when a
crest or trough lies inside the range.

**Data structures.** Intervals are ordinary `Expr` trees; endpoint comparison
(`interval_endpoint_cmp`) is done **exactly** — both endpoints are numericalised to 200-bit
MPFR and compared with `mpfr_cmp`, never through the tolerant comparison heads, since a
~1-ULP tolerance there once kept the inner of two candidate endpoints and broke containment.

**Complexity / limits.** Canonicalisation is `O(k log k)` in the number `k` of ranges. The
general certifier is depth-capped (`IV_CERTIFY_MAX_DEPTH`) so oscillatory derivative chains
(Bessel) terminate symbolically rather than looping. `Sqrt` is lowered to `Power[·, 1/2]` at
parse time, so it threads through the power hook, not a `Sqrt` guard. The inclusion property
is fuzzed by `tools/interval_fuzz.py` (`make check-interval`).
