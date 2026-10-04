---
source: src/list/subdivide.c
---
**Algorithm.** `builtin_subdivide` returns `n + 1` equally spaced points
spanning an interval, endpoints included (`n` counts the parts, not the points).
The three surface forms are normalised to one `(min, max, n)` triple — `[n]`
implies `0..1`, `[max, n]` implies `0..max` — and point `i` (0-based) is
`min + i (max - min) / n`. Descending intervals need no special case: when
`max < min` the span is negative and the points descend. `n` must be a positive
machine integer (capped at `SUBDIVIDE_MAX_N = 10^6`); anything else leaves the
call unevaluated.

**Exactness.** Two guarantees. The endpoints are *copied, never computed*
(element 0 is a copy of `min`, element `n` of `max`), so no representation change
can touch them. Interior points are each derived directly from their index `i`,
never as `previous + step`, so nothing accumulates error. For integer endpoints
inside `SUBDIVIDE_MAX_ENDPOINT` the interior point is `(min*n + i*(max-min))/n`
reduced once by `make_rational` (so whole points print as integers alongside
rationals: `Subdivide[10, 4]` is `{0, 5/2, 5, 15/2, 10}`); everything else
(bigint, rational, symbolic) builds a `Plus`/`Times`/`Power` tree and lets the
evaluator do the arithmetic, keeping exact input exact.

**Buffer fast path.** When either endpoint is a machine `Real` the whole result
is one `float64` array: `ndbuild_open_f64` fills it with `min + i*step` directly
(the endpoints written from the inputs, not computed), avoiding 10⁶ evaluator
entries — `np.linspace`-grade. A sub-threshold or packing-off case computes the
same doubles into a plain `List`. Only a machine `Real` is contagious; an MPFR
or symbolic endpoint keeps the exact/general path. `ATTR_PROTECTED`.
