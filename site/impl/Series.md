---
references:
  - "Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003)."
source: src/calculus/series.c
---
**Algorithm.** Series computes a truncated power series via recursive *series
algebra*. `builtin_series` (`ATTR_HOLDALL`) parses each spec — full form
`{x, x0, n}` or leading-term form `x -> x0` (via `parse_series_spec`) — threads a
`List` first argument, and handles the multivariate form by expanding left-to-
right (each coefficient of the outer series is recursively expanded in the next
variable). Inexact inputs are rationalised then numericalised back
(`internal_rationalize_then_numericalize`).

`do_series_single` does the work. It evaluates `f`, returns it verbatim if free of
`x`, optionally Apart-decomposes a rational function in `x` (each partial-fraction
term hits the cheap monomial binomial path instead of a Newton inversion), pulls
out a symbolic `x^alpha` prefactor when expanding at 0 or Infinity, and handles
expansion at Infinity by substituting `x -> 1/u`. It expands to a padded internal
order (`order + pad`, pad = 12 for numeric `x0`, 2 for symbolic `x0` to keep
symbolic convolution from exploding) and then truncates back to the user O-term;
the leading-term form widens the O-term to the first non-zero coefficient.

The recursion is `series_expand(e, ctx)`: a subexpression free of `x` becomes a
constant series; `x` itself becomes the identity series; `Plus`/`Times`/`Power`
combine child series with `so_add`/`so_mul`/`so_inv`/`so_pow_int`; the elementary
heads `Exp, Log, Sin, Cos, Tan, Sinh, Cosh, Tanh` and the inverse trig/hyperbolic
family compose their known Taylor kernels (`kernel_coefs` →
`so_compose_scalar_kernel`) with the inner argument's series, with reciprocal
heads (`Sec`/`Csc`/`Cot`/...) rewritten via `rewrite_reciprocal_head` and a large
set of branch-point / at-infinity identities for the inverse functions. Any
unrecognised head falls back to naive Taylor `series_taylor_via_D`: coefficients
`a_k = (D^k f at x0)/k!`, capped at `MAX_NAIVE_ORDER` = 20 and bailing out
(`has_infinity`) at singularities.

**Data structures.** The internal `SeriesObj` struct holds the expansion variable
`x`, expansion point `x0`, an owned array of coefficient `Expr*`, the leading
exponent numerator `nmin`, the (exclusive) O-term exponent numerator `order`, and
a common exponent denominator `den` (>= 1) — so coefficient `i` multiplies
`(x - x0)^((nmin+i)/den)`, supporting Laurent and Puiseux (fractional-exponent)
series. `so_rescale`/`so_align_den` reconcile denominators before arithmetic.
`SeriesObj` is converted to the user-facing `SeriesData[x, x0, {coefs}, nmin,
nmax, den]` head by `so_to_expr`.

**Complexity / limits.** Kernel-path heads are exact and fast; the naive D-path is
capped at order 20 and fails at true branch points where derivatives blow up
(Puiseux at such points is out of scope except for the explicit inverse-function
branch-point handlers). Symbolic expansion points cap the internal pad tightly to
avoid `O(N^2)` symbolic coefficient blow-up.
