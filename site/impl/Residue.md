---
source: src/calculus/residue.c
---
**Algorithm.** `builtin_residue` handles `Residue[f, {z, z0}]` by reading the
coefficient of `(z - z0)^-1` out of the Laurent expansion of `f`.
`residue_compute` first classifies the integrand with `residue_is_rational_in`
(is `Together[f]` a ratio of polynomials in `z`?). For a rational `f` it tries a
simple-pole fast path, `residue_simple_pole`: with `P/Q = Together[f]` and `Q`
having a simple zero at `z0` (`Q(z0) == 0` via `PossibleZeroQ`, `Q'(z0) != 0`),
the residue is `P(z0)/Q'(z0)`, which bypasses the series inverter entirely.
Failing that, `residue_shift_form` substitutes `z -> z0 + w` and `Expand`s the
shifted numerator and denominator *separately* (keeping them coprime), then
`residue_extract` runs `Series[..., {w, 0, 0}]` and reads the `w^-1` term.
Transcendental integrands (`Cot`, `Zeta` near its pole, unknown `f[z]`) skip the
shift and expand directly about `z0`, so the series engine can use its built-in
knowledge of the function's Laurent series there.

**Data structures.** `Expr` trees throughout, driven by `eval_and_free` wrappers
(`residue_eval1`/`residue_eval2`) that build and evaluate `Together`,
`Numerator`, `Denominator`, `Expand`, `D`, `ReplaceAll`, and `Series`. The
expansion returns a `SeriesData[z, z0, {coefs}, nmin, nmax, den]`;
`residue_extract` reads the coefficient at index `-1 - nmin`. A fresh local
`w = Residue\`$w` carries the shifted expansion.

**Complexity / limits.** The simple-pole path is one differentiation plus two
substitutions; the series path costs a Laurent expansion whose order
`residue_extract` raises adaptively (capped at 256) until the `-1` term is
resolved. A residue is defined only for an ordinary Laurent expansion
(`den == 1`): a fractional-power (Puiseux) expansion signals a branch point, so
`Residue[1/Sqrt[z], {z, 0}]` is left unevaluated. Only the two-argument form is
handled; fewer than two arguments emit `Residue::argm`. The separate-shift
design is what keeps an algebraic pole location (a `z0` that is a sum of
radicals) tractable, since re-cancelling the shifted ratio can spawn a second
spelling of the same algebraic number and blow the series inversion up
combinatorially.
