---
references:
  - "B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26."
  - "W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11 — the `rf`/`rd` duplication loops."
  - "M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions* (Dover, 1964), ch. 17."
  - "DLMF §19.25 — Legendre integrals in terms of Carlson's symmetric forms."
source: src/special_functions/elliptic.c
---
**Algorithm.** `builtin_elliptick` tries, in this order: the pole (`m == 1` →
`ComplexInfinity`), an exact zero (`π/2`), infinity (`K(m) ~ (π − i log 16m)/(2√m)`,
so the magnitude vanishes and the answer is `0`), then the two lemniscatic singular
values `m = 1/2` and `m = −1` as closed forms in `Γ` (`ell_k_closed_half`,
`ell_k_closed_minus_one`). Only then does an inexact argument reach a number:
`elliptic_machine_k` evaluates `K(m) = R_F(0, 1−m, 1)` in `double` and the result is
returned as an `EXPR_REAL`, so `Precision[EllipticK[0.5]]` is `MachinePrecision` and the
value packs. If the kernel declines — `m ≥ 1`, where the value is a pole or genuinely
complex — `flint_num_elliptic_k` answers through Arb's `acb_elliptic_k`.

**Order is load-bearing, twice.** The pole is tested *before* the numeric path because
Arb returns a non-finite ball there and the bridge maps that to `NULL`, which would leave
`EllipticK[1.]` unevaluated instead of `ComplexInfinity`; `ell_is_one` therefore also
knows `EXPR_MPFR`, so `SetPrecision[1, 30]` and `1.` cannot disagree about being at a
pole. Conversely the numeric path runs *before* the remaining exact reductions, so an
inexact argument gets an inexact answer (`EllipticK[0.]` is `1.5707963267948966`, not the
exact `π/2`) — the reductions stay underneath as the fallback for spellings Arb cannot
read, such as an exact `Pi/2` upper limit.

**Data structures.** Plain `double` scalars: `carlson_rf` runs the duplication
`x,y,z → (x+λ)/4` with `λ = √x√y + √y√z + √z√x`, stopping when every relative deviation
from the mean is below `EC_ERRTOL_RF = 0.0025`, then adds the fifth-order tail in the
`E₂`/`E₃` invariants. The tolerance is not cosmetic: the tail is fifth order, so the error
goes as `(q/A)⁶` and the previous `0.01` bought 1e-12 where the comment claimed 1e-16.
The arbitrary-precision path is `acb_t` throughout, with the bridge's accuracy ladder
(`nb_eval1`) doubling the working precision while `acb_rel_accuracy_bits` falls short of
the digits requested.

**Complexity / limits.** Seven duplication steps at most, so the machine path is `O(1)`
with ~3 ulp error and no allocation; the `Listable`/NDArray path runs it element-wise
across cores. The machine kernel is real-domain only (`m < 1`); everything else — complex
`m`, `m > 1`, and every request above machine precision — is Arb's, which is rigorous but
single-threaded. `Series` at `m = 0` comes from a dedicated kernel rather than
Taylor-via-`D`, because `d/dm K` at `m = 0` is `(π/2 − π/2)/0`.
