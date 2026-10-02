---
references:
  - "B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26."
  - "W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11."
  - "DLMF §19.25.5 — `F(φ|m) = s R_F(1−s², 1−m s², 1)`, `s = sin φ`."
source: src/special_functions/elliptic.c
---
**Algorithm.** `builtin_ellipticf` runs the numeric path first (inexact in, inexact out),
then the exact reductions `F(0, m) = 0`, `F(φ, 0) = φ`, `F(φ, ∞) = 0`,
`F(π/2, m) = K(m)`, and finally the oddness fold in the amplitude. The machine kernel is
`elliptic_inc_real(φ, m, want_E = false)`:

1. reduce the amplitude to the principal strip — `k = ⌊φ/π + 1/2⌋`, `r = φ − kπ`, so
   `|r| ≤ π/2`;
2. evaluate `F(r|m) = s · R_F(cos²r, 1 − m s², 1)` with `s = sin r` (DLMF 19.25.5, after
   using `R_F`'s homogeneity to clear the `csc²` scaling);
3. restore the shift with the quasi-period `F(φ + kπ | m) = F(φ|m) + 2k K(m)`, which calls
   `elliptic_machine_k` — so a shifted amplitude declines whenever the complete form does.

Outside the real principal domain (`1 − m sin²φ < 0`, or a complex amplitude, which
`EllipticF[ArcSin[z], m]` produces as soon as `|z| > 1`) the kernel declines and
`flint_num_elliptic_f` answers through Arb's `acb_elliptic_f`, which carries the
quasi-period and the branch placement itself.

**One spelling carries eight digits.** The first argument of `R_F` is `Cos[r]^2` and must
be computed as `c*c`, never as `1 - s*s`: near `r = π/2` that subtraction's absolute error
is the size of the true value, so `R_F`'s first argument arrives with ~100% relative
error. Measured against the 30-digit path at the same `double`:
`EllipticF[π/2 − 1e−8, 0.99]` was wrong by **2.7e-08** relative (1.2e8 ulp) against
6.7e-16 at a generic amplitude; it is now 1.7e-15 worst case over a π/2 approach ladder.
A single `fmax`-style clamp (`if (a < 0) a = 0`) covers the one case squaring cannot — an
`r` at which `cos r` underflows to zero, where a rounding-negative `a` would make
`carlson_rf` decline at exactly the amplitude where `F` is simply `K`.

**Data structures.** `double` scalars and the shared `carlson_rf` duplication loop; `acb_t`
on the Arb path. `EllipticF` is registered as a binary (`REG_B`) ND kernel, so a packed or
visible `NDArray` pair runs element-wise through the same `double` code.

**Complexity / limits.** `O(1)` per element, ~3 ulp; `Compile[]` lowers it at scalar and
rank-1 shapes. `D[EllipticF[φ,m], m]` is a closed form (checked against a central
difference at 30 digits before it was allowed in), unlike the `EllipticPi` parameter
derivatives, which stay inert.
