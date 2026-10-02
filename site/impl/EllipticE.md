---
references:
  - "B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26."
  - "W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11."
  - "DLMF §19.25.7 — `E(φ|m) = s R_F(...) − (m/3) s³ R_D(...)`."
source: src/special_functions/elliptic.c
---
**Algorithm.** The head is arity-overloaded, and `builtin_elliptice` dispatches on argument
count to `elliptice_complete` (`E[m]`) or `elliptice_incomplete` (`E[φ, m]`); a wrong count
emits `EllipticE::argt`. Both run the numeric path first, then their reductions —
`E[0] = π/2`, `E[1] = 1`, `E[∞] = ComplexInfinity` for the complete form; `E[0, m] = 0`,
`E[φ, 0] = φ`, `E[π/2, m] = E[m]` and the oddness fold for the incomplete one.

The kernels are Carlson's:

- complete — `E(m) = R_F(0, 1−m, 1) − (m/3) R_D(0, 1−m, 1)`, with `m == 1` answered
  directly as `1`;
- incomplete — `elliptic_inc_real(φ, m, want_E = true)`, the same principal-strip
  reduction `EllipticF` uses, plus the `R_D` term, and the quasi-period closed with
  `E(m)` in place of `K(m)`.

**`R_F` and `R_D` are one loop.** Every caller that wants `R_D` wants `R_F` at the same
arguments, and the two recurrences walk an identical `xₘ, yₘ, zₘ` sequence, so
`carlson_rf_rd` computes both from one duplication — three square roots per step instead
of six. That is what pays for the tighter `EC_ERRTOL_RD = 0.0015`: `EllipticE` over 10⁶
elements is back under its pre-fix cost with 18× the accuracy (106 ulp → 6).

**`E[φ, 1]` is `Sin[φ]` only on the principal strip.** `E(φ|1) = ∫₀^φ |cos t| dt`, which is
`Sin[φ]` for `|φ| ≤ π/2` and not beyond it. Applied unconditionally the rule made
`EllipticE[2, 1]` answer `Sin[2] = 0.909297` where the value is `2 − Sin[2] = 1.090703` —
a jump of 0.18 against its own neighbour at `m = 1 − 10⁻¹⁸`. The reduction is now gated on
`ell_in_principal_strip`, so a symbolic amplitude stays symbolic rather than wrong, and
`N[]` routes it to Arb.

**Data structures.** `double` scalars through `carlson_rf_rd`; `acb_t` with the bridge's
accuracy ladder above machine precision. Both arities are registered ND kernels — unary
(`REG_U`) for the complete form, binary (`REG_B`) for the incomplete — so each rides the
packed/NDArray buffer and lowers inside `Compile[]`.

**Complexity / limits.** `O(1)` per element, ~6 ulp, no allocation. Real principal domain
only: `m > 1` for the complete form and `1 − m sin²φ < 0` for the incomplete one decline to
Arb. `Interval[]` threads by certified monotonicity (`E` decreasing in `m` below 1);
`Series` at `m = 0` uses the dedicated kernel `(π/2) Σ aₖ mᵏ/(1−2k)`.
