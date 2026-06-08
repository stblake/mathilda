# Task: Close gaps between Mathematica's `Gamma` and ours (Tiers A + C)

Scope chosen by user: **Tiers A + C** (Gamma-internal closed forms + Gamma numerics).
Tier B (Erf/Erfc/PolyGamma sister functions, series) is explicitly out of scope.

## Tier A — self-contained Gamma closed forms (no new functions)

- [ ] **A1. `Gamma[n, z]` for positive integer `n ≥ 1`** → finite closed form
      Γ(n,z) = e^{-z} Σ_{k=0}^{n-1} (n-1)!/k! · z^k.
      - Replaces the current `n==1` special case (n=1 ⇒ just E^-z).
      - Runs for symbolic/exact z (so `Gamma[2,3]→4/E^3`, `Gamma[2,x]→(1+x)/E^x`);
        inexact z keeps falling to the mpfr numeric path (returns a Real).
      - Integer coefficients via GMP (BigInt-safe); guard n to a sane cap.
- [ ] **A2. `D[Gamma[a, z], z] = -z^(a-1) E^-z`** in `deriv.c`
      - `if (h == SYM_Gamma && n == 2)` block, LucasL-style drop-zero-arms.
      - dz arm: `-z^(a-1) Exp[-z] · D[z,x]`.
      - da arm: `Derivative[1,0][Gamma][a,z] · D[a,x]` (generic; PolyGamma is Tier B).
      - 1-arg `D[Gamma[z],z]` left to the generic `Derivative[1][Gamma][z]` path
        (closed form needs PolyGamma = Tier B).

## Tier C — Gamma numerics (complex / arbitrary precision)

- [ ] **C0. Complex-MPFR toolkit** (static, local to `gamma.c`)
      cmplx struct {mpfr_t re, im}; init/clear/set, add/sub/mul/div, exp, log,
      sin, cpow-via-exp(b·log). Pairs-of-mpfr_t (no MPC lib available).
- [ ] **C1. Arbitrary-precision complex `Gamma[z]`** via Spouge's approximation
      - Runtime-computable coefficients for any precision (unlike fixed Lanczos).
      - Reflection for Re(z) < 1/2. Replaces the "stays symbolic" limitation.
      - Machine complex keeps the fast double Lanczos path.
- [ ] **C2. Complex incomplete `Gamma[a, z]`** (machine + arbitrary precision)
      - Lower series γ(a,z)=z^a e^-z Σ z^n/(a)_{n+1} for small |z|;
        Lentz continued fraction for Γ(a,z) when Re(z)>0 & |z| large.
      - Γ(a,z) = Γ(a) − γ(a,z) using C1 for Γ(a) when needed.
      - Works at 53 bits (machine complex) and input precision (mpfr complex).

## Files
- [ ] `src/gamma.c` — A1, C0, C1, C2.
- [ ] `src/calculus/deriv.c` — A2.
- [ ] `tests/test_gamma.c` — new cases for each gap.
- [ ] `docs/spec/builtins/special-functions.md` — update Gamma section.
- [ ] `docs/spec/changelog/2026-06-08.md` — note the gap-closing.

## Verification
- [ ] `make -j` clean (`-std=c99 -Wall -Wextra`).
- [ ] `gamma_tests` + `deriv`-related tests pass.
- [ ] valgrind vs Sin[1.0] baseline — no new Mathilda-src leaks.
- [ ] REPL spot-checks vs Mathematica reference values.

## Review

**Status: complete.** All Tier A + C gaps closed and verified.

Files touched:
- `src/gamma.c` — A1 (`gamma_incomplete_int` finite closed form, replacing the
  `n==1` special case), C0 (`gcx` mpfr_t-pairs complex toolkit:
  add/sub/mul/div/exp/log/sin/pow), C1 (`gcx_gamma` Spouge + `gamma_mpfr_complex`),
  C2 (`gcx_lower_series` / `gcx_upper_cf` / `gcx_inc_gamma` /
  `gamma_mpfr_inc_complex`), shared `gamma_complex_result` (Real ≤53 bits, MPFR
  above).
- `src/calculus/deriv.c` — A2 (`SYM_Gamma && n==2` block).
- `tests/test_gamma.c` — new groups: arbitrary complex, integer-incomplete
  closed forms, derivatives, complex incomplete; symbolic group updated.
- `docs/spec/builtins/special-functions.md`, `docs/spec/changelog/2026-06-08.md`.

Verification:
- `make -j` clean under `-std=c99 -Wall -Wextra`; all `gamma_tests` groups pass;
  `deriv*`/`integrate_derivdivides` suites pass.
- C1 cross-checked against the independent machine Lanczos path (reflection +
  large args agree). C2 cross-checked via the recurrence
  `Γ(a+1,z) = a Γ(a,z) + z^a e^-z` to ~30 digits, and via the tiny-imaginary
  limit against the real `mpfr_gamma_inc` (both series and CF branches).
- Valgrind: gamma/deriv driver totals byte-for-byte identical to the `Sin[1.0]`
  baseline (12,800 B/400 + 3,720 B/56 — documented dyld/Accelerate noise); zero
  gamma.c / deriv.c / gcx_ frames in any lost stack.

Deliberate limits (Tier B, out of scope): `Gamma[1/2, z] = Sqrt[Pi] Erfc[..]`,
`D[Gamma[z], z] = Gamma[z] PolyGamma[0, z]`, and `Series[Gamma, ..]` still need
`Erf`/`Erfc`/`PolyGamma`, which are not yet implemented. Exact non-integer or
exact-complex incomplete forms stay symbolic, matching Mathematica.

Known cosmetic: arbitrary-precision complex prints `+ -0.498*I` (the MPFR
printer defect tracked in `MPFR_bugs.md`), not introduced here.
