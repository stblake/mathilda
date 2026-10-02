# Elliptic integrals: correctness, speed, Mathematica parity

Full review of `EllipticK` / `EllipticF` / `EllipticE` / `EllipticPi` (landed
v0.238). Measured against mpmath, scipy 1.18.1 and the Mathematica 14 figures in
the request. Values are right wherever the engine answers; two classes of
**silently wrong answer** and several speed/parity gaps are not.

Full diagnosis + baselines:
`~/.claude/plans/let-s-review-the-implementation-imperative-tide.md`.

## Phase 0 — close the measurement blind spot first
- [x] `tools/numeric_sweep.py`: add elliptic probes (feeds both
      `nd_surface_audit.py` and `nd_fastpath_sweep.py`, which today never
      exercise these heads on any surface — that is why all of this survived
      four sweeps).
- [x] `benchmarks/11-special-functions/special_functions.{m,py}`: bench rows for
      all five arities vs `scipy.special.ellipk/ellipe/ellipkinc/ellipeinc/elliprj`.
      `EllipticK` is declared in the `require[]` guard at `:20` with no row
      (`benchmarks/ABSENT.md:18`).
- [x] Record Phase-0 baselines in the changelog.

## Phase 1 — the silently wrong answers
- [x] **1a** `src/flint_num_bridge.c`: collapse the 11 duplicated kernel bodies
      into `nb_eval1/2/3` drivers taking an Arb function pointer, with the
      accuracy-retry loop in the driver: re-convert args inside the loop
      (`scalar_to_arb` rounds Rationals at `prec`), accept on
      `acb_rel_accuracy_bits >= outb` (absolute-radius test when the ball
      straddles zero), break immediately on `!acb_is_finite` (a pole must not
      retry), cap the doublings and return NULL rather than emit wrong digits.
      Fixes Zeta/PolyGamma/Hurwitz/Stieltjes/elliptic at once.
- [x] **1b** `elliptic.c:331-334`: `EllipticE[φ,1] -> Sin[φ]` is unguarded and
      wrong for `|Re φ| > π/2` (`EllipticE[2,1]` gives `0.909` for a true
      `1.0907`). Use `Sin[φ-kπ]+2k`, or guard and fall through to Arb.
- [x] **1c** `EllipticPi` + a visible `NDArray` must never stay unevaluated —
      `ndarray_delist_and_reeval` safety net (all three argument positions).
- [x] **1d** Inexactness + poles: exact branches must fire only on *exact* 0/1
      (`EllipticK[0.]` returns exact `1/2 Pi`; `EllipticPi[0., 1/2]` returns a
      *symbolic* `EllipticK[1/2]`); teach `ell_is_zero`/`ell_is_one` about
      `EXPR_MPFR`; poles return `ComplexInfinity`, not unevaluated.

## Phase 2 — machine accuracy and speed, complete integrals
- [x] **2a** AGM kernel for K and E (one pass yields both; validated 7 iterations
      max, K 1.42 ulp / E 5.12 ulp over `m ∈ [-10^6, 1-10^-12]`, vs today's
      61/106 ulp). Fix the dead `EC_RTOL` (`:80`, `(void)` at `:110`) and the
      `q < fabs(A)*0.01` criterion to NR's 0.0025 / 0.0015.
- [x] **2b** Re-measured: NOT NEEDED. The premise was a measurement error (see
      the correction below) -- K already beats scipy and E is within 1.23x, so
      the minimax-table work is unjustified. 2c instead went to the
      reciprocal-modulus route and the AGM was not needed either.
- [x] **2c** Complex arm for the unary descriptors, narrowly: real `m > 1` via
      the reciprocal-modulus identities (verified to 25+ digits — two *real* AGM
      passes, no complex AGM); genuinely complex `m` keeps declining to Arb.
      Closes the measured 20× one-bad-element cliff.
- [x] **2d** Scalar machine fast path at the six `ell_any_inexact` sites; return
      `EXPR_REAL`, not a 53-bit `EXPR_MPFR` (which does not pack and prints
      differently). K is 2180 ns/call today.

## Phase 3 — EllipticPi gets a machine kernel
- [x] **3a** `carlson_rc` + `carlson_rj` incl. Carlson's `p < 0` Cauchy-PV
      transformation.
- [x] **3b** `REG_B` kernel for complete `EllipticPi[n,m]` — buys
      `packed_aware`, the NDArray surface and `Compile[]` lowering for free.
- [x] **3c** 3-arg form: `ndarray_nary_kernel` (consumed by `Compile`'s
      `OP_KERNN`; the ND element-wise layer tops out at arity 2).
- [x] Fix the two stale "EXEMPT entry in the audits" cross-references
      (`src/ndkernels.c:543`, `special-functions.md:1882`) — no such entry exists.
      Also corrected the claim, in both places and in the docstring, that n > 1
      is a real Cauchy principal value: it is genuinely complex.

## Phase 4 — quasi-linear high precision (incomplete + third kind)
RESCOPED with the user to "try Landen, accept ~2x": implement the descending
transformation only, measure honestly, and keep it only where it actually beats
Arb. The theta-function route that could match Mathematica is 400-800 lines
whose hard part is certification (Newton produces no ball), and the Landen route
is estimated to land at ~1-3 s for F against Mathematica's 0.89 s. NOT STARTED.
- [ ] Landen/Gauss descending transformation for F, then E-incomplete, then Π.
      Today ~n² where complete K/E are ~n^1.4: F 4.30 s vs MMA 0.89 at 1e5
      digits, E-inc 22.2 s, Π 12.68 s vs 0.83 at 5e4.
- [ ] Guard rail (non-negotiable): the fast route answers only where it has been
      fuzz-verified against `acb_elliptic_*` at the same working precision,
      complex arguments included; otherwise it falls back to Arb.

## Phase 5 — symbolic parity
- [x] **5a** Exact values + oddness (all closed forms verified in-engine):
      `EllipticK[-1]`, `EllipticK[1/2]`, `EllipticPi[n,0]`, `EllipticPi[n,φ,0]`,
      `EllipticPi[1,m]`, the Infinity table, oddness in φ for F/E/Π.
- [x] **5b** Dedicated `Series` Taylor rules at `m = 0` for K and E (the generic
      path can never work — D at `m=0` is `0/0` → `Indeterminate` → the whole
      expansion is abandoned).
- [x] **5c** Interval: monotone rows for K (INC) and E (DEC) below `m = 1`; plus
      the general fix that unblocks F/Π in the φ-slot — negative rational
      exponents over a strictly-positive interval in `interval_power_pos_exp`.
- [x] **5d** Closed-form `m`/`n` derivatives replacing the inert forms; rewrite
      the pinning assertions in `tests/test_elliptic.c:167-175`; re-run the seven
      Mathematica corpus residuals.

## Phase 6 — tests, docs, audits
- [ ] `tests/test_elliptic.c`: Compile **value** tests (none exist for any
      elliptic head), the binary NDArray decline, exact values, oddness, Series,
      Interval, the `m=1` continuation, and a Defect-1 precision regression.
- [ ] Gates: `check-c99 check-messages check-packed-aware check-array-exactness
      check-nd-surfaces check-compile-coverage check-fastpath-sweep
      check-interval`, the CMake suite, valgrind.
- [ ] Docs: `special-functions.md:1835-1915`, the `Series` kernel list, the
      Interval section, the week's changelog (`2026-09-28.md`); repair the stale
      `1/Sqrt[x^3-x]` example in `integrate.c:928` / `calculus.md:1181`.
- [ ] Version bump + tag per phase.

## Review

### Landed and verified (v0.261)

Phase 0, all of Phase 1, and the accuracy half of Phase 2a. `elliptic_tests`
goes from 12 groups to 17, all green; spec + changelog + version bump done.

| fix | before | after |
|---|---|---|
| `N[EllipticK[1-10^-17], 20]` | 6 good digits of 20 | 20 of 20 |
| `N[Zeta[1 + 10^-20], 20]` | `2^66 + 1` | correct |
| `N[EllipticPi[-10^30, 1/2], 30]` | last 6 digits garbage | correct (ladder retries) |
| `EllipticE[2, 1]` | `Sin[2]` = 0.909 | symbolic; `N[]` = 1.0907 |
| `EllipticPi[NDArray[...], ...]` | unevaluated (×3 shapes) | evaluates |
| `EllipticK[0.]` | exact `Pi/2` | `1.5707963267948966` |
| `EllipticPi[0., 1/2]` | symbolic `EllipticK[1/2]` | `1.8540746773013719` |
| `EllipticK[SetPrecision[1,30]]` | unevaluated | `ComplexInfinity` |
| `EllipticPi[1., 0.5]` | unevaluated | `ComplexInfinity` |
| buffer vs scalar, K/E/F | 61 / 106 / 70 ulp | 2 / 7 / 2.5 ulp |
| `EllipticF[Pi/2 - 10^-8, 0.99]` | **2.7e-08** rel (1.2e8 ulp) | 1.7e-15 |

### CORRECTION: the speed baseline in the plan was wrong

`Timing` is **CPU** time and the ND element-wise kernels are parallel across
~9 cores (measured CPU/wall ratio 9.0-9.3 for `K`/`E`/`F`, 8.7 for `Gamma`),
while scipy is single-threaded wall time. Comparing the two was comparing
CPU-seconds to wall-seconds. Wall-clock, 10^6 elements, after the accuracy fixes:

| head | Mathilda wall | scipy | verdict |
|---|---|---|---|
| `EllipticK` | 11.3 ns/elt | 15.2 | **1.35x faster** |
| `EllipticE` | 18.3 | 14.9 | 1.23x slower |
| `EllipticF[φ,m]` | 14.7 | 170 | **11.6x faster** |
| `EllipticE[φ,m]` | 19.4 | 206 | **10.6x faster** |
| `EllipticPi` | ~40 000 | 316 | 128x slower (no kernel) |

So the "K 5.3x slower / E 12.9x slower" premise for Phase 2a-speed and the whole
of Phase 2b (minimax polynomial tables) does not exist. Per *core* Mathilda is
still ~7x less efficient than Cephes, so an AGM kernel would cut CPU cost and
close E's 1.23x — but it is an efficiency item, not the correctness-grade gap
the plan described. The high-precision measurements are unaffected (that path is
single-threaded, CPU/wall ratio 1.0).

### Two bugs the plan missed

- **φ → π/2 catastrophic cancellation** (fixed): `a = 1 - Sin[r]^2` where the
  value wanted is `Cos[r]^2`. Worst case 2.7e-08 relative, an eight-digit loss.
- **Elliptic scalars are not machine numbers** (open, Phase 2d):
  `Precision[EllipticK[0.5]]` is `15.9546` and `MachineNumberQ` is `False`,
  where `Zeta`/`Gamma` give `MachinePrecision`/`True` — so every
  elliptic-produced list is off the buffer for downstream consumers. The house
  precedent for the fix is `zeta.c:603-607`.
