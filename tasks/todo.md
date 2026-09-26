# Task: libm inf/nan must never leak into a Mathilda result

## Bug report
`Gamma[256.]` returned `Gamma[256.0]` (unevaluated). Broader class: elementary /
special functions on machine-real args whose true value over/underflows the IEEE
double range leak libm artifacts — `inf.0`, `-inf.0`, `nan`, a spurious `0.`, or
stay unevaluated. User directive: **`inf.0` and `nan` are not Mathilda
constructs; they must never be returned.**

## Root cause
`N[f[exact]]` already promotes machine overflow/underflow to an extended-exponent
MPFR real (numeric.c:1248 retry) — so every affected head already has a correct
MPFR path. The bug is that a machine-`EXPR_REAL` argument takes a *separate*
double-only fast path that calls libm directly and never promotes.

## Design (two complementary layers + enforcement)
- **L1 — universal safety net (invariant).** In `eval.c` at the builtin-return
  site: if a `ATTR_NUMERICFUNCTION` result carries a non-finite `EXPR_REAL`
  (`±inf`/`nan`, incl. inside `Complex[...]`), rewrite it to a real construct:
  `+inf → Overflow[]`, `-inf → -Overflow[]`, `nan → Indeterminate`,
  complex-nonfinite → `Overflow[]`/`Indeterminate`. Reads only the *result*
  (avoids the builtin-ownership-of-`res` problem); fires only on the cold
  non-finite path (zero hot-path cost). Guarantees the invariant for ALL heads.
- **L2 — per-function promotion (correctness).** Each affected head, in its
  machine-real branch (args still in scope), detects the degenerate double and
  recomputes via its existing MPFR kernel at 53-bit precision → extended-exponent
  `EXPR_MPFR`, matching `N[f[exact]]`. Shared helper for the common
  `get_approx`+libm pattern; inline MPFR for compound heads.
- **L3 — enforcement sweep.** Machine-input companion to `test_numeric_stress.c`:
  for each (head,arg), assert `f[x.]` is never inf/nan and, when the HP oracle is
  finite, `f[x.] ≈ N[f[x],30]`. Beyond-MPFR extremes must be `Overflow[]` (not
  inf.0). Ratcheting known-gaps like the existing harness.

## Affected heads (audit)
overflow→inf/nan: Exp, Power, Sinh, Cosh, Factorial, Binomial(nan), Erfi,
Gamma(2-arg), BernoulliB, EulerE, BarnesG, Hyperfactorial(nan), ExpIntegralEi.
unevaluated: Hypergeometric1F1/HypergeometricPFQ.
underflow→0: Exp, Power, Sech, Csch, Erfc, Gamma(2-arg), BesselK.
Already correct (mirror these): Gamma(1-arg, now fixed), Factorial2, Pochhammer,
Beta, LogGamma, BesselI, AiryAi/Bi, Sinh/CoshIntegral, Fibonacci, LucasL.

## Plan / checklist
- [x] Reproduce; locate machine-real paths; validate promotion mechanism on Gamma[1-arg].
- [x] L1: safety net in eval.c + `numeric_result_has_nonfinite`/`numeric_sanitize_nonfinite` (handles REAL + MPFR; not attribute-gated).
- [x] L2 helpers: `numeric_promote_machine_call` (re-eval with args→MPFR, 2-pass precision, demote-if-fits), `numeric_machine_real_or_promote`, `numeric_promote_result_if_degenerate`.
- [x] L2 wiring: Exp, Power, Sinh/Cosh/Sech/Csch, Factorial, Binomial, Gamma(2-arg, local keep-MPFR), Erfc, Erfi, ExpIntegralEi, BernoulliB, EulerE, BarnesG, Hyperfactorial, Hypergeometric1F1/PFQ, BesselK.
- [x] numericalize retry: also trigger on Overflow[]/Indeterminate and extended-range MPFR (fixes net×N interaction AND N[Exp[1000]] accuracy).
- [x] Tests: numeric_largearg (updated the case that codified the old inf.0), numeric_stress classifier reads Overflow[]/Underflow[] function form. Full battery 31/31; controls unchanged.
- [x] valgrind: no leaks from new code (only macOS objc/dyld startup noise, stacks all in system libs).
- [x] Full suite: 500 pass; 6 failures PROVEN pre-existing (byte-identical on stashed clean main): crc_corpus, dsolve_stress, intrischnorman, moebiusmu, primenu (documented probabilistic ECM), risch_rde_tower; dsolve_tests/dsolve_corpus slow/crashy pre-existing. Only numloop + numeric_largearg needed updating (they had codified the old inf.0).
- [ ] Commit + tag v0.209 — awaiting user go-ahead (on `main`, shared tree).

## Review
- Root cause: machine-`EXPR_REAL` args took a double-only fast path; overflow/underflow leaked `inf.0`/`nan`/`0.`/unevaluated. `N[f[exact]]` already promoted, so every head had a working MPFR path — the fix routes the machine branch into it.
- Two layers: (L1) evaluator net guarantees no inf/nan ever escapes (→ Overflow[]/-Overflow[]/Indeterminate); (L2) per-head promotion gives the correct finite extended-exponent value. Machine result when in double range (Binomial[1000.,500.]→MachineNumberQ True), MPFR when past it (Gamma[256.]), Overflow[] past MPFR range.
- Pre-existing, NOT introduced (verified against stashed main): `Sum[Exp[c. k],{k,1,n}]` with large machine c hangs / gives 0.0 (geometric-sum path). Out of scope; separate bug.
- 17 heads fixed + universal net. `check-*` audits untouched (no packed/compile surface change).
