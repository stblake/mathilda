# NSum — Numerical Summation (implementation plan)

Goal: implement `NSum[f, {i, imin, imax (, di)}, ...]` (plus multidimensional
`NSum[f, {i,...}, {j,...}, ...]`) in `src/numerical_calculus/nsum.c`, modelled on
Keiper 1992 "The N functions of Mathematica" §3–4 and Wolfram Language semantics.
Full state-of-the-art method suite, machine + arbitrary (MPFR) precision, no leaks,
extensive unit tests.

References: `Keiper 1992 - NFunctions.pdf` (SequenceLimit/Wynn ε §3, NSum/Euler–
Maclaurin/EulerSum §4); Cohen–Villegas–Zagier "Convergence Acceleration of
Alternating Series" (2000) for the AlternatingSigns method; Takahasi–Mori
double-exponential (exp-sinh) quadrature for the Euler–Maclaurin tail integral.

## Codebase facts that shaped the design (verified)
- **No general real-axis `NIntegrate`** — `SYM_NIntegrate` is only a *method name*
  for the Cauchy-contour path in NResidue/ND. The Euler–Maclaurin tail integral
  must come from a NEW self-contained quadrature (chosen) or symbolic `Integrate`.
- **Acceleration kernels are `static` in `nlimit.c`** (`nl_richardson_machine`,
  `nl_wynn_machine`, `nl_richardson_mpfr`, `nl_wynn_mpfr`). Extract into a shared
  module so NSum and NLimit share one implementation (no duplication).
- `BernoulliB[n]` exists (machine + MPFR) — usable for EM corrections.
- `D` will be applied **iteratively** (n times) for f^(2k-1); do not assume a
  `D[f,{x,n}]` nth-order form.
- Existing `quadrature.{c,h}` is a periodic-trapezoidal *contour* engine — not a
  half-line integrator; the DE quadrature is a new sibling module.
- Symbols already interned: AccuracyGoal, PrecisionGoal, WorkingPrecision,
  WynnDegree, MachinePrecision, EulerSum, SequenceLimit, NIntegrate. Need to add:
  NSum, NSumTerms, NSumExtraTerms, VerifyConvergence, EulerMaclaurin,
  AlternatingSigns (and reuse EulerSum/SequenceLimit as method aliases).
- Localization (Block-style) reuses the `NlBind` snapshot/set/restore pattern from
  `nlimit.c` (save OwnValues+attrs, bind var to a numeric value, eval, restore).

## New / modified files
New:
- `src/numerical_calculus/seqaccel.{c,h}` — shared Wynn-ε + Richardson kernels
  (machine `double _Complex` and MPFR real/imag pairs), pure numeric (no Expr).
- `src/numerical_calculus/dequad.{c,h}` — double-exponential (exp-sinh) quadrature
  for ∫_N^∞ f(x) dx, machine + MPFR, complex-valued integrand, sample-callback
  style mirroring `quadrature.h`. Reusable groundwork for a future real NIntegrate.
- `src/numerical_calculus/nsum.{c,h}` — `builtin_nsum`, option parsing, method
  dispatch, localization, multidimensional recursion.
- `tests/test_nsum.c` (primary), `tests/test_seqaccel.c`, `tests/test_dequad.c`.

Modified:
- `src/numerical_calculus/nlimit.c` — call shared `seqaccel_*` kernels (behavior-
  preserving refactor; guarded by `nlimit_tests` before/after).
- `src/core.c` — declare + call `nsum_init()` next to `nlimit_init()`.
- `src/sym_names.{c,h}` — intern NSum, NSumTerms, NSumExtraTerms, VerifyConvergence,
  EulerMaclaurin, AlternatingSigns.
- `src/info.c` — NSum docstring (terse; per `symtab_set_docstring`).
- `tests/CMakeLists.txt` — add the 3 new `src/*.c` to COMMON_SRC (~line 277) and the
  new test executables (~line 349 block).
- `docs/spec/builtins/<numerical-calculus page>.md` + weekly changelog
  `docs/spec/changelog/2026-06-08.md` (Monday of current ISO week).

Attributes: `NSum` = `ATTR_HOLDALL | ATTR_PROTECTED` (HoldAll so the summand and the
iterator specs are held; we localize the index ourselves). Not Listable.

## Methods (the algorithms)

1. **Direct** — finite sum with a small term count (≤ DIRECT_MAX, e.g. ~ a few
   hundred or `2*NSumTerms + NSumExtraTerms`): bind i, numericalize f(i), accumulate
   at WorkingPrecision. Exact small finite sums.

2. **"WynnEpsilon" / "SequenceLimit"** (general fallback) — build partial sums
   S_N..S_{N+M} (M ≈ `NSumExtraTerms`, ≥ 2·WynnDegree+2) starting after `NSumTerms`
   leading terms; extrapolate the partial-sum sequence with `seqaccel_wynn_*`.
   Robust for alternating & many series; complex-capable.

3. **"EulerMaclaurin" / "Integrate"** (best for monotone, slowly-converging
   positive series) — sum the first `NSumTerms` terms explicitly; let N = imin +
   NSumTerms. Tail:
   `Σ_{i=N}^∞ f(i) ≈ ∫_N^∞ f dx + f(N)/2 − Σ_{k≥1} B_{2k}/(2k)! · f^(2k-1)(N)`.
   - ∫_N^∞ via the new `dequad` exp-sinh quadrature (refine h until PrecisionGoal/
     AccuracyGoal; doubling-comparison error estimate).
   - f^(2k-1)(N) by iterated symbolic `D` then numericalize at N; `B_{2k}` via
     `BernoulliB`. Asymptotic series: cap K, stop at the smallest term.
   - Fallbacks: if `D` blows up / fails, or quadrature won't converge → fall back to
     WynnEpsilon. (Matches Keiper: Integrate then degrade.)

4. **"AlternatingSigns"** (alternating series, genuinely SOTA) — Cohen–Villegas–
   Zagier algorithm on the magnitudes a_k (with f(imin+k) = (-1)^k a_k). For p-bit
   WorkingPrecision use n ≈ p/1.76 + guard terms; ~1.76·n bits from n terms.
   Machine + MPFR; complex a_k allowed.

5. **Automatic dispatch** — sample first ~10–15 terms numerically: strictly
   alternating signs → AlternatingSigns; positive & monotonically decreasing →
   EulerMaclaurin; else → WynnEpsilon. Finite small → Direct.

6. **VerifyConvergence** (infinite limits only, default True) — ratio test on
   |a_{k+1}/a_k|: >1 ⇒ divergent (return unevaluated + `NSum::div` message);
   ≈1 inconclusive (proceed, per Keiper); <1 proceed. `VerifyConvergence->False`
   skips the test (faster).

7. **Step `di`** — reindex i = imin + di·j (j = 0,1,...; jmax = (imax−imin)/di or ∞);
   summand g(j) = f(imin + di·j); then a standard single-var sum in j.

8. **Multidimensional** `NSum[f, s1, s2, ...]` — recursive: the outer single-var
   driver's "summand" is the held expression `NSum[f, s2, ..., opts]`. Binding the
   outer index as an OwnValue (Block) lets inner bounds depend on it (e.g.
   `{k,1,n}`); evaluating the inner NSum returns a number. Falls out of HoldAll +
   localization; options propagate to the inner call.

## Options & defaults
WorkingPrecision (MachinePrecision → machine path; digits → MPFR, bits =
digits_to_bits), NSumTerms (15), NSumExtraTerms (precision-scaled; ~12 machine),
WynnDegree (1), Method (Automatic | EulerMaclaurin/Integrate |
WynnEpsilon/SequenceLimit | AlternatingSigns | EulerSum), VerifyConvergence (True),
AccuracyGoal (Infinity), PrecisionGoal (Automatic = WP−10 digits). Unknown option ⇒
return NULL (unevaluated), matching the NLimit pattern.

## Implementation phases (each ends with its tests green + valgrind clean)
- **P0** Wire skeleton: `nsum.{c,h}`, `nsum_init`, sym_names, CMake COMMON_SRC,
  docstring stub. `NSum[1/i^2,{i,1,5}]` via Direct only. Build + smoke REPL.
- **P1** Extract `seqaccel.{c,h}`; refactor `nlimit.c` to use it; `nlimit_tests`
  pass unchanged (regression guard). Add `test_seqaccel.c`.
- **P2** WynnEpsilon method (partial sums + seqaccel) for infinite sums; machine +
  MPFR; complex. Tests: geometric, alternating, Keiper π/4 example.
- **P3** `dequad.{c,h}` exp-sinh quadrature (machine + MPFR), standalone
  `test_dequad.c` against known integrals (∫_1^∞ 1/x^2 = 1, etc.).
- **P4** EulerMaclaurin method (explicit terms + dequad tail + Bernoulli +
  iterated D corrections + asymptotic stop). Tests: Zeta(11/10), Zeta(51/50),
  Σ1/i^2 from 100..10^6.
- **P5** AlternatingSigns (CVZ) machine + MPFR; high-precision alternating test
  (Keiper `(-1)^x/(1+(x-12)^2)`, WP→30).
- **P6** Automatic dispatch + VerifyConvergence ratio test (divergence message).
- **P7** Step `di` + multidimensional recursion (Keiper examples incl. dependent
  inner bound `{k,1,n}`).
- **P8** Docs (spec page + changelog 2026-06-08), final valgrind sweep, full
  scoped test run.

## Testing (tests/test_nsum.c — compare in-language, never parse printed reals)
Use the established in-language comparison pattern (eval `NSum[...] - exact` and
check |result| < tol via the evaluator), e.g. against `Exp[-5]`, `Zeta[11/10]`,
`Zeta[51/50]`, `N[Pi/4]`, `N[Pi^2/6]`, `Sinh[Pi]/Pi`. Cover: finite small/large,
infinite, alternating (machine + WP→30), complex summand
(`Log[x]/x^(2+2I)`), step di (`1/2^i,{i,0,∞,2}` = 4/3), multidim (independent and
dependent inner bound), Fibonacci reciprocal, each explicit Method, VerifyConvergence
True/False, divergence (`2^i` ⇒ message + unevaluated/ComplexInfinity-style),
WorkingPrecision precision check. Add `nsum.c`/`seqaccel.c`/`dequad.c` to
`tests/CMakeLists.txt` COMMON_SRC (else `*_tests` fail to link `nsum_init`).

## Memory & standards
Builtin ownership contract (return new Expr; never `expr_free(res)`; NULL-out reused
args). `mpfr_clear` every temp; free Sr/Si arrays. Strict C99 (no `M_PI`/POSIX
without guards). Valgrind diff against a `Sin[1.0]` baseline to ignore dyld noise.

## Risks / mitigations
- nlimit refactor regression → guard with `nlimit_tests` pre/post.
- Symbolic high-order `D` explosion → cap order, stop at smallest EM term, fall back
  to Wynn.
- exp-sinh on slowly-decaying integrands (e.g. ζ near 1) → bounded refinement +
  AccuracyGoal/PrecisionGoal termination; degrade to Wynn if non-convergent.
- CVZ assumes near-monotone magnitudes → only used when signs strictly alternate.

## Review (completed 2026-06-13)
All phases landed and verified.

- **seqaccel.{c,h}** — extracted Wynn-ε + Richardson kernels (machine + MPFR);
  `nlimit.c` refactored onto them; `nlimit_tests` pass unchanged.
- **dequad.{c,h}** — exp-sinh half-line quadrature (machine + MPFR, complex).
  Two subtleties resolved: (1) MPFR `t`-range must scale generously for slow
  algebraic tails (ζ near 1) — `tmax = 6.5 + 1.1 ln(bits)`; (2) the tail
  truncation must stop at the term *minimum* and break only on a ≥4× blow-up,
  so an integrand hitting a roundoff noise floor (e.g. `Log[a/b]` evaluated as
  `Log a − Log b`) doesn't get its noise amplified by the growing weight.
- **nsum.{c,h}** — Direct / WynnEpsilon / EulerMaclaurin / AlternatingSigns
  (CVZ) / Automatic; step `di`; multidimensional via inner-NSum; large finite
  via difference-of-tails; VerifyConvergence ratio test → ComplexInfinity. EM
  correction loop capped by derivative **expression size** so `1/Fibonacci[i]`
  can't explode.
- Registered (`core.c`, `sym_names.{c,h}`), docstring (`info.c`), docs
  (`calculus.md` + changelog `2026-06-08.md`), tests (`test_nsum.c`, 13 groups),
  CMake wired. Build clean `-Wall -Wextra`; `leaks` shows no NSum-attributable
  leaks.

Accuracy vs Keiper reference values: ζ(2)=2e-16, ζ(11/10)=2e-15 (machine),
ζ(51/50)=2.6e-27 (WP40), Wallis log-sum=4e-8, alternating peaked sum to 27
digits (WP30, AlternatingSigns), complex `Log[x]/x^(2+2I)` exact, multidim
1.14434 / 0.770188 exact, large finite Σ_{100}^{10^6} 1/i²=0.0100492.

Known soft spot: a fast-rising-then-falling alternating summand like
`(-5)^i/i!` needs `NSumTerms -> 25` to clear the peak (default 15 → 9e-11);
matches Keiper's NSumTerms guidance.
