# M43 — §2.2.24 corpus (Problems 2301–2400) + deferred-class attempts

## Plan
Add Nasser Abbasi's §2.2.24 (2301–2400) to the DSolve corpus; measure; drive
coverage with general root-cause fixes. User opted in to ALSO attempt the deferred
research-grade classes (Abel Invariant Rational / special-fn Riccati), each behind
numeric self-verify (decline-safe). Invariants: **0 FAIL, 0 regression**. Leave
uncommitted. Section: 100 scalar (51 IVP), 0 systems; mixed 1st-order + 2nd-order.
Milestone M43, version 0.140 → 0.141.

## Tasks
### Stage 0 — ingest + measure
- [ ] 1. Generate `DSolve_test_status/DE_examples_2224.m` (converter, saved HTML)
- [ ] 2. Add `dsolve_corpus_2_2_24_tests` gate in `tests/CMakeLists.txt`
- [ ] 3. Build `dsolve_corpus_tests`; run full baseline (capture TSV)
- [ ] 4. Generate `reports/2.2.24.{md,tsv}` via `dsolve_corpus_report.py`
- [ ] 5. Triage the non-PASS list into buckets (already-solvable-but-broken vs deferred)

### Stage 1 — proven root-cause gap closure (0 FAIL, 0 regression)
- [x] 6. Baseline 90/100, 1 FAIL (2329), 9 U. Fixes landed:
      - **2329 FAIL** (correctness): IVP fit picked complex C[1]=I*Pi (wrong branch,
        satisfied IC not ODE). Fix in dsolve_common.c: collapse ConditionalExpression
        family to principal (C[_]->0) before the numeric candidate check + a scoped
        final numeric-verify gate on first-order scalar non-radical fits. → 2329 PASS.
      - **Homogeneous degree gate** (dsolve_homogeneous.c): decline exp-log-invert
        explicit Root of degree>=4 → clean implicit first integral (fixes unverifiable
        quintic-Root hang; helps LinearCoefficients recursion).
      - **Lagrange affine-ratio gate** (dsolve_lagrange.c + dsolve_lincoeff.c
        predicate): Lagrange defers the linear-coefficients class (spins on its
        integrating-factor solve) → LinearCoefficients solves it. → 2335 PASS.
      - **Riccati SeriesData gate** (dsolve_riccati.c): decline when linearised ODE
        returns a series (mapback hangs) — robustness.
- [x] 7. Regression sweep 1 (2329+homog+fit): 0 FAIL, 0 regression, +4 improved
      (221 2->1, 222 7->5, 226 5->3, 227 6->5). Sweep 2 (＋Lagrange) in progress.
      Final sweep with all 4 fixes pending.

### Stage 2 — deferred classes (attempted per user opt-in; decline-safe)
- [x] 8. 2a special-fn Riccati: CONFIRMED already solved (Airy 2348, Bessel 2358 PASS)
      via existing linearize+specialform. Residue 2349/2350/2351 (y'=e^{-t^2}+y^2 ->
      u''=e^{-t^2}u, non-elementary) + 2327 (symbolic Erf, needs 2nd-order operator
      factoring) = bounded declines.
- [x] 9. 2b Abel: ALL Abel-tagged cases in 2.2.24 already PASS (2353 via Chini/Abel;
      2335 via lincoeff; 2339/2343 via exact). No Abel gap in this section — no AIR
      work needed here.
- [x] 10. 2c `y=_G(x,y')`: 2352/2355/2356 genuinely non-elementary (e.g.
      y'=e^{-t}+Log[1+y^2]) — bounded declines (same class prior sections left).
      Residue after Stage 1+2: 8 U, all sympy=False, all bounded declines, 0 FAIL.

### Stage 3 — tracking + docs + version
- [ ] 11. Regenerate report; lower CMake baseline to new non-PASS count
- [ ] 12. Update STATUS.md + README row
- [ ] 13. Add M43 entry to DSOLVE_PLAN.md (+ method status lines if new method)
- [ ] 14. Changelog M43 in docs/spec/changelog/2026-09-07.md; new builtin docstring/docs if any
- [ ] 15. Bump src/version.h → 0.141
- [ ] 16. Gates: dsolve ctest suites, make check-c99, valgrind spot-check, REPL spot-checks
- [ ] 17. Rebuild code-review graph; write Review section here

## Review

**Outcome: §2.2.24 = 92/100, 0 FAIL, 0 crash, 0 regression.** Baseline was 90/100 with
1 FAIL. Final fix set is **two** root-cause fixes (an intermediate homogeneous degree-gate
and a Riccati SeriesData gate were tried and REVERTED — see below):

1. **2329 FAIL → PASS (correctness), `dsolve_common.c`.** The IVP constant-fitter collapsed a
   `ConditionalExpression` integer-family to its principal member (`C[_]→0`) only at the FINAL
   fit, not while CHOOSING among Solve's inverse branches. `Solve[Sinh[C]==0,C]` returns an odd
   `Iπ+2Ikπ` family and an even `2Ikπ` family; the candidate loop substituted a generic
   non-integer for the family index → both looked complex-nonzero → defaulted to `args[0]` (odd,
   wrong) → shipped `(1−t²)/2` (meets IC, fails ODE; correct `(t²−1)/2`). Fix: new
   `ds_collapse_principal` collapses each candidate to its principal member BEFORE the numeric
   `ds_branch_num_ok` check, plus a scoped final numeric-verify gate (first-order scalar,
   non-radical) dropping any confidently-wrong fit. Bonus: §2.2.1/2/6/7 each improved 1–2.

2. **2335 solved, `dsolve_lagrange.c` + `dsolve_lincoeff.c`.** Every affine-ratio equation also
   matches Lagrange as `y=x F(y')+G(y')` (F rational), but Lagrange's integrating-factor solve
   spins uninterruptibly (`TimeConstrained` cannot bound an inner `Integrate`). New predicate
   `dsolve_is_linear_coefficients_form` lets Lagrange defer the class to LinearCoefficients
   (solves in ~1 s via explicit Root). Genuine d'Alembert (polynomial F) unaffected.

**Dead ends (reverted, lesson learned):** first blamed the 2335 hang on homogeneous
verification and added a degree-≥4 `Root`→implicit gate in `dsolve_homogeneous.c` (broke
`t_stress_homogeneous`/`t_stress_lincoeff`, which verify explicit degree-4/5 Roots fine) and a
Riccati SeriesData gate (regressed 2346/2347, which pass via a concrete-coefficient
series-Riccati). Isolating the ACTUAL hang (Lagrange, via per-method timing) showed the
Lagrange gate ALONE fixes 2335 — the homogeneous Root form verifies fine (~1 s). **Lesson:
pin down the exact hanging method before adding gates; a net-improved count can hide a swap.**

**Residue 8** (all `sympy=False`, bounded declines, no wrong answers): 2304 (non-elementary IF
integral), 2327 (symbolic Erf-Riccati, needs symbolic 2nd-order operator factoring),
2349/2350/2351 (`y'=e^{−t²}+y²`, non-elementary), 2352/2355/2356 (`y=_G(x,y')`).

**Gates:** all DSolve corpus ctests 0 FAIL / at-or-below baseline (§2.2.1–24 + §2.1.2);
`dsolve_stress_tests` PASS; `dsolve_tests` SIGALRMs at the PRE-EXISTING `t_rischnorman` Abel
hang (49 tests pass first, 0 assertion failures); `make check-c99` PASS; valgrind: 0 serious
errors, 0 leaks in new code (only the documented pre-existing engine per-call baseline).
Left UNCOMMITTED per request.
