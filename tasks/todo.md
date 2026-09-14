# M44 — §2.2.25 corpus (Problems 2401–2500), push for 100/100

## Plan
Add Nasser Abbasi's §2.2.25 (2401–2500) to the DSolve corpus; measure; drive
coverage toward **100/100** with general root-cause fixes only. Section: 100
scalar ODEs (~32 IVP), 0 systems; 2nd-order-linear-heavy (Airy/Hermite/Legendre/
Chebyshev/Laguerre/Gegenbauer/Bessel, Euler–Cauchy, regular-singular Frobenius,
Emden–Fowler nonlinear, const-coeff nonhomogeneous incl. one `sec(t)` forcing).
Invariants: **0 FAIL (0 wrong answers), 0 regression**. Real methods only — no
hacks/overfit; a genuinely research-grade case is a *documented decline*, never a
faked pass. Leave uncommitted. Milestone M44, version 0.141 → 0.142.
Probe of 11 anticipated-hard cases: 10 already solve, 1 declines (#9 parabolic
cylinder + cos forcing).

## Tasks
### Stage 0 — ingest + measure
- [x] 1. Generated `DSolve_test_status/DE_examples_2225.m` (100 recs, 32 IVP, 0 sys)
- [x] 2. Spot-check found + fixed THREE converter bugs (all root-caused in
      `tools/latex_ode_to_mathilda.py`):
      - **indvar juxtaposition**: `t` glued to digit/letter (`ty`,`2t`) missed →
        defaulted to `x` (8 recs in §2.2.25; also 6 latent in §2.2.24). Fallback
        re-scan on implicit-mult-separated body.
      - **autonomous parameter as indvar**: `y'=k(a-y)(b-y)` picked param `a`
        (ABORT); now unique-candidate rule → fresh `x` (#2498; also §2.2.24-2327).
      - **piecewise sentinel shred + `\begin{array}[]{cc}` optional arg**: `PWFORCE0`
        split by implicit-mult → garbage; now whitespace-tolerant expand + optional
        `[]` strip (#2487 → clean `Piecewise[...]`).
      Byte-identity: §2.2.20–23 IDENTICAL; §2.2.24 improved by 7 correct recs
      (regenerated in place). §2.2.1–19/§2.1.2 are tex4ht-sourced, out of scope.
- [x] 3. Built harness. Baselines: **§2.2.25 = 95/100** (0 FAIL/crash, 5 UNEVAL);
      **§2.2.24 re-measured = 91/100** (converter fix: 2327 U→P, but 2353/2357 revealed
      as honest UNEVAL — M43's 92 included 2 false-passes on t-as-constant equations).
- [ ] 4. Generate `reports/2.2.25.{md,tsv}` via `dsolve_corpus_report.py`
- [ ] 5. Triage non-PASS: (a) converter artifact (b) engine gap (c) FAIL (d) research-grade

### Stage 1 — root-cause gap closure toward 100/100 (0 FAIL, 0 regression)
- [x] 6. No FAILs to fix — §2.2.25 had 0 FAIL, 0 crash from the baseline.
- [x] 7. Converter artifacts closed (3 root-cause fixes, see item 2).
- [x] 8. **VoP fractional-power Simplify-hang fix** (`dsolve_common.c`): the final
      `ds_simplify(yp)` hangs on a fractional-power×exp answer (`Simplify[t^(5/2)E^(-2t)]`
      itself spins); the body is already auto-eval clean, so skip Simplify for a
      fractional-power answer (`ds_has_fractional_power`). No nested TimeConstrained.
      → 2406 solves; §2.2.25 95→**96/100**.
- [x] 9. Probed all residue. §2.2.25 residue 4 (2409/2410/2444/2477) are ALL
      sympySolved=False (parabolic-cylinder, clean-basis-VoP, slow transcend series,
      non-elem IF) — genuine research-grade/nonelementary, documented declines, no
      hacks. §2.2.24 residue 9 (2353 Abel, 2357 nonelem Bernoulli IF, + M43's set).
- [~] 10. Re-measured after fix (§2.2.25=96). Full §2.2.x regression sweep RUNNING.

### Stage 2 — register + docs + version
- [x] 11. Saved `reports/2.2.25.{tsv,md}` + re-saved `reports/2.2.24.{tsv,md}`; added
      `dsolve_corpus_2_2_25_tests` gate (baseline 4); §2.2.24 gate 8→9.
- [x] 12. STATUS.md §2.2.25 block + §2.2.24 M44 correction + M44 milestone-log line;
      README.md row + §2.2.24 note.
- [x] 13. `src/version.h` 0.141→0.142; DSOLVE_PLAN.md M44 entry; changelog M44 section.

### Stage 3 — verification
- [x] 14a. `ctest -R dsolve_corpus_2_2` — **100% (25/25) passed, 0 regression** across
      §2.2.1–25 (incl. new §2.2.25 baseline 4, updated §2.2.24 baseline 9). §2.2.25 gate
      re-confirmed on the final rebuilt binary.
- [x] 14b. §2.1.2 gate (**1000 cases**) PASSED, 0 regression; dsolve_stress_tests +
      dsolve_m34_stress_tests PASSED. **dsolve_tests SIGALRM** — CONFIRMED **pre-existing**:
      reverted my change to clean HEAD (M43), rebuilt, ran → also exits 142 (cumulative
      >120s `alarm(120)` on the Risch–Norman-heavy suite, unrelated to VoP). Fix restored.
- [x] 15a. `make check-c99` PASSES. Rebuilt `./Mathilda`; REPL spot-checks (2406 solves,
      2487 piecewise, symbolic Bessel/Chebyshev) all good.
- [x] 15b. valgrind: VoP fix is **leak-neutral** — fractional-power path (2406) and
      integer path both lose the identical 13,496 B (the documented inherited
      Integrate/Solve-engine per-call leak); the read-only helper adds no allocations.
- [x] 16. Staff-engineer self-review done (diff is exactly helper + gated Simplify, no
      debug residue; read-only traversal, interned-pointer compare). Graph refresh skipped:
      incremental diff is vs HEAD~1, would not index uncommitted work; left uncommitted.

## Review

**Outcome.** M44 adds §2.2.25 (Problems 2401–2500) to the DSolve corpus at **96/100,
0 FAIL, 0 crash**, with **0 regression** across every sibling gate (§2.1.2 + §2.2.1–24).

**Root-cause fixes (no overfit, no hacks):**
- **`dsolve_variation_of_parameters` fractional-power Simplify hang** (`dsolve_common.c`,
  general): resonant fractional forcing closes to an elementary answer but the final
  `ds_simplify` spun on `t^(p/q) E^(a t)`; skip it for a fractional-power answer (already
  auto-eval clean). Fixes 2406 and the whole class. Leak-neutral, C99-clean.
- **Three converter root-cause bugs** (`latex_ode_to_mathilda.py`): indvar juxtaposition,
  autonomous-parameter-as-indvar, piecewise-sentinel/optional-arg. §2.2.20–23 byte-identical;
  §2.2.24 regenerated (7 latent wrong-equation records corrected — a genuine integrity fix).

**Honest ceiling.** The §2.2.25 residue of 4 (2409/2410/2444/2477) are all `sympySolved=False`
— parabolic-cylinder + forcing, clean-basis VoP, slow transcendental series, non-elementary
integrating factor — genuinely beyond current CAS reach (SymPy fails them too). §2.2.24 is a
truer 91/100 (M43's 92 counted 2 mistranscribed-equation passes; 2327 now genuinely solves).

**Deliverables.** New `DE_examples_2225.m` + `reports/2.2.25.{md,tsv}`; regenerated
`DE_examples_2224.m` + its reports; gates `dsolve_corpus_2_2_25_tests` (baseline 4) and
`dsolve_corpus_2_2_24_tests` (8→9); STATUS.md/README.md/DSOLVE_PLAN.md/changelog updated;
version 0.141→0.142. Left uncommitted per request.

**Pre-existing (not addressed, not caused here):** dsolve_tests SIGALRM (Risch–Norman suite
>120s), and the §2.2.24-2353 Abel / research-grade residue (M13-deferred).

## Review
_(to be filled in on completion)_
