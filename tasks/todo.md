# DSolve M5 — NormalForm + Kovacic (Cases 1 & 2) + Frobenius/PowerSeries

Milestone M5: second-order linear variable-coefficient ODEs `y'' + P(x)y' + Q(x)y == 0`.
Confirmed scope: full slice + Kovacic Case 2 (Case 3 gated); Frobenius auto-fallback.

## Phase 1 — NormalForm (substrate + inspection builtin)  ✅ DONE
- [x] Add `dsolve_second_order_PQ(P, &Pc, &Qc)` to dsolve_common.{c,h} (shared P/Q extractor)
- [x] Add `dsolve_normal_form(Pc, Qc, xvar, &recovery)` → r = Pc^2/4 + Pc'/2 − Qc (guard ∫Pc/2)
- [x] New `dsolve_normalform.c`: `DSolve\`NormalForm[eqn,y,x]` → {r, recovery}
- [x] Refactor dsolve_specialform.c to use dsolve_second_order_PQ
- [x] Wire: dsolve_normalform_init in dsolve.c init; append .c to tests/CMakeLists COMMON_SRC
- [x] Build main + smoke-test NormalForm at REPL (Bessel/const/Airy all correct)

## Phase 2 — Frobenius / PowerSeries (auto last-resort fallback)  ✅ DONE
- [x] New `dsolve_frobenius.c`: classify x0 (ordinary/regular-singular/irregular)
- [x] Ordinary → two PowerSeries; regular-singular → Frobenius (indicial quadratic via analyze_roots)
- [x] Log solution for equal roots (d/ds derivative method); Puiseux x^(p/q) via x^r*SeriesData
- [x] FIXED substrate dsolve_verify_body: substitute D[body,{x,k}] directly (Derivative[Function[SeriesData]]=0 bug)
- [x] Register DSolve`PowerSeries + DSolve`FrobeniusSeries; cascade LAST (after autonomous)
- [x] Wire dsolve.c + CMakeLists; build + test (ordinary/regsing distinct+log/decline all verify O[x^k])
- [x] Repointed t_declines_unsupported to irregular Exp[1/x] (Sin[x]y now solves by design)

## Phase 3 — Kovacic (Cases 1 & 2)  ✅ DONE
- [x] New `dsolve_kovacic.c`: reduced form r via NormalForm; rationality gate (PolynomialQ)
- [x] Case 1a: Riccati ansatz ω'+ω²=r over pole structure (poles + ∞ poly), Solve coeff system
- [x] Case 1b: polynomial r via √r-poly-part + degree-bounded P (apparent singularities, e.g. x^2+3)
- [x] Case 2 (degree-2 algebraic): σ-ansatz D'+2σD=0, ω=(σ±√D)/2, NUMERIC back-sub verify
- [x] z1=Exp[∫ω] (guarded), z2=z1∫1/z1² ; recovery y=w·z; realify complex (Cosh+Sinh->E)
- [x] Register DSolve`Kovacic; cascade after specialform/before reduce_order; wire + CMakeLists
- [x] Fixed double-free (Denominator consumes rt) + numeric_verify (substitute C[1],C[2] basis)
- [x] Verified: 1+x^2→exp(x^2/2), x^2+3→x exp(x^2/2), Case2 x^(-1/4)exp(x^(3/2)/3), Bessel/Exp[1/x] decline

## Phase 4 — Tests, docs, gates
- [x] Unit tests in test_dsolve.c (14 new: NormalForm×3, Kovacic×6, Frobenius×5) — all pass
- [x] Stress: tests/test_dsolve_m5_stress.c (Kovacic fwd-generator 4 families + Frobenius) + add_test — passes in 6s
- [x] FIXED realify hang: Erfi[I z]->I Erf[z] rewrite (zero_test can't sample Erfi at complex irrational arg)
- [x] docs/spec/builtins/calculus.md + changelog/2026-08-31.md; flipped M5 in DSOLVE_PLAN.md
- [x] make check-c99 (exit 0); valgrind (no invalid access; leak == engine baseline); rebuilt graph
- [x] Updated tasks/lessons.md (3 lessons) + 3 durable memories + MEMORY.md index

## Review

M5 delivered three DSolve methods for 2nd-order linear variable-coefficient ODEs:

1. **NormalForm** (`dsolve_normalform.c` + substrate `dsolve_second_order_PQ` /
   `dsolve_normal_form`) — `{r, w}` reduction `z''==r z`; refactored
   `SpecialFunctionForm` onto the shared extractor.
2. **Kovacic** (`dsolve_kovacic.c`) — Riccati/undetermined-coefficient search
   (transparent + verifiable, unlike the exponent tables). Case 1 (rational ω),
   Case 1b (polynomial r apparent singularities, `x²+3 → x·exp(x²/2)`), Case 2
   (degree-2 algebraic, numeric-verified). Case 3 declines. `PolynomialQ` gate.
3. **Frobenius/PowerSeries** (`dsolve_frobenius.c`) — auto last-resort fallback;
   ordinary + regular-singular (incl. Log for equal roots) truncated SeriesData.

Substrate: `dsolve_verify_body` now substitutes `D[body,{x,k}]` directly (fixes
SeriesData verify). Three bugs found & fixed en route: double-free (Denominator
consumes rt), numeric_verify symbolic-in-C, zero_test hang on Erfi[I·irrational].

Verification: `dsolve_tests` (14 new unit tests) + `dsolve_m5_stress_tests` (50
generator cases) both green via ctest; `make check-c99` clean; valgrind shows no
invalid access and leak == Integrate/Solve/Simplify engine baseline (13.6KB, ≈
the DSolve-free baseline of 13.8KB).

Deferred (future M5b): Kovacic Case 3, ExactODE, OperatorFactor/DFactor,
Sturm–Liouville EigenvalueProblem; rational-r apparent singularities in Kovacic
(only polynomial-r apparent singularities handled); IVP-fitting of series
solutions. Known limitation: a very-messy raw coefficient (e.g. unexpanded
`(x+1/x)^2`) can make the substrate verify slow via zero_test on Erf/Erfi — clean
inputs and the internal r-simplification avoid it.
