# DSolve §2.1.3 — solve the 12000.org "unsolved" ODE set

Plan: `/Users/user/.claude/plans/let-s-devise-algorithmic-methods-glowing-toast.md`
Baseline (157 cases, current binary): 40 SOLVED · 100 UNSOLVED · 14 TIMEOUT · 3 CRASH.
Corpus harness + per-case verdicts in the session scratchpad.

## Work items (P0–P4, comprehensive)

- [x] **P0a** Fixed the `#144` segfault — heap double-free in Risch–Norman
      `enumerate_monomials` cap-cleanup (`intrischnorman.c`, NULL-out-before-free).
      ASAN-clean, 18-case generalized suite + negative control + unit regression.
- [~] **P0b** Anti-hang: timeouts 14 → 5 (crashes 3 → 0 real; #49/#50 were
      corpus-conversion syntax errors, not crashes). Fixed:
      #1,15,16,27,52,58,59,71,76,94 (+#11 regression repaired).
      - SeparableReduced: gate to F rational-in-y (`dsolve_sepreduced.c`).
      - LinearFirstOrderSystem: reject nonlinear A(Y) (`dsolve_linsys.c`).
      - Integrate linearity + Expand[TrigReduce] (`integrate.c`) — fixes the
        trig-sum/power hang class; μ = PowerExpand[Simplify[Exp[∫p]]]
        (`dsolve_common.c`). Both with 15–20 case anti-overfit suites (100%).
      - Remaining 5 (all pre-existing): #123,#132 (triangular-glue verify) → P1;
        #19,#26,#114 (Lie/Exact/Abel Expand blowup on transcendental-of-y) → P3.
- [x] **P1a** Higher-order coupled systems: `dsolve_sys_reduce.c` (state
      augmentation → `Y'=AY+b` → `dsolve_linsys_assemble`). 2nd-order systems
      **3 → 50 SOLVED** of 58. Perf: size/content-adaptive `dsolve_linsys_tidy`
      (Expand for large/mixed-spectrum bodies) + Expand the VoP integrand (forced
      systems). Suite green; P1 unit tests added.
      Remaining: #37/38/91/132 (slow forced mixed-order, ~17s), #40/#43 (symbolic
      coefficients → symbolic Jordan), #151 (singular leading matrix → P1b),
      #72 (corpus mis-conversion).
- [ ] **P1b** Operator-determinant elimination fallback (singular leading
      matrix / DAE #151) inside `dsolve_sys_reduce.c`.
- [x] **P2a** Const-A + variable-forcing first-order systems: dropped the
      `b(x)`-variable rejection in `dsolve_linsys_solve` (VoP handles it).
      Solves #7,10,67,156.
- [x] **P2b** Wider variable-coeff systems. (1) Fixed scalar-factor miss:
      `dsolve_linsys_extract_Ab` now admits a t-dependent leading coeff (`t x'+y=0`)
      → #33,35,56,73,74,78,109,110. (2) New `dsolve_linsys_commutative.c` for the
      2×2 commutative class A=a(t)I+b(t)K0 (rotation/hyperbolic/nilpotent via
      Φ=Exp[∫a]Exp[(∫b)K0]) → #48,69,80,126,127. Corpus **101 → 114 SOLVED**,
      zero regressions (#47 3×3-symbolic solves in ~15s, borderline). 13-case
      anti-overfit + `t_system_varcoeff`; DSOLVE_PLAN §1e + changelog.
      Remaining var-coeff: #3,68,108,128 (genuinely non-commutative), #30/31
      (arbitrary f,g → non-elementary ∫), #53/54/55 (corpus mis-conversion).
- [x] **P3** First-order linearizable tail: new `dsolve_linearizable.c`
      (`u=φ(y)` for φ∈Log/Exp/Sin/Cos/Tan → linear/Bernoulli in u; runs before
      Exact). Solves #15,16,26,27,76,77,114,115 (+ the log family stays solved).
      Corpus **93 → 101 SOLVED**, timeouts 10 → 8, zero regressions. 17-case
      anti-overfit suite + `t_linearizable_first_order`; DSOLVE_PLAN §1a + changelog.
      Remaining tail: #18 (u=Log[1+y], not in table), #143 (solvable-for-x),
      #19/#144 (sqrt / rational-symmetry — not transcendental-of-y).
- [x] **P4** Higher-order singletons: `dsolve_linear_normalize` (`dsolve_common.c`)
      — clear denominators (rational RHS → Euler form) + divide by coefficient GCD
      (common factor → const-coeff), called by const-coeff/euler/undetcoeff,
      depth-gated to the outermost call (pinned entries now count a level) so
      OperatorFactor recursions stay fast. Solves #95 (Euler `x³y'''-24y=24x`),
      #96 (const-coeff `y'''+2y''-y'-2y=1/x` → ExpIntegralEi). #12 (symbolic-order
      Bessel with arbitrary a(m)) left. 12-case anti-overfit + `t_linear_coeff_normalization`;
      DSolve/stress/m5 suites green; c99 clean; DSOLVE_PLAN §1b + changelog.
- [ ] **P1b** Operator-determinant elimination (singular leading matrix #151;
      symbolic-coefficient systems #40/#43). (NOT DONE.)
- [x] **P5 (nonlinear systems)** New `dsolve_autosys.c` — 2D autonomous (incl.
      nonlinear) via phase-plane: orbit `dy/dx=g/f` → reconstruct `x(t)` from
      `x'=f(x,Y(x))`, orbit constant renumbered to C[2]. Solves #52,58,59,64,65,
      135,139. Gated to rational orbits (radical → decline; also sidesteps a
      pre-existing NULL-deref in radical Risch-Norman `risch_squarefree_t`).
      Corpus **116 → 123 SOLVED**, zero regressions, ASAN-clean. 10-case
      anti-overfit + `t_system_autonomous`; DSOLVE_PLAN §1e + changelog.

- [x] **P6 (robustness)** Fixed the pre-existing NULL-deref in
      `risch_squarefree_t` (`risch_canonical.c`): the initial `c=p/g` / `ppg=p'/g`
      field-division chain went unguarded, so a NULL (on a non-polynomial-in-`t`,
      e.g. radical, argument) flowed into `rc_ddt`/`rc_sub_expand` → function node
      with a NULL arg → `evaluate()` SIGSEGV. Now declines cleanly (guards mirror
      the loop body). The `AutonomousSystem` radical gate is consequently
      downgraded to a pure **speed** gate (crash no longer possible). 10-case
      radical-input suite `test_squarefree_radical_no_crash` in
      `test_risch_canonical.c` (negative control: 5/10 SIGSEGV on reverted code);
      `check-c99` clean; DSolve/stress/risch suites green; no corpus regression.

## Verification
- [x] Re-run the 157-case corpus harness. **Corpus: 40 → 93 SOLVED**,
      100 → 52 UNSOLVED, 14 → 10 TIMEOUT, 3 → 0 CRASH (2 "ERR" are corpus
      conversion syntax errors, not crashes). **Zero regressions.**
- [x] `tests/test_dsolve.c`: added `t_rischnorman_enum_cap_no_crash`,
      `t_trig_coeff_linear_first_order`, `t_system_higher_order`;
      `test_integrals.c`: trig-linearity assertions. All suites green.
- [x] ASAN: fixed 2 real leaks in new code (`dsolve_sys_reduce` `t` arg,
      `try_linearity` `sum`); residual forced-system leak = documented pre-existing
      FLINT `rat_canon` 56 B/call.
- [ ] `make check-c99` (not yet run this session); valgrind full pass.

## Review (checkpoint after P0 + P1 + P2a)

Delivered, all verified (2nd-order systems **3 → 50 SOLVED**):
- P0a segfault (#144 Risch-Norman double-free); P0b integration linearity +
  method gates + μ canonicalization (timeouts 14 → 10, all remaining pre-existing);
  P1 higher-order coupled systems (`dsolve_sys_reduce.c` state augmentation +
  adaptive tidy + VoP-integrand expand); P2a const-A variable forcing.
- Anti-overfit suites (100%): 18-case Risch-Norman, 20-case trig integrals,
  15-case trig-coeff linear; DSOLVE_PLAN.md §1e + weekly changelog updated.

Not done (clearly scoped follow-ups): P2b (variable-coefficient first-order
systems — needs symbolic MatrixExp, ~18 cases), P3 (first-order linearizable
scalar tail, ~12 cases), P4 (higher-order singletons), P1b (elimination for
singular-leading-matrix / symbolic-coefficient systems).

---

# DSolve Tier 1 — bug fixes via general methods

Plan: `/Users/user/.claude/plans/the-following-are-bugs-rosy-pelican.md`

## Tier 1 work items

- [x] **T1** First-order integrating-factor returns integral form
      (`dsolve_common.c` `dsolve_linear_factor_solve`) — fixes A1 (iterlim →
      Erf closed form), A2 (wrong answer → integral form), A3 (Bernoulli
      fall-through). Also: undetcoeff now verifies L[y_p]==T decidably; verify
      guard against zero_test E^a-E^b false-negative. Suites green.
- [x] **T2** Legendre recognizer in `dsolve_specialform.c` + `LegendreP` D-rule
      in `deriv.c` — fixes B5, B-legendre2, C4 (Riccati handoff). Also: Riccati
      back-map skips Simplify on special-function bodies (was a hang). Suites green.
- [ ] **T3** Exponential-argument Bessel branch in `dsolve_specialform.c` —
      fixes B10.
- [ ] **T4** Missing-lower-derivatives order reduction (new method) — fixes D19.
- [ ] **T5** Symmetric-square recognizer (new method) — fixes E16 (Airy),
      E17 (Bessel).
- [ ] **T6** Anti-hang guards: Kovacic numeric-coeff guard + operator_factor
      complexity gate — B7, E18 stop hanging.
- [ ] **T7** Gauss-2F1 symbolic-c gate + zero-test special-function decline —
      fixes B8.

## Verification
- [ ] Per-case REPL probes (with timeouts) match expected forms; no hang, no
      `$IterationLimit`.
- [ ] `tests/test_dsolve*.c` regression green; add focused cases.
- [ ] `make check-c99`; docs/spec + weekly changelog; docstrings/attributes.

## Review — Tier 1 COMPLETE

All 7 work items landed; every reported Tier-1 case is fixed and all DSolve
suites (`dsolve_tests`, `dsolve_stress_tests`, `dsolve_m5_stress_tests`) plus
adjacent suites (`zero_test`, `trigexp_zero`, `fullsimplify`, `reduce`, `solve`)
are green. `make check-c99` clean. Final verification: 8 headline cases return
correct closed forms with **0** `$IterationLimit` messages.

Case outcomes:
- **A1** `y'+x y==Exp[3x]` → `Erf` closed form (was `$IterationLimit`).
- **A2** `y'+y==Q[x]` → `E^{-x}(∫E^x Q dx + C[1])` (was WRONG `Q[x]+C[1]E^{-x}`).
- **B5 / B-legendre2** → `LegendreP/Q[3/2 · , x]`, `[3/4, x]` (was series).
- **C4** Riccati → Legendre-based `u=w'/w` (was series).
- **B10** → `BesselI/K[0, (2/5)e^{5x/2}]` (was series).
- **D19** → `C3 − C2/(3(4x+7C1)^{3/4})` (was inert).
- **E16** → Airy products; **E17** → Bessel products (were inert/series).
- **B8** → `Hypergeometric2F1[a,b,c,x]` + 2nd soln, symbolic c (was hang).
- **B7, E18** → terminate cleanly (series/inert; closed form is Tier 2).

Files: new `dsolve_symsquare.c`, `dsolve_lower_reduce.c`; edited `dsolve.c`
(wiring), `dsolve_specialform.c` (Legendre + exp-Bessel + normal-form-Bessel +
Gauss symbolic-c), `deriv.c` (LegendreP D-rule), `dsolve_common.c` (integral
form + verify normalization), `dsolve_undetcoeff.c` (verify y_p), `dsolve_riccati.c`
(special-fn Simplify guard), `dsolve_kovacic.c` / `dsolve_operator_factor.c` /
`dsolve_factorable.c` / `dsolve_autonomous.c` (anti-hang guards), `zero_test.c`
(symbolic-exponent guard), `tests/CMakeLists.txt`, docs + changelog.

Memory: valgrind on the new cases shows the new code is leak-clean; the one
extra 1,832-byte block (E16) is the documented pre-existing FLINT
`rat_canon`/`Together` leak (`ratcanon.c:877`, via `FactorTerms`), not new code.

Tier 2 roadmap (not built): B9 (incomplete-Gamma 2nd soln), B12 (Sech²/assoc-
Legendre), B13 (Integrate `x^p trig(Log x)`), D14/D15 (autonomous implicit
quadrature / Abel-2nd-kind), B7/E18 closed forms, A6, B11.

## Progress log
- T1 done (A1 Erf, A2 integral, undetcoeff verify, verify normalization). Suites green.
- T2 done (Legendre B5/Bleg2/C4; LegendreP D-rule; Riccati special-fn Simplify guard). Green.
- T3 done (exp-arg Bessel B10). Green.
- T4 done (missing-lower-deriv reduction D19; new method wired). Green.
- T5 done (symmetric-square E16 Airy / E17 Bessel; normal-form Bessel branch;
  new method wired; CMake COMMON_SRC updated). Green.
- T6 done (anti-hang: Kovacic numeric-coeff gate, operfactor pole gate,
  factorable derivative-degree gate, autonomous missing-x pre-gate). B7 & E18
  terminate. Suites green.

## Follow-up — General Cauchy-Euler solver (post-Tier-1)

Reported hang `DSolve[x^2 y'' + y == x^2]` (complex-root inhomogeneous Euler)
plus two more requests: `x^2 y'' − 2x y' + 2y == x^2 e^x` (ExpIntegralEi) and
`(x+1)^2 y'' − 3(x+1) y' + 3y == x^2` (shifted centre). Root cause: the old
`dsolve_euler.c` used x-domain variation of parameters unconditionally, and its
final `ds_simplify(yp)` hung on the trig-of-Log products with irrational
frequencies that only complex roots produce.

Fix (`src/calculus/dsolve_euler.c`, full rewrite): a fully algorithmic solver —
- centre detection `b = x − n c_n/c_n'` (handles shifted `(x−b)^k`);
- **real roots** → `(x−b)^r` basis + x-domain VoP (keeps `ExpIntegralEi` etc.);
- **complex roots** → reduce to constant-coefficient via `(x−b)=e^t`, solve with
  the const-coeff engine, map back `t→Log[x−b]` (no trig-of-Log Simplify hang).

Tests: `t_euler_inhomogeneous_complex` + `t_euler_regression_corpus` (In[1]–In[11]
by residual). All DSolve suites green; `make check-c99` clean. Memory stable over
a 250× solve loop (`MemoryInUse` 11.0→11.5 MB, history/memo only).
