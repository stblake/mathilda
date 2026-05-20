# Risch-Norman Implementation Tracker

Implementation per `plans/RISCH_NORMAN_PLAN.md`.  Reference:
`parallel_risch/pmint.maple`.

## Phase 1 — Plumbing & skeleton ✅

- [x] Create `src/intrischnorman.{c,h}` with stub `builtin_rischnorman`.
- [x] Wire `intrischnorman_init()` into `integrate_init()`.
- [x] Extend `src/integrate.c` dispatcher to fall through to RischNorman.
- [x] Create `tests/test_intrischnorman.c` with scaffolding tests.
- [x] Register test binary in `tests/CMakeLists.txt`.
- [x] Add `src/intrischnorman.c` to `COMMON_SRC`.

## Phase 2 — Trig→Tan, indet collection, substitution maps ✅

- [x] `convert_to_tan(f, x)` — full Weierstrass on Sin/Cos/Tan/Cot/Sec/Csc and hyperbolic siblings.
- [x] `decot_rec(f)` — undo evaluator's 1/Tan → Cot rewrite to keep field generators pure.
- [x] `convert_sincos_to_tan(f, x)` — partial Weierstrass for the post-hoc verifier (Sin/Cos/Sec/Csc only).
- [x] `pythagorean_rewrite(f)` — Sec[u]² → 1+Tan[u]², Csc[u]² → 1+Cot[u]², etc.
- [x] `collect_indets_closed(ff, x, &si, &n)` — atoms closed under one round of D, with Cot/Coth canonicalized to Tan/Tanh.
- [x] `build_substitution_maps(si, n, &lin, &lout, &vars)` — fresh pmint$v_i symbols.  `pm_sub_map_to_rule_list` emits companion reciprocal rules for Tan ↔ Cot pairs.
- [x] Testable surfaces: `Integrate`Helpers`PMConvertToTan`, `PMSincosToTan`, `PMPythagoreanRewrite`, `PMCollectIndets`, `PMSubstMap`.

## Phase 3 — Vector field, splitFactor, deflation, monomials ✅

- [x] `build_vector_field(li, n, lin, x, &l, &q)`.
- [x] `apply_d(f, vars, l, n)` — Σ l[k] · D[f, vars[k]].
- [x] `split_factor(p, vars, l, n, &s, &h)` — pmint.maple:80-90.
- [x] `deflation(p, vars, l, n)` — pmint.maple:92-98.
- [x] `enumerate_monomials(vars, nv, dg, &monoms, &n)` with `PMINT_MAX_MONOMIALS` cap.
- [x] Testable surfaces: `Integrate`Helpers`PMVectorField`, `PMApplyD`, `PMSplitFactor`, `PMDeflation`, `PMEnumerateMonoms`.

## Phase 4 — Candidate ansatz + linear system + RowReduce solve ✅

- [x] `build_candidate(monomials, nm, &cand_num, &A_names, &nA, &unknowns)`.
- [x] `walk_coefficient_table` — recursive coefficient extraction over multivariate.
- [x] `linear_builder_add_row` — extracts linear-form rows for each monomial.
- [x] `solve_linear_undet(equation_numer, vars, nvars, unknowns, nunk, &sub, &status)` — RowReduce + decode.
- [x] `try_integral_full` — polynomial-equation form (no Together inside hot path).
- [x] 8 Q-rational integrands verified: Exp[x], Exp[2x], x Exp[x], x² Exp[x], Log[x], x Log[x], Sin[x] Exp[x], Cos[x] Exp[x].

## Phase 5 — Log candidates + getSpecial + K=I retry ✅

- [x] `my_factors(p, over_Qi, vars, nv, &factors, &n)` — Factor[p] over Z or Factor[p, Extension -> I] via Trager.
- [x] `get_special_all(si, vars, n, &darboux, &flags, &count)` — Darboux polys for Tan (1+Tan²), Tanh (1±Tanh), LambertW.
- [x] Integral-flagged Darboux multiplied into cden's `s` part.
- [x] `try_integral_full` extended with `Σ B_j Log[g_j]` candidate and polynomial-equation form.
- [x] K=I retry loop in `rischnorman_integrate` (re-factors candlog polys over Q(I)).
- [x] 10 log-bearing integrands verified: 1/x, Sin[x], Cos[x], Tan[x], Cot[x], 1/(1+Exp[x]), Exp[x]/(1+Exp[x]), 1/(x Log[x]), Log[x]², (1+Log[x])/(x Log[x]).

## Phase 6 — Dispatcher polish + post-hoc verification ✅

- [x] Wall-clock budget via `setitimer` + SIGALRM + sigsetjmp (PMINT_WALL_CLOCK_SEC = 30s).
- [x] Polish docstrings; updated `docs/spec/changelog/2026-05.md`.
- [x] Extended test corpus to 50+ tests across all phases.
- [x] 10 extended-corpus successes: Exp[a x], x³ Exp[x], Log[x]³, Sinh[x], Cosh[x], Tanh[x], Sin[x]², Cos[x]², Tan[x]², 1/(1+x²).
- [x] 5 known-fail cases bubble back cleanly: 1/Log[x], Exp[x²], Log[x] Cos[x], Exp[x]/Log[x], Exp[Sin[x]].
- [x] Full regression: every other `tests/` binary still passes.

## Deferred follow-ups

- [ ] **Post-hoc verifier in production path** — currently relies on the algebra; a verifier that reduces trig-half-angle output against Sin/Cos integrands would catch any silent bugs.  Test corpus exercises this manually via explicit Weierstrass rule lists; production path skips it.
- [ ] **Multivariate Q(I) factoring** — `Factor[_, Extension -> I]` is univariate-only via Trager.  Unlocks integrands like `1/(1 + Sin[x]^2)` (requires multivariate Q(I) factoring of the candlog denominator).
- [ ] **BesselJ Darboux specials** — Bronstein 2004's worked example uses a Bessel-function differential field.  Requires `BesselJ` to be implemented in Mathilda first.
- [ ] **LambertW integration** — Darboux poly is wired (`get_special_all`), but no test corpus yet.
- [ ] **Trig-half-angle output cleanup** — pmint's results live in Tan[x/2] form, which is mathematically correct but less readable than the equivalent Sin/Cos form.  A post-process Tan-half-angle → Sin/Cos pass would improve user output.
- [ ] **Algebraic-function integrands** (`Sqrt[1+x^3]`, etc.) — out of scope for pmint's heuristic; would require a separate algebraic-integration module (Trager's algorithm).
