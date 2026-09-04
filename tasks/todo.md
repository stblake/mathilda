# DSolve M10 Lie — hang fix + `abaco2_unique_general` (§4.4.2) + `chi` (CPC 101)

Continuation of the M10 Lie point-symmetry work. Three items (user request):

## Phase A — fix the pre-existing hang (undefined-function omega) ✅ DONE
Fixed & verified: pinned Lie A 0.26s / B 1.1s (was 10s / 82s+hang); all 3 DSolve
suites green + new t_lie_undefined_function_declines; check-c99 clean; valgrind
leak-flat (identical totals to a rational-ODE decline baseline). Fix = node-budget
guards on the hot-path helpers (lie_ratsimp/free_of_var/is_zero/sep_xfactor), a
structural `lie_lit_zero` polynomial zero-test replacing the hanging general zero_test,
lie_check fast poly pre-check + full test gated on `lie_has_undefined_function`,
lie_first_integral declines undefined-function integrands, and the classifier ansatze
(abaco1_product/function_sum/abaco2_similar) skipped when omega has an undefined fn.

## (original Phase A notes)
Root cause found: `Tan[ArcTan[y]+F[x²+y²]]` (undefined `F`) makes the fast helpers
`lie_free_of_var` / `lie_is_zero` do a polynomial `Expand` that blows up
(pre-Expand num LeafCount 217K → post-Expand 1M+) inside `function_sum` (and would
recur in `abaco1_product` / `abaco2_similar`). The heuristic never reaches
`abaco2_unique_unknown`, which would otherwise handle the paper's Eq-70 class.
- [ ] Add a capped node-counter `lie_too_big(e, budget)` (O(budget), short-circuits).
- [ ] Guard `lie_free_of_var` (bail→"not free of") and `lie_is_zero` (bail→"not zero")
      before the `Expand`. Bailing is safe: Lie is heuristic; a declined branch never
      gives a wrong answer and the correct heuristic still runs; back-sub verifies all.
- [ ] Remove the temporary MATHILDA_LIE_DEBUG probes.
- [ ] Verify: user's ODE + paper Eq-70 no longer hang; paper Eq-70 finds `[y,-x]` via
      `abaco2_unique_unknown` (declines at the non-elementary quadrature — correct, no
      inert head); existing Lie tests unchanged.

## Phase B — `abaco2_unique_general` (§4.4.2, Cheb-Terrab & Roche 1998)
Two parts:
- [ ] **B1 — §4.4.1 general ("differential invariant of order zero"), Eqs 73–81.**
      Extend `lie_abaco2_unique_unknown`: for each kernel M with R=M_y/M_x, also try the
      candidates `[-R, 1]` (Eq 75, pattern [f(x)g(y),1]) and `[1, -1/R]` (Eq 78, pattern
      [1, f(x)+y h_x/h]) — no separability required. Catches family (77) / Kamke 433.
- [ ] **B2 — §4.4.2 Case I / Case II, Eqs 82–90.** New `lie_abaco2_unique_general_cand`
      (reuse `lie_run_with_inverse` for the [G(y),F(x)] inverse pattern). φ=Log[ω];
      A=φ_xy, B=φ_yy+φ_y², C=φ_xx−φ_x²; D (Eq 85). D==0 → Case I (E2==0,E3==0,E1≠0;
      η=Exp[∫…dy] free of x, ξ=−4A³η/E1). D≠0 → Case II (Eq88==0,E5==0,E6==0,E4≠0;
      η=Exp[∫…dy], ξ=−E4η/D). Gate each by `lie_check` + `lie_first_integral`.
- [ ] Wire after `abaco2_unique_unknown`, before `bivariate`.
- [ ] Tests: forward-generate from invariant family (62) y'=(f_x/g_y)(K(f+g)−1).

## Phase C — `chi` (CPC 101 1997, 5th algorithm), Eqs 9–10
η=ξΦ+χ; since S(ξ,ξΦ+χ)=S(0,χ), χ solves the linear PDE χ_x+Φχ_y−Φ_yχ=0 (Eq 10) for
ANY ξ → take symmetry [0, χ]. Ansatz: χ = Σ_k p_k(x,y)·basis_k, basis_k = building
blocks from ω (functions / non-integer powers / their reciprocals) + 1; p_k degree-d
polys in x,y with undetermined coeffs. Substitute into Eq 10, treat basis functions as
independent generators, CoefficientList → NullSpace → χ. Gate + integrate.
- [ ] First cut: polynomial + single-function-block basis; document richer-basis future.

## Gates (each phase)
- [ ] `dsolve_tests`, `dsolve_stress_tests`, `dsolve_m5_stress_tests` green
- [ ] `make check-c99` clean; per-solve valgrind leak-flat (modulo known FLINT rat_canon)
- [ ] docstring, `docs/spec/builtins/calculus.md`, `docs/spec/changelog/2026-08-31.md`,
      `docs/design/dsolve_lie_symmetry.md`, `DSOLVE_PLAN.md` updated

## Review — DONE (2026-09-04)

**Delivered (both verified):**
1. **Hang fix.** Root cause: on an `ω` with an undefined function of both variables the
   quadrature classifiers balloon (217k→>1M nodes under `Expand`) and the general
   `zero_test`/`Integrate` hang on the arbitrary-function atoms. Fix = node budget
   (`LIE_EXPR_BUDGET=6000`) on `lie_ratsimp`/`lie_free_of_var`/`lie_is_zero`/
   `lie_sep_xfactor`; structural `lie_lit_zero` polynomial zero-test replacing the
   hanging `ds_is_zero` in the fast helpers; `lie_check` = fast poly test + full test
   only when `!lie_has_undefined_function`; `lie_first_integral` declines
   undefined-function integrands; `abaco1_product`/`function_sum`/`abaco2_similar`
   skipped on undefined-function `ω`. Pinned Lie 0.26s/1.1s (was 10s/hang); valgrind
   leak-flat; `t_lie_undefined_function_declines`.
2. **Order-zero extension of `abaco2_unique_unknown`** (§4.4.1 Eqs 73–81): non-separable
   candidates `[-R,1]`/`[1,-R]`/`[1,-1/R]`. Catches Kamke 433 → `x−Sqrt[x²+xy+a]==C[1]`.
   `t_lie_abaco2_order_zero`.

3. **`chi` (CPC 101 5th algorithm)** — the `η=ξω+χ` reformulation with a rich
   transcendental-atom basis for `χ`. Solves Kamke 357 → `−x+Log[x]Sec[y[x]]==C[1]`
   (~1.2s, verified), the first genuinely-transcendental `χ` beyond `bivariate`. Trig `ω`
   now skips the rational/algebraic classifiers (→ `chi`); `abaco1_simple` uses bounded
   `lie_ratsimp`. `t_lie_chi`; valgrind leak-flat. See [[project_lie_chi_rich_basis]].
4. **All DSolve methods REPL-callable + docstrings.** Added pinned `DSolve`<Name>`
   builtins for `DecoupleSystem`/`TriangularSystem`/`LinearFirstOrderSystem`/
   `PDELinearFirstOrder` (`dsolve_method_builtin_system`/`_pde`); rewrote the `DSolve`
   docstring to catalog all methods. `t_sys_pde_pinned_methods`.

**Deliberately NOT done (findings, discussed with rationale):**
- **`abaco2_unique_general` (§4.4.2 Case I/II)** — documented exemption. The paper calls
  the closed-form route "very inefficient, if not just impractical"; every
  `[F(x),G(y)]`-symmetric ODE I could construct from the invariant family is already
  solved by Riccati/Bernoulli/kernel methods, so it cannot be tested non-vacuously, and
  its "huge expressions" would trip the robustness node budget anyway.
  (was: `chi` — now DONE, see item 3 above.)

**Gates:** all 3 DSolve ctest suites green; `make check-c99` clean; valgrind leak-flat.
