# DSolve M10 Lie — `function_sum` (L2) + `abaco2_unique_unknown` (L3)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-humble-magpie.md`

## function_sum (§4.2) — completes L2
- [ ] `lie_function_sum_cand` extractor: `R=D[1/ω,{x,2}]`; decline if R≡0; x-factor of
      `D[1/R,y]` via `lie_sep_xfactor`; `ξ=1/(xfac·R)`, `η=0`
- [ ] `lie_function_sum` = `lie_run_with_inverse(...)` (inverse pattern free)
- [ ] wire into `dsolve_lie_try` after `abaco1_product`
- [ ] REPL: find isolating ODE(s), verify back-sub == 0, confirm earlier heuristics decline
- [ ] unit test `t_lie_function_sum` + register
- [ ] stress `t_stress_lie_function_sum` + register

## abaco2_unique_unknown (§4.4.1) — opens L3 abaco2_unique
- [ ] `lie_collect_kernels` tree-walk (non-int powers + non-arith fn-apps of both vars, dedup)
- [ ] `lie_abaco2_unique_unknown`: per-kernel `R=M_y/M_x`; sep x-factor X;
      candidates `[X,-X/R]` and `[-R/X,1/X]`, gated + integrated
- [ ] wire into `dsolve_lie_try` after `abaco2_similar`
- [ ] REPL: isolating ODE (arbitrary-fn / non-integer-power kernel), verify, isolation
- [ ] unit test `t_lie_abaco2_unique_unknown` + register
- [ ] stress `t_stress_lie_unique_unknown` + register

## Docs
- [ ] `docs/spec/builtins/calculus.md` heuristic list
- [ ] `docs/spec/changelog/2026-08-31.md` two sections
- [ ] `docs/design/dsolve_lie_symmetry.md` staging ✅ + §4.2/§4.4.1 notes
- [ ] `DSOLVE_PLAN.md` tick both
- [ ] docstring + file-header staging comment in `dsolve_lie.c`

## Gates
- [ ] REPL spot-checks green
- [ ] `dsolve_tests` + `dsolve_stress_tests` green, wall-clock clear of watchdog
- [ ] `make check-c99` clean
- [ ] valgrind per-solve leak-flat (own code; pre-existing FLINT leak out of scope)
- [ ] rebuild `./Mathilda` + code-review-graph

## Review — DONE (2026-09-04)

Both heuristics landed in `src/calculus/dsolve_lie.c` (+152 lines), chained cheapest-first:
`abaco1_simple → linear → abaco1_product → function_sum → abaco2_similar →
abaco2_unique_unknown → bivariate`.

**Key correction:** §4.2's classifying quantity is `ω·∂²ₓ(1/ω) = F''/(F+G)` (rational —
the leading `ω` cancels the transcendental part of `1/ω`), NOT `∂²ₓ(1/ω)` alone (my
first read of the OCR'd Eq 27). Same as SymPy's `odefac*(1/odefac).diff(x,2)`. With that
fix `function_sum` solves its family cleanly and fast.

**Findings:**
- `function_sum` ODEs are intrinsically transcendental (Log/ArcTan) — the paper (§3)
  even notes these patterns don't help Kamke ODEs. Tests use invariant-family members
  built with F=1/x,G=y. The F=1/x² ArcTan member verifies but is ~7 s, so omitted.
- `abaco2_unique_unknown` value is the non-integer-power cases (integrate elementarily,
  e.g. `(x/y)(x²+y²)^(1/3)` → `[1/x,−1/y]`). Genuinely-arbitrary-function ODEs find the
  symmetry but the `∫F/(1+F)` quadrature is non-elementary → decline (correct).
- **Pre-existing hang** (out of scope, NOT introduced): `Tan[ArcTan[y]+F[x²+y²]]`
  (undefined `F`) hangs in `abaco2_similar`/earlier — confirmed by disabling both new
  heuristics; still hangs. Avoided undefined-function ODEs in tests.

**Gates (all green):** isolation confirmed by per-heuristic attribution (both new
heuristics are the actual solvers, non-vacuous); `dsolve_tests` 30.6 s,
`dsolve_stress_tests` 47 s, `dsolve_m5_stress_tests` 8 s; `make check-c99` clean;
per-solve valgrind — `function_sum` leak-flat, `abaco2_unique_unknown` adds only the
pre-existing 56 B/call FLINT `rat_canon` leak (`tasks/flint_ratcanon_leak.md`). Auto
cascade (unpinned `DSolve[]`) solves both new classes — real coverage, not just pinned.
Docs: changelog, calculus.md, design doc §4.2/§4.4.1, DSOLVE_PLAN.md all updated.
