# DSolve — ExactODE (§1c higher-order exact linear equations)

Implement `DSolve`ExactODE` for linear ODEs of order ≥ 2 whose left side is a
total derivative `L[y] == d/dx(M[y])`: integrate once to the first integral
`M[y] == ∫g + C[n]` and recurse into the scalar cascade. Mirrors
`dsolve_reduce_order.c` (recurse + reuse `extract_applied`); constant `C[n]` is
contiguous with the sub-solve's `C[1..n-1]`, so no renumbering.

## Implementation
- [x] src/calculus/dsolve_exactode.c — three-function contract; b[j] recurrence,
      exactness test `a0 == b0'`, forcing ∫g guard, reduced-eqn recursion
- [x] src/calculus/dsolve.c — extern decls, DS_EXACTODE enum, "ExactODE" map,
      cascade slot after euler / before specialform, pinned case, init chain
- [x] tests/CMakeLists.txt — add dsolve_exactode.c to dsolve source block

## Tests
- [x] tests/test_dsolve.c — t_method_exactode, t_exactode_more (2nd + inhomog +
      3rd-order), t_exactode_declines (non-exact + first-order), auto-dispatch
- [x] tests/test_dsolve_stress.c — exact_ode_ok forward generator (total-deriv
      construction over 6 (b1,b0) pairs + inhomogeneous)

## Docs
- [x] docs/spec/builtins/calculus.md — DSolve`ExactODE row after EulerCauchy
- [x] docs/spec/changelog/2026-08-31.md — ## DSolve — ExactODE entry
- [x] DSOLVE_PLAN.md — §1c [✓] ExactODE; update M5 future-work note

## Gates
- [x] make -j build clean
- [x] REPL spot-checks (2nd/3rd-order, inhomogeneous, pinned, non-exact decline,
      ?docstring, Airy/Euler regressions) — all pass
- [x] ctest -R dsolve (all 3 targets green, 34.9s)
- [x] make check-c99 (exit 0)
- [x] valgrind leak-check (1× and same-eqn 6× both 13,440 + 6,312 B — no per-call leak)
- [ ] commit + push

## Review — DONE

`DSolve`ExactODE` (§1c) solves linear ODEs of order ≥ 2 whose left side is a total
derivative `L[y] == d/dx(M[y])`: the first-integral coefficients come from the
recurrence `b_{n-1}=a_n`, `b_{k-1}=a_k−b_k'`, the exactness test is the leftover
`a_0 == b_0'`, and it integrates once to `M[y] == ∫g + C[n]` then recurses into
the scalar cascade on the order-(n−1) equation. New file
`src/calculus/dsolve_exactode.c` (mirrors `dsolve_reduce_order.c`: recurse into
`DSolve`, reuse `extract_applied`); constant `C[n]` is contiguous with the
sub-solve's `C[1..n-1]` (no renumbering); iterated exactness is free via the
recursion (3rd-order doubly-exact reduces twice).

**Verified:** `x y''+y'==0 → C[1]+C[2]Log[x]`; inhomogeneous `x y''+y'==x`;
3rd-order `x y'''+y''==0` (3 constants); pinned `DSolve`ExactODE`; non-exact
(Airy) + first-order both decline; Airy/Euler-Cauchy regressions intact.
3/3 dsolve ctest green; check-c99 exit 0; no per-call valgrind leak.

**Placement:** after EulerCauchy, before SpecialFunctionForm (elementary
order-reduction preferred over special-function/Kovacic/series). Reference
Airy/Bessel/Euler/Kovacic equations are all non-exact, so their pinned tests are
untouched.

**Follow-ups (out of scope):** OperatorFactor/DFactor (differential-operator
factoring); integrating-factor (adjoint) exactness `L*[μ]==0`; nonlinear
total-derivative detection; Sturm–Liouville EigenvalueProblem (§1f).
