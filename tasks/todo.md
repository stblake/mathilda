# DSolve M6 — PDE Charpit + factorable 2nd-order (lower-order terms)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-purring-cherny.md`

## Tasks

- [x] 1. Extended `dsolve_pde2.c`: symbol factoring with lower-order terms →
      `Σ e^{-m_i v1} C[i][ξ_i]` (distinct/repeated λ; A=0 swap; pure-mixed).
      Backward compatible (m=0 = current output; wave/Laplace regressions pass).
- [x] 2. Substrate: `PDERelation[Ψ]` path in `dsolve_common.{c,h}`; factored the
      shared implicit-diff core `pde_implicit_residual_ok`.
- [x] 3. `src/calculus/dsolve_pdecharpit.c` — PDECharpit standard forms I/II/III.
- [x] 4. Wired into `dsolve.c` (extern decls, cascade after Clairaut, init).
- [x] 5. Wired `dsolve_pdecharpit.c` into `tests/CMakeLists.txt` COMMON_SRC.
- [x] 6. Unit tests `t_pde2_lower_order`, `t_pde_charpit`.
- [x] 7. Stress `t_stress_pde2_factorable` (5), `t_stress_pde_charpit` (7).
- [x] 8. 3 dsolve ctest suites green (184/37/9); `make check-c99` PASS.
- [x] 9. valgrind: no new-code frames; delta is pre-existing Integrate/Solve leak.
- [x] 10. Docs: calculus.md (Charpit row + pde2 row), changelog, DSOLVE_PLAN.md.

## Review

Two closed-form deferred M6 items landed, each back-substitution verified.

**PDECharpit** (`dsolve_pdecharpit.c`) — first-order fully nonlinear PDE
`F(v1,v2,u,p,q)==0` by Charpit's method, three standard forms: F(p,q) →
`u=C[1]v1+q v2+C[2]`; F(u,p,q) → implicit `∫du/P == v1+C[1]v2+C[2]`; separable
`f(v1,p)==g(v2,q)` → `u=∫P dv1+∫Q dv2+C[2]`. Explicit forms reuse the PDEBranches
verify; the implicit F(u,p,q) form uses the new PDERelation substrate path
(implicit-function-rule verify with arbitrary constants).

**PDELinearSecondOrder + lower-order terms** (in place) — factors the full symbol
`A ξ²+B ξη+C η²+D ξ+E η+F` into two first-order operators → exponential-damped
`Σ e^{−m_i v1} C[i][v2+λ_i v1]`. Distortionless telegraph, damped/convection,
pure-mixed. Backward compatible (D=E=F=0 → principal-part form). Declines the
non-factorable general telegraph (needs Bessel).

Worked examples verified end-to-end:
- `u_tt − c² u_xx + a u_t + (a²/4)u == 0` → `e^{−a t/2}(C[1][x+ct]+C[2][x−ct])`
- `u_xy + u_x == 0` → `C[1][y] + e^{−y}C[2][x]`
- `p²+q²==1` (Charpit I), `p q==u` (II, implicit), `p²−q²==x−y` (III)

Deferred (documented future): inhomogeneous 2nd-order forcing; non-factorable
general telegraph (Bessel); Charpit general integrable-combination + singular
solutions; genuinely-quasilinear non-conservation.
