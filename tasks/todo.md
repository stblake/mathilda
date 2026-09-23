# TODO: DSolve M59 — Inactive primitive + implicit-first-integral autonomous companion

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-binary-hartmanis.md`
Unblock M58's deferred lever: return `∫dy/p == x + C[1]` (inert Integrate) for non-elementary
stage-2 quadratures. Key finding: `Inactive[Integrate][f,x]` ALREADY stays inert (no spin);
only D-FTC + registration are missing.

## Part 1 — Inactive / Activate primitive
- [ ] `SYM_Inactive`, `SYM_Activate` (sym_names.h/.c); register `Inactive`/`Activate` (Protected) in core.c.
- [ ] `Activate[expr]` = `expr //. Inactive[h_]->h` then re-eval (builtin or .m rule).
- [ ] deriv.c: FTC for `Inactive[Integrate]` — var matches → f; var differs → D-under-integral
      (df==0 → RETURN 0, else Inactive[Integrate][df,u]); definite Leibniz mirror.
- [ ] REPL: `D[Inactive[Integrate][1/p[y],y],y]==1/p[y]`; `D[...-x,x]==-1`; Activate round-trips.

## Part 2 — autonomous implicit companion
- [ ] Factor `ar_reduce(P,&n,&pfun)->pbody` out of dsolve_autonomous_try (explicit byte-identical).
- [ ] `dsolve_autonomous_implicit_try`: build G=Inactive[Integrate][1/pbody,Ysym]/.Ysym->y[x] - x;
      NUMERIC self-verify (reconstruct y'..y^(n) from chain, check original ODE ~0); return bare G.
- [ ] Wire dsolve.c: extern + 2nd cascade slot (dsolve_run_implicit) + pinned explicit-then-implicit.
- [ ] `make -j` + `make check-c99` green.

## Verify + release
- [ ] REPL: 263/264/267/268/1167/1168 → {{Inactive[Integrate][1/p,y[x]]-x==C[1]}}, back-sub ~0.
- [ ] Order-2 + 1143 explicit UNCHANGED (ar_reduce factoring behavior-preserving).
- [ ] Units + stress; all DSolve stress suites green.
- [ ] CLEAN corpus run (no concurrent Mathilda): §2.1.2 delta, 0 FAIL, lower baseline 628.
- [ ] Docs (Inactive/Activate spec, AutonomousReduction, changelog, DSOLVE_PLAN M59, STATUS).
- [ ] version.h 0.174→0.175; commit `; v0.175`; tag v0.175.

## Review
M59 complete (v0.175). Realised M58's deferred "big lever" — and it delivered **+9** deterministic
(the whole order-3 autonomous residue 263/264/267/268/1167/1168 + dups 710/711/712), 0 FAIL.

Two parts. (1) `Inactive`/`Activate` primitives: the decisive finding was that `Inactive[Integrate][f,x]`
is ALREADY inert (an unevaluated compound-head fixed point, like `Derivative[n][f][x]`), so no
evaluator surgery — just `ATTR_PROTECTED` registration + a D-FTC rule (`D[Inactive[Integrate][f,u],u]==f`,
deriv.c) + `Activate`. (2) `dsolve_autonomous_implicit_try`: where the M58 explicit method declines a
non-elementary stage-2 quadrature, return `Inactive[Integrate][1/p,y[x]]-x==C[1]` — built and verified
(via FTC) WITHOUT running Integrate, so `y'==Sqrt[y Log y+…]`'s 45 s spin → 0.0 s. A numeric self-verify
(reconstruct y'..y⁽ⁿ⁾ from the reduction chain) is the real correctness gate (dsolve_run_implicit's
verify is vacuous at order n≥2). M58 stage-1 factored into shared `ar_reduce`; explicit path unchanged.

§2.1.2 clean re-run **591 → 595, 0 FAIL** (deterministic +9; `3rd_high_reducible` 8→17 PASS). The 7 P→U
losses are all x-dependent `_with_linear_symmetries` timing-cluster cases (M59's autonomous companion
never touches them — verified 439/877/1153 are x-dependent); the 1 crash (1153) is the known
intermittent fork SIGSEGV (clean 3/3 in isolation). Gate baseline 628 → 619. inactive_tests +
m58/m59 stress + core/deriv/integrate regression suites green; check-c99 green.

Watch item (non-blocking): 1168's numeric self-verify takes 6.2s (messy `p=Sqrt[(C+cy³)/y]`) — under
the 8 s budget (passed), improvable by substituting numeric constants before building the chain.

Future lever: upgrade separable/fos/chini/exact implicit paths to emit inert first integrals for
their own non-elementary quadratures (they currently decline), now that Inactive exists.
