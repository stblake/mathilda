# Task: Reduce over Reals — general elementary-real-function sign diagram

## Goal
Make `Reduce[..., x, Reals]` (and the 2-arg default where the statement forces
Reals) solve univariate statements built from polynomials, `Abs`, real radicals,
rational poles, `Log`, inverse-trig, and isolated `Floor`/`Ceiling`/`Round`/`Mod`
— generally and soundly. Target: the 15 cases (0–14) in the approved plan
(`~/.claude/plans/wise-gathering-simon.md`).

## Design (one engine + one head table)
- **Engine A** (`reduce_realdiag.c`): univariate real sign diagram tolerating
  radical/pole atoms — breakpoints from `Solve` roots ∪ poles ∪ domain boundaries;
  truth via domain gate + pole gate + substitution-evaluate; numeric-sign fallback
  for transcendental (Pi) breakpoints. Fallback after poly `reduce_univar` declines.
- **Preprocess B + domain table C** (`reduce_realfn.c`): `Abs` sign-split,
  `Mod`→`Floor` + integer-part isolation, and the head→(domain constraint, boundary)
  table.
- **Dispatch** (`reduce.c`): detect real-fn of x → preprocess + route to Reals →
  general-engine fallback.

## Steps
- [x] 1. Extract emission `rru_emit_sign_diagram` (shared) from `reduce_univar.c`.
- [x] 2. `reduce_realfn.{c,h}`: head-domain table (`reduce_real_domain_collect`),
        `reduce_stmt_has_realfn`.
- [x] 3. `reduce_realdiag.{c,h}`: `reduce_univar_general` (breakpoints, domain/pole
        gate truth oracle, numeric-sign fallback).
- [x] 4. `reduce_realfn.c`: `reduce_realfn_preprocess` (Abs split, Mod/Floor).
- [x] 5. Wire `reduce.c` dispatch (detector + preprocess + force reals + fallback).
- [x] 6. Build; run cases 0–14; 14/15 exact (case 14 x>0 vs x>=0: sound 0/0 point).
- [x] 7. Tests: `test_reduce.c` decline case updated + `test_real_functions` added;
        corpus records + verifier green (118/118); decline soundness checks.
- [x] 8. valgrind (new code leak-clean) + `make check-c99` (pass).
- [ ] 9. Docs: `REDUCE_PLAN.md` Phase 9, `docs/spec/builtins/`, changelog 2026-08-25.

## Review

Implemented Phase 9: a general univariate real sign-diagram engine for elementary
real functions. Two new modules (`reduce_realfn.{c,h}` preprocessing + domain
table, `reduce_realdiag.{c,h}` engine), a shared emission helper in
`reduce_real_util.c`, a branch-cut-transcendental `Together` guard in
`reduce_atom.c`, and dispatch wiring in `reduce.c`.

14/15 target cases match Mathematica exactly; case 14 returns the sound `x>0`
where MMA lists `x>=0` (at x=0 the expression is literally `0/0+0/0`). Also fixed
pre-existing gaps: `Reduce[Sqrt[x-1]==2]`→`x==5`, `Reduce[Abs[x]<1,Reals]`→`-1<x<1`.

Bugs found & fixed during impl: (1) `expr_new_function` arity-3 on a 2-element
Plus array (garbage read → wrong Mod result + segfault); (2) `Solve` returning the
identity `{{}}` hid isolated polynomial-factor roots → added `collect_factor_roots`;
(3) `Together` applied complex-log identities unsound over ℝ (`Log[x^2]-2Log[-x]`
→ `-2 I Pi`) → skip Together for branch-cut transcendentals.
