# DiffUnderInt: finite-domain Feynman families (3 types)

## Goal
Make `Integrate[..., Method -> "DiffUnderInt"]` solve three parametric families:
- **F1 power-log**: `Log[1+c x^p]/(x Sqrt[1-x^(2p)])` on {0,1} → `(Pi^2/8 - ArcCos[c]^2/2)/p`
- **F2 secant-radical**: `Sec[2x] Log[1+c Sqrt[1-Tan[x]^2]]` on {0,Pi/4} → `Pi^2/8 - ArcCos[c]^2/2`
- **F3 tangent-power**: `Csc[2x]^2 Log[1+Tan[x]^a]` on {0,Pi/4} → `(Pi Csc[Pi/a]-a)/4`

## Tasks
- [ ] Read the target regions of integrate_diffunderint.c precisely
- [ ] `strip_re_im` helper + wire into `relation_apply` (Re[a]>0 parsing)
- [ ] `insert_feynman_param(canon, p)` helper (single Log[1+W]→Log[1+p W])
- [ ] `normalize_power_sub` (u=x^p) — F1 front-end
- [ ] `normalize_tan_half` (rule-based t=Tan[x]) — F2/F3 front-end
- [ ] `inner_arccos_family` (Form A + Form B → A ArcCos[q]/Sqrt[1-q^2])
- [ ] `stage_finite_feynman` (F1+F2): normalize→insert→D→inner→closed-form G→verify→eval@1
- [ ] `stage_tangent_power` (F3): normalize→emit digamma→subst-identity+a0=3 verify
- [ ] Control-flow edit in `integrate_diffunderint_try` (drop np==0 bail, add stages)
- [ ] Header overview comment update
- [ ] Tests: 3 targets + corpus subset + emit-guard cross-checks
- [ ] Docs: calculus.md, tutorial, changelog 2026-08-31.md
- [x] Build, run tests, REPL spot-check, no-hang check, valgrind

## Review (done)
- All 3 user integrals return correct closed forms (Method -> "DiffUnderInt"):
  I1=π²/12, I2=π²/8, I3=(π Csc[π/a]−a)/4 (incl. Re[a]>0 variant).
- Held-out corpus: 21/21 cases across F1 (p,c), F2 (c), F3 (a) all Simplify to 0.
- New tests: test_finite_power_log / _secant_radical / _tangent_power /
  _arccos_emit_guard — all pass; full DiffUnderInt suite + regression suites
  (residue, newton_leibniz, dispatch, ramanujan, symmetry, beta, principalvalue)
  all PASS.
- Clean build (no warnings), make check-c99 clean, valgrind: 0 lost blocks trace
  to new code (193 baseline blocks are libobjc/flint/eval noise).
- Key implementation notes: parameter introduction (insert_feynman_param),
  emitted ArcCos inner (engine unevaluated/wrong), digamma reflection for F3 with
  a0=3 anchor, has_noninteger_var_power guard + tangent_power-first ordering to
  avoid the t^a × radical Simplify hang, PowerExpand for the p=1/2 (u²)^(1/2) case,
  Re[]/Im[] stripping in relation_apply.
- Known limitation (pre-existing, documented): F3 needs explicit Method; under
  Automatic an earlier method hits the hang-prone *indefinite* integral first.
