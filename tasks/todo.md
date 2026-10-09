# Conditional identities under Assumptions — general additions

Plan: `~/.claude/plans/let-s-test-the-implementation-bright-knuth.md`
Rule: fully general + sound only; NO one-off hacks (false zeros are the worst outcome).

## Phase 1 — Class 1 (clean, reuse existing provers)
- [ ] **B** inverse-of-forward family (`ArcTanh[Tanh[e]]->e` etc.) in
      `simp_assume_rewrite.c`, gated `prov_re`/`prov_nn`. → hyp #92.
- [ ] **D** whole-input ComplexExpand under reality in `simp_assume_rewrite.c`. → trig #83.
- [ ] **A** radical-product combine `Sqrt[a]Sqrt[b]->Sqrt[ab]` under provable
      nonneg (`denest_is_nonneg`), `simp_denest.c` + wire ctx in `simp_search.c`.
      → hyp #66, #68; enables #74.

## Phase 2 — Class 2 region engine (sound derivative-constancy)
- [ ] Shared budgeted Reduce-entailment helper (from `refine.c`).
- [ ] `simp_region_identity_is_zero(e,ctx)` new `simp_region_identity.c`:
      ∇f≡0 (exact) + real box + Reduce branch-cut-domain entailment + exact
      sample-point value via extended invtrig certifier. → #75/#76/#77 (trig+hyp).
- [ ] extend `simp_invtrig_combo_is_zero` to ArcTanh/ArcSinh/ArcCosh constants.
- [ ] #101 best-effort (trig square recognition + Reduce `Cos>=Sin` on 0<x<Pi/4).

## Cross-cutting
- [ ] Assumption harness over both corpora; **adversarial soundness tests**
      (no-assumption / wrong-region must NOT reduce; non-identities never certify).
- [ ] No regression (corpora 93/101 & 85/93 bare; suites; check-messages/c99; valgrind).
- [ ] version bump + changelog + docs; commit/tag.

## Already pass under assumptions (verified): trig #100, hyp #93.

## Review (DONE — v0.329)

**Outcome: 11 of 14 conditionals now reduce under Assumptions, via general sound
mechanisms (no hacks). Closed: #66,#68,#74 (radical combine), #92 (inverse-of-
forward), #83 (ComplexExpand), #75 trig+hyp, #76 hyp, #77 hyp (region engine +
hyperbolic-addition recognizer). + #100/#93 already. Residual: trig #76/#77
(ArcSin/ArcCos addition — radical range inequality beyond Reduce), trig #101.**

Built: Reduce-entailment bridge (`simp_assume.c` assume_reduce_entails/_nonneg);
radical-product combine (`simp_denest.c`); inverse-of-forward hyperbolic family
(`simp_assume_rewrite.c`); region derivative-constancy engine + hyperbolic const
& addition certifiers (`simp_builtins.c`). Reduce enters Simplify only on the
decline branch, budgeted (var cap).

Verified: both corpora bare unchanged (93/101, 85/93); 19 regression suites;
14-case adversarial soundness (no false zeros — no-assum/wrong-region/non-identity
all stay); valgrind clean; check-messages/check-c99. Tests:
`test_invtrig_simplify.c::{test_simplify_under_assumptions,
test_simplify_assumptions_soundness}`. Docs + changelog; version 0.329.
