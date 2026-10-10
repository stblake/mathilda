# FullSimplify log→inverse-function collapse class (`logexp.m`)

Plan: add the full class of "Log of a trig/hyperbolic cofunction combination
collapses to a single inverse/linear form" identities to FullSimplify's rule tier.
Broadest scope: dual class (8) + half-angle (4) + coefficient-tolerant variants.

## Steps

- [ ] Build Mathilda; probe prerequisites:
  - [ ] confirm `ArcCoth` is a registered head (else use `ArcTanh[Sech[x]]`)
  - [ ] probe pipeline-form risk: what does `Simplify[Log[Sec[x]+Tan[x]]]` etc. produce?
- [ ] Numerically self-verify all 12 identity signs (`Chop[N[lhs-rhs]]` at several points)
- [ ] Create `src/internal/simp/transforms/logexp.m` (12 rules + coefficient variants)
- [ ] Wire manifest: `ModuleFile[Log]` in `src/internal/simp/FullSimplify.m`
- [ ] Rebuild; end-to-end check every rule in every configuration (per-rule matrix)
- [ ] Add corpus rows to `tests/fullsimplify_corpus.m` (every configuration)
- [ ] Add negative/"unchanged" assertions to `tests/test_fullsimplify.c`
- [ ] Build & run `test_fullsimplify` + `test_fullsimplify_corpus` (+ regression)
- [ ] Docs: `docs/spec/builtins/simplification.md` family list
- [ ] Changelog: `docs/spec/changelog/2026-10-05.md`
- [ ] Version bump `src/version.h` 0.339 → 0.340
- [ ] Rebuild code-review-graph
- [ ] Commit `; v0.340` + tag `v0.340`

## Review

Done. New library `src/internal/simp/transforms/logexp.m` registered under `Log`
via one `ModuleFile[Log]` manifest line in `FullSimplify.m`. Covers the full
Gudermannian log→inverse-function class: Sec/Tan, Csc/Cot, Cosh/Sinh, Coth/Csch
(both signs) + four half-angle forms, each two-term rule coefficient-tolerant in
both factored and distributed surface forms with a positivity guard.

Key findings during implementation:
- Base Simplify does NOT rewrite `Sec/Tan` to fractions (patterns match directly),
  and already collapses `Log[Cosh±Sinh]` to `±x` on its own.
- The pipeline factors positive sums (`2Sec+2Tan -> 2(Sec+Tan)`) but leaves
  differences distributed — so each pair needs BOTH a factored and a two-coefficient
  distributed rule. A single optional coefficient can't express the non-unit minus.
- All 12 identity signs verified numerically (`Chop[N[lhs-rhs]] == 0`).
- No integrate regression: integrators using FullSimplify produce the
  `1/2(Log[2+2Sin]-Log[2-2Sin])` surface form, which these patterns don't match;
  `Integrate[Sec[x],x]` byte-identical before/after.

Verification: corpus 39/39 (30 new configs); `fullsimplify_tests`,
`simp_tests`, `simplify_hang_tests`, `trig_tests`, `integrate_ramanujan_tests`,
`integrate_diffunderint_tests` all green. Docs + changelog updated; version 0.340.
