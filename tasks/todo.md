# DSolve M22 — §2.2.2 corpus + coverage waves

## Wave 0 — Corpus infrastructure
- [ ] Fix `tools/latex_ode_to_mathilda.py :: is_condition_row` (anchor to LHS)
- [ ] Download §2.2.2 HTML into DSolve_test_status (or keep source note) & regenerate `DE_examples_222.m`
- [ ] Regenerate `DE_examples_221.m` (must be byte-identical → proves no regression)
- [ ] Register `dsolve_corpus_2_2_2_tests` in `tests/CMakeLists.txt`
- [ ] Build harness + measure true 20s-forked baseline → `reports/2.2.2.{tsv,md}`
- [ ] Add §2.2.2 block to STATUS.md + README contents row

## Wave 1 — Homogeneous correctness (latent wrong-answers; do first)
- [x] Fix 1A: reduced function via `F(1,v)` not `F(x,v·x)` (dsolve_homogeneous.c) — fixes 117
- [x] Fix 1B REVERTED (too broad, broke 121/144); replaced with in-method numeric verify `homog_num_wrong` — fixes 112
- [x] Log-gate homog_exp_log_invert + reject `$rad` placeholder leaks — fixes 118
- [x] Verify 112/117/118 (117 explicit-verified; 112/118 implicit) + 121/144 restored

## Wave 2 — Exact transcendental
- [x] Refactor potential into `exact_potential()`; add `dsolve_exact_implicit_try`
- [x] Wire implicit fallback after explicit Exact (dsolve.c :343 + pinned)
- [x] Linearizable Bernoulli-shape recursion gate (fixes 141 pre-Exact hang)
- [x] Verify 140/141/182/195 solved

## Wave 3 — FirstOrderSubstitution implicit fallback
- [x] Add `dsolve_fos_implicit_try`; wire after explicit FOS (dsolve.c :386 + pinned)
- [x] Verify 159 (inert-integral implicit, as Mathematica); 165 = slow-explicit residue
- [x] Bernoulli mixed-radical gate → 107 solves via Homogeneous

## Regressions found & fixed (from ctest/corpus)
- [x] Fix 1B (shared extract) reverted → 121/144 + §2.2.1 FAIL restored
- [x] `ds_contains` pointer bug in Log-gate → log-family `(x+2y)/(2x+y)` explicit restored (dsolve_stress)
- [x] Linearizable gate moved to recursion-only branch → `y'=y(e^x+Log y)` restored (dsolve_tests)

## Residue + docs
- [x] Residue documented (133/160/165/170/175-178) in STATUS.md
- [x] DSOLVE_PLAN.md M22 entry + docs/spec/changelog + reports/2.2.2.{tsv,md}
- [x] All 8 dsolve unit+stress ctests green; check-c99 green
- [ ] Final corpus re-measure (§2.2.2/§2.2.1/§2.1.2 on final binary) — running

## Review
**Result: §2.2.2 82→92/100 (0 FAIL, +10). §2.2.1 held at 96 (0 FAIL). Converter fix
benefits all sections (byte-identical §2.2.1 regen). 3 solver waves + Bernoulli gate,
all via verified implicit substrate. 3 self-inflicted regressions caught by ctest/corpus
and fixed. Delivery: uncommitted (user reviews).**
