# Mellin Integrate — assumptions correctness + general extensions

Plan: `~/.claude/plans/in-the-following-examples-serialized-wave.md`

## Part 0 — Assumptions correctness  [DONE, verified]
- [x] 0a. `def_has_undefined_function` → x-aware; In[6–8] now exact Gamma-ratios; protective case holds
- [x] 0b. Refine-augmented conjunct-wise strip discharge; In1/2/3/20 collapse to bare value

## Part 1 — Route A reductions
- [x] SinIntegral → 1F2 rule (In[14]) verified (check→0)
- [ ] StruveH → 1F2 rule (In[17]) [rule added; needs 3a registration to pass gate]
- [x] extend cheap guard; removed Erfc reduction (Route B owns it)

## Part 2 — Route C recognizers  [DONE, verified]
- [x] rec_besselk (In[9]) — numeric x-check 1.0, monomial √x works
- [x] rec_airy (In[12]) — anchor s=1→1/3, NIntegrate x-check 0.2887
- [x] wire into try_recognizers

## Part 3 — Route B + new special functions  [DONE, verified]
- [x] 3a. StruveH head (new module + deriv skipped-by-design + docstring + SYM_); In17 verified vs NInt 1.912
- [x] 3a. ExpIntegralE head (expint CF/series numeric, D[E_n]=-E_{n-1}, E_0); E1(2),E2(1.5) verified
- [x] 3b. rec_ibp IBP fallback (sv>0 lower + kernel-strip upper); In11/15/16 verified (Ci vs NInt -13.398)

## ALL 20 cases close end-to-end [verified]

## Tests / docs / verify
- [ ] Unit tests: 10 working + 10 now-closed + numeric anchors (∫Ai=1/3, ∫Erfc=1/√π, ∫E1=1)
- [ ] docs/spec/builtins/calculus.md + special-functions.md; changelog 2026-10-05.md
- [ ] version.h bump per substantive commit
- [ ] Build, run integrate_ramanujan_tests, re-run user's 20 cases, numeric-surface audits, check-messages, check-c99, valgrind

## Review

**Result: all 20 stress cases close (was 10). No regressions.**

Changes:
- `integrate.c`: `def_has_undefined_function` now x-aware (one early-out via
  `depends_on_var`); unblocks pFq (In6–8), repairs the red 1F1/2F1 tests.
- `integrate_ramanujan.c`: `discharge_strip` (Simplify→Refine, per-atom,
  all-or-nothing collapse); `rec_besselk`, `rec_airy`; `rec_ibp` IBP fallback;
  `reduce_to_hypergeometric` gained Si/StruveH rules (dropped Erfc — Route B
  owns it).
- New special functions: `src/special_functions/expintegrale.{c,h}` (expint
  CF/series numeric), `struveh.{c,h}` (1F2 numeric); SYM_ names; `ExpIntegralE`
  deriv rule in `deriv.c`; docstrings in `info.c`; init in `core.c`; CMake list.

Verification:
- `integrate_ramanujan_tests` green (21 tests incl. formerly-red 1F1/2F1, new
  `test_mellin_assumptions` + `test_mellin_special_functions`).
- No regressions: residue/beta/intrep/newton_leibniz/line/principalvalue/
  symmetry/deriv/unknown test suites all PASS.
- Caught + fixed a self-regression: Refine discharge had partially dropped a
  proven conjunct (sector `1/(1+x^n)`), making the ConditionalExpression strip
  incomplete — fixed to all-or-nothing collapse.
- Audits: check-c99, check-messages, check-packed-aware, check-array-exactness
  all green. Numerics cross-checked vs NIntegrate (Airy, Ci, StruveH) and known
  anchors (∫Ai=1/3, ∫Erfc=1/√π, ∫E₁=1).

Lesson captured: CMake test lib has an explicit source list — a new `src/**`
file must be added to `tests/CMakeLists.txt` (it does not auto-glob like the
makefile).
