# Task: DSolve §2.2.20 + §2.2.21 corpus coverage (M40, M41)

Plan: `~/.claude/plans/let-s-continue-our-implementation-functional-papert.md`

## Phase 0 — Converter migration to LaTeXML (shared prerequisite)
- [ ] `strip_array`: strip `\begin{array}[]` optional-arg `[]`
- [ ] `parse_table`: add LaTeXML branch (auto-detect no `id='TBL-'`); header-row
      column mapping; multiline (`re.S`) array alttext; yield {n,tex,classif,sympy}
- [ ] Sanity: convert both pages to scratch, round-trip parse all 200 records

## Phase 0 result
- [x] Converter: strip_array `[]` fix + `_implicit_mult` (juxtaposition spacing)
      + LaTeXML `parse_table` branch. Both files convert 100/100, round-trip clean.

## Phase 1 — M40: §2.2.20 (Problems 1901–2000)  — baseline 96/100, 0 FAIL
- [x] Generate `DE_examples_2220.m` (66 general + 34 IVP)
- [x] Baseline: 96 PASS / 4 UNEVAL (1916,1941,1964,1976), 0 FAIL, 0 crash
- [ ] Register `dsolve_corpus_2_2_20_tests` (baseline 4)
- [ ] `reports/2.2.20.{md,tsv}`
- Residue = Kovacic churn (Case-1 Riccati ds_solve on regular-singular series
  ODEs, 9-18s → harness 8s timeout). Investigated a bounded fix; NO safe
  discriminant exists (the churn IS Liouvillian-existence, cf. 1822 which needs
  the same quadratic-pole term). Bounded declines, documented. NOT a wave.

## Phase 2 — M41: §2.2.21 (Problems 2001–2100)  — baseline 96/100, 0 FAIL
- [x] Generate `DE_examples_2221.m` (100 general)
- [x] Baseline: 96 PASS / 4 UNEVAL (2003,2005,2006,2080), 0 FAIL, 0 crash
- [ ] Register `dsolve_corpus_2_2_21_tests` (baseline 4)
- [ ] `reports/2.2.21.{md,tsv}`

## Phase 3 — Dashboards, docs, verification
- [x] `tests/CMakeLists.txt`: `dsolve_corpus_2_2_20_tests` + `_2_2_21_tests`, baseline 4 each
- [x] `reports/2.2.20.{md,tsv}` + `reports/2.2.21.{md,tsv}`
- [x] STATUS.md: §2.2.20 / §2.2.21 blocks + wave-history M40
- [x] README.md: DE_examples table rows + LaTeXML source note
- [x] DSOLVE_PLAN.md: M40 entry (incl. Kovacic-churn investigation record)
- [x] Changelog `docs/spec/changelog/2026-09-07.md`; version 0.137→0.138 (src/version.h)
- [x] (no new builtin — no docs/spec/builtins change needed)
- [x] Memories: LaTeXML migration + Kovacic-churn-no-safe-gate (+ MEMORY.md index)
- [x] Verify: both new ctests PASS (baseline 4); 2218/2219 regression clean; check-c99 green; 0 FAIL

## Review

**Outcome.** Two corpus sections added: **§2.2.20 (1901–2000) 96/100** and
**§2.2.21 (2001–2100) 96/100**, both **0 FAIL, 0 crash**. Consistent with the
series-heavy §2.2.5 (99) and §2.2.18/19 (96/95).

**Key work.** The site migrated tex4ht→LaTeXML, breaking the converter (found 0
cells). Ported `latex_ode_to_mathilda.py` to the LaTeXML layout (new `parse_table`
branch + `_implicit_mult` juxtaposition-spacing + `strip_array` `[]` fix); LaTeX→Mathilda
core unchanged; 200/200 records convert + round-trip.

**Kovacic churn (not fixed — deliberately).** The 8 non-PASS are non-Liouvillian
regular-singular ODEs whose Kovacic Case-1 solve churns 9–18 s > the 8 s harness
window before Frobenius (correct series) runs. Tried a nested-`TimeConstrained`
bound (no-nest hazard, >90 s regression) and structural pre-gates (all regress the
pinned `2.2.19-1822` / `t_m39` — the churn *is* Liouvillian-existence, no cheap
discriminant). Reverted to 0-regression clean main; documented as bounded declines
(the M39/1823 class). Memory written so it isn't re-attempted.

**No C/behavior change** (Kovacic reverted). Deliverable = converter + 2 corpora +
2 ctests + docs. Version 0.138.
