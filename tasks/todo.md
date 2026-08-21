# Task: Advanced Diophantine-equations tutorial (multi-page section)

Plan: `~/.claude/plans/we-need-to-add-pure-plum.md`. Section placed right after
`07-solutions-of-equations.md`, modelled on `integration-methods/`.

## Steps
- [x] Build `./Mathilda` (v0.075 engine) — already current
- [x] Read primary sources: SOLVE_INTEGERS.md, benchmarks 87 + 88 reports
- [x] Create folder `site/docs/tutorials/diophantine-equations/`
- [x] Write `index.md` (landing: domain, two guarantees, card grid)
- [x] Write `linear-and-pell.md`
- [x] Write `quadratic-forms.md`
- [x] Write `cubics-and-thue.md`
- [x] Write `exponential.md`
- [x] Write `famous-power-sums.md`
- [x] Write `performance.md` (vs sympy exp.87, vs PARI/GP exp.88)
- [x] Create subfolder `.pages`
- [x] Register: edit `site/docs/tutorials/.pages` (after 07)
- [x] Register: add card to `site/docs/tutorials/index.md`
- [x] Generate every `Out[]` via `verify_tutorial.py transcript`; NO timings in checked blocks
- [x] Verify each page: `verify_tutorial.py check <page>` → OK (all 7)
- [x] `make docs-build` (mkdocs --strict) → 0 warnings, section built
- [x] Changelog note in docs/spec/changelog/2026-08-17.md

## Review

Added a 7-page advanced tutorial section
`site/docs/tutorials/diophantine-equations/` (index + 6 pages), placed after
`07-solutions-of-equations.md`, mirroring the `integration-methods/` model.

Verification:
- All 48 checkable `In[]/Out[]` transcripts pass `verify_tutorial.py check`
  against the live binary (index 1, linear-and-pell 14, quadratic-forms 8,
  cubics-and-thue 11, exponential 4, famous-power-sums 10; performance is
  tables/prose only). One paste typo (extra `}`) caught and fixed.
- `mkdocs build --strict` completes with **0 WARNING-level messages** — every
  relative link (siblings, `../07-...`, `../../documentation/...`) resolves. The
  only INFO notes are pre-existing broken links in the auto-generated
  `documentation/` reference pages, unrelated to this change.

Content: each family page follows the house advanced-tutorial skeleton
(numbered `##` sections, `### References`, `\(...\)`/`\[...\]` math, `Next:`
link). Famous results covered with exact runnable syntax: Lander–Parkin, Frye,
taxicab 1729, Euler brick, Ramanujan–Nagell, Booker three-cubes, FLT, Catalan,
Markov, Brocard. Performance page reproduces the two existing benchmark suites
(vs sympy: 19/19, 20–100× faster, 11 sympy NotImplementedError; vs PARI/GP:
0-wrong parity + one PARI-incompleteness catch, honest speed/coverage gap).
Timings live only in prose/tables (never in verified `Out[]` blocks) since
`Timing[]` is nondeterministic.

Out of scope (flagged to user): the generated
`documentation/solutions-of-equations/Solve.md` still describes the *old*
`Integers` post-pass and understates the engine — a separate docstring/regen
task.

Preview locally: `make docs-serve` → http://127.0.0.1:8000/tutorials/diophantine-equations/

---

## Follow-up (2026-08-21): mod-9 impossibility + Solve.md refresh

Two follow-on requests from the user.

### A. Sum-of-three-cubes mod-9 impossibility (engine)
- [x] `si_solve_three_cubes_mod9` in `src/solve/solveint_cubes.c`; prototype in
      `solveint_internal.h`; dispatched after the FLT short-circuit in
      `solveint.c`. Returns `{}` for `x^3+y^3+z^3==k` with `k ≡ ±4 (mod 9)`, no
      bound needed; gated to unit `±1` cube coefficients.
- [x] Unit tests `test_three_cubes_mod9` (14 cases: k≡4/5, negative k, signed
      variants, bounded box; guards for reachable residue, squares, coeff≠±1).
      Full `solve_integers_tests` suite green.
- [x] gcc build clean (incl. `-Werror=unused-function`); `check-c99` clean;
      `make check-diophantine-heldout` → "No silent wrong answers".
- [x] Tutorial `cubics-and-thue.md` §2 rewritten: mod-9 global proof (`==4 → {}`)
      + Booker bounded search; unbounded reachable residue (`==3`) still declines.
      Re-verified (12/12).

### B. Solve.md reflects the integer engine (docs)
- [x] Expanded the `Solve` docstring (`src/solve/solve.c`) to describe the
      Diophantine engine (HNF, Pell, quadratic/ternary forms, Mordell, Thue,
      three-cubes+mod9, exponential, power-sum MITM; `{}`=proof, else decline).
- [x] `docs/spec/builtins/solutions-of-equations.md`: scoped the univariate
      "post-pass" paragraph; added the mod-9 method to the three-cubes section.
- [x] Changelog entry (docs/spec/changelog/2026-08-17.md).
- [x] Regenerated `site/docs/documentation/` via `make docs` (writes the mkdocs
      Solve.md; NOTE: `frontend/public/refpages/.../Solve.md` is a separate
      generator, left untouched).
- [ ] Re-run `mkdocs build --strict` after regen; confirm Solve.md updated.
