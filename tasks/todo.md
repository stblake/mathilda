# Documentation refresh: numerical calculus + special functions — 2026-06-14

Goal: tutorials for the numerical-calculus routines and for the special
functions, list both pillars in the main-page "What's inside", and bring all
generated docs up to date with every built-in.

## Plan

- [ ] 1. Spec file `docs/spec/builtins/numerical-calculus.md` — H2 sections for
      `ND`, `NIntegrate`, `NSum`, `NProduct`, `NLimit`, `NSeries`, `NResidue`
      (drives a new site category; was falling into "Other & Advanced").
      Include verified `In[]/Out[]` example blocks.
- [ ] 2. `docs/spec/builtins/special-functions.md` — add a `## LogGamma` H2 so
      `LogGamma` is filed under Special Functions (was Other & Advanced).
- [ ] 3. `Mathilda_spec.md` — add a "Numerical calculus" row to the built-in
      reference table.
- [ ] 4. Changelog note in `docs/spec/changelog/2026-06-08.md`.
- [ ] 5. New tutorial `site/docs/tutorials/07-numerical-calculus.md`.
- [ ] 6. New tutorial `site/docs/tutorials/08-special-functions.md`.
- [ ] 7. `site/docs/tutorials/.pages` + `index.md` — register both tutorials.
- [ ] 8. `site/docs/index.md` — add "Numerical calculus" and "Special functions"
      cards to "What's inside"; refresh function count.
- [ ] 9. Regenerate the site: `python3 site/generate.py`.
- [ ] 10. Verify generated output + tutorial transcripts.

All example transcripts captured by running `./Mathilda` directly.

## Review — DONE 2026-06-14

All 10 steps complete and verified.

- **New doc category.** `docs/spec/builtins/numerical-calculus.md` now owns
  `ND`, `NIntegrate`, `NSum`, `NProduct`, `NLimit`, `NSeries`, `NResidue`. The
  rich authoritative sections were **moved** out of `calculus.md` (not copied —
  no duplication); `calculus.md` keeps the symbolic routines + FindRoot/FindMin/
  FindMax. Added verified `In[]/Out[]` example blocks to the NProduct/NIntegrate
  sections (they had none) and fixed a line-wrapped `NSum` example that blocked
  example mining. Registered the category in `Mathilda_spec.md`.
- **Special functions.** Promoted `LogGamma` from `###` to `##` in
  `special-functions.md`; the category now generates 21 pages (was 6 — stale).
- **Tutorials.** Added `07-numerical-calculus.md` (19/19 transcripts match the
  build) and `08-special-functions.md` (34/34). Registered both in
  `tutorials/.pages` and the tutorials index.
- **Landing page.** Added "Numerical calculus" + "Special functions" cards to
  "What's inside"; bumped the count to ~435.
- **Regeneration.** `python3 site/generate.py` → 435 pages, 1111 verified
  examples, **0** in "Other & Advanced". 24 categories. Doc centre landing,
  per-category indexes, `.pages` nav, and `builtins.json` all refreshed.
- **Validation.** `mkdocs build --strict` passes (no broken links). Changelog
  note added to `docs/spec/changelog/2026-06-08.md`.

No `src/` changes — documentation only. `src/external/` untouched.
