# M63 — DSolve corpus §2.2.36 (3501–3600) + the parameter-as-indvar converter repair

Next hundred in the `DSOLVE_PLAN.md` campaign: upstream §2.1.36
(`Ch2.S1.SS36.htm`), Problems 3501–3600, 100 records / 100 scalar / 9 IVPs.

Converting it turned up a bug with reach far beyond the new section: the
converter promotes a lone **parameter** letter to the independent variable, so
`3570` (`y'' - 2a y' + a² y == 0`) became an ODE *in `a`*. Measured over all 36
committed corpora, **17 records are silent wrong equations of that class** — one
(`§2.2.16-1534`) already written down in `README.md` as an unexplained residue.
A mis-transcribed record scores PASS or UNEVAL against an equation the book never
asked, so this is a correctness bug in the gated corpora.

Two commits, each bumped and tagged (user's call).

## Commit 1 — converter fix + the 17-record repair

- [x] `detect_symbols`: restrict the step-3 lone-letter adoption to `y` (+ the
      existing Greek path) — every other lone Latin letter is a parameter
- [x] `detect_symbols`: under a `_missing_x` classification, take the indvar only
      from a function-argument position, else the fresh-letter fallback (needs
      `classif` plumbed through `convert_row` ← `main`)
- [x] Verify in isolation: old-vs-new converter on the SAME fetched page, **19**
      sections. Exactly 13 records move, all known victims; the other 10
      sections (incl. §2.2.1/20/21/24/25/29/30/35 as negative controls) are
      byte-identical
- [x] Patch the **17** record lines into the **9** committed `DE_examples_*.m`
      (record lines only — keep upstream LaTeX drift out of the diff). The 17th
      is §2.2.28-2789, found by the new audit, not by the hand scans
- [ ] Re-measure those 9 sections before/after, same machine, back to back;
      update `argv[3]` in `tests/CMakeLists.txt` + `reports/*` + `STATUS.md`
- [x] `tools/check_corpus_indvar.py` + `make check-corpus-indvar`: green on the
      repaired tree, 24 findings on the pre-fix copies
- [x] `README.md`: the stale `Ch2.S2.SSN.htm` fetch URL, the "cannot be
      regenerated" note, and an audit-first instruction
- [ ] v0.255, changelog, `DSOLVE_PLAN.md` M63 entry, tag

## Commit 2 — the §2.2.36 wave

- [x] Generate `DE_examples_2236.m`; validate (100 records, all parse, trap greps
      clean, indvar census x/z/t only, `check-corpus-indvar` green)
- [ ] `add_test(dsolve_corpus_2_2_36_tests)` + `STATUS.md` block + `README.md` row
- [ ] Baseline: fork-per-case, twice, per-case identical; bucket report
- [x] The `3521`/`3527`/`3599` root cause, **re-diagnosed by measurement**: these
      are NOT a missing-answer case. `DSolve\`Linearizable` already solves them
      from the original equation in 0.11 s, three cascade slots after
      `Separable` — which eats the whole 8 s budget first on an `Integrate` of
      its own SAMPLED integrand. So the fix is the mixed-angle `TrigExpand`
      normalisation in `sep_find_split`, applied to `F` BEFORE sampling (a retry
      cannot work: `TrigExpand` is a no-op inside a denominator), gated on a
      mixed-angle kernel and guarded on the rewrite actually eliminating it
      (§2.1.2-1134 matches the gate, cannot be helped, and must stay byte-identical)
- [x] `tests/test_dsolve_m63_stress.c`, five families, negative controls + a
      latency bound (the fix is a latency property, so an answer-only test would
      pass before AND after)
- [ ] Regression: isolated-worktree A/B over the touched sections + the §2.1.2
      master corpus; `check-c99`, `check-messages`, valgrind
- [ ] v0.256, docs, `STATUS.md` (also add M62's missing wave-history bullet), tag

## Review

(filled in as the work lands)
