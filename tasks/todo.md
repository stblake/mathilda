# Task: Review & broaden Book §4.1 "Arithmetic" (+ document 5 arithmetic builtins)

Plan: /Users/user/.claude/plans/let-s-review-section-4-1-bubbly-magpie.md

## Part 1 — Document 5 arithmetic builtins & publish reference pages
- [ ] 1a. Add docstrings in src/info.c: Mod, Quotient, QuotientRemainder (near GCD cluster), IntegerPart, FractionalPart (after Round)
- [ ] Rebuild `make -j`; sanity `?Mod` etc.
- [ ] 1b. site/impl/<Name>.md (5) + site/overlays/<Name>.md (4 new; Mod exists) + generate.py CATEGORY_OVERRIDES (IntegerPart, FractionalPart)
- [ ] 1c. python3 site/generate.py; curated-diff revert of drift; mkdocs --strict clean
- [ ] 1d. cd book && make links → \B{} resolves for the 5

## Part 2 — Broaden & polish book/chapters/math/arithmetic.tex
- [ ] A. Fix trim fallout (intro wording, transition ending "types of number")
- [ ] B. New subsection: Complex numbers (+ examples/arithmetic/complex.m)
- [ ] C. New subsection: Rounding and integer division (+ rounding.m, division.m)
- [ ] D. New subsection: Rational parts and cleanup (+ parts.m)
- [ ] E. Interval: IEEE-1788 comparison (+ references.bib ieee1788)
- [ ] F. Prose fixes (15.95 derivation; N[1001!] wide-exponent note)
- [ ] G. Closing forward-pointing paragraph
- [ ] Index concept entries for new material

## Verification
- [ ] make examples green (adjust inputs, never outputs)
- [ ] make check-links = 0
- [ ] make pdf clean (makeindex reruns)
- [ ] changelog docs/spec/changelog/2026-08-24.md

## Review — COMPLETE

Root cause found: `Mod`, `Quotient`, `QuotientRemainder`, `IntegerPart`,
`FractionalPart` were registered but had **no docstring**, so the site (which
discovers pages by parsing `symtab_set_docstring`) never generated pages and `\B{}`
could not link them. Fixed the docstrings → regenerated the site → made them book-linkable
→ broadened §4.1 to use them.

Done & verified:
- **C** (`src/info.c`): 5 docstrings added; behaviour unchanged (attributes + fast paths
  already present). Semantics re-verified at the REPL (floored Mod/Quotient, offset
  3-arg, toward-zero IntegerPart, signed FractionalPart, identity n·Q+Mod=m).
- **Site**: 5 impl + 4 overlay inputs, 2 CATEGORY_OVERRIDES; `python3 generate.py`
  produced 5 new pages + 5 builtins.json entries; drift reverted to a curated diff;
  `mkdocs build --strict` rc=0.
- **Book link table**: `make links` → all 5 now resolve.
- **Book §4.1**: intro/transition trim fallout fixed; 3 new subsections (Complex
  numbers; Rounding & integer division; Rational parts & cleanup); IEEE 1788 comparison
  + `ieee1788` bib entry; prose fixes (15.95 derivation, wide-exponent N[1001!]); closing
  paragraph. 4 new verified example files.
- **Verification**: `make examples` (61 transcripts green), `make check-links` (0
  unlinked), `make pdf` clean (no undefined refs / no "No reference link"; ieee1788
  cited; 6 new concept index entries pp.37–41; TOC shows 4.1.1–4.1.8). Widest overfull
  hboxes are pre-existing (other chapters).
- **Docs**: changelog `docs/spec/changelog/2026-08-24.md`; ROADMAP §4.1 scope refreshed.

Note (not actioned): CONTEXT §6 lists callouts as boxed environments, but the
implementation is the `\calloutnote` footnote macro — a one-line CONTEXT fix, separate.
