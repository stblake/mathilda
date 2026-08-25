# Book: insert Chapter 2 "Compiling and Running Mathilda"

Insert a new chapter between Ch.1 (About) and the Introduction, renumbering the rest
so filename number = chapter position (the book's convention; cross-refs are all
label-based so they auto-update).

## Plan

- [ ] Renumber chapters 12→13 … 2→3 with `git mv` (descending to avoid collisions).
- [ ] `git mv` `examples/02-introduction` → `examples/03-introduction` and
      `figures/02-introduction` → `figures/03-introduction`; sed the intro file's
      `02-introduction/` paths → `03-introduction/` (\pair/\mtranscript/\plotcell).
- [ ] Fix the "next chapter" forward-reference in `01-about.tex`.
- [ ] Update the `\include` list in `TheMathildaBook.tex` (insert 02-building at pos 2).
- [ ] Add a reusable `shell` listing style + environment to `mathilda.sty`.
- [ ] Write `chapters/02-building.tex` (label `ch:building`). Sections:
      obtain the source; prerequisites & toolchain (GCC-not-clang, required vs
      optional libs, per-OS install, USE_* graceful-degradation table);
      building (make, leaner builds, common failures); first session & banner;
      how the REPL works (In/Out/%, ; suppression, multiline `\`, readline
      history, ?help/Names, Quit/Ctrl-D); scripts & run modes (-file, bare file,
      --help/--version, NDJSON pipe); building & running the tests; forward ptr.
- [ ] Add verified examples `examples/02-building/{version,first-session,back-reference,help}.m`
      (only pipe-round-trippable expressions; shell commands as non-verified listings).
- [ ] Concept `\index{}` entries; update `ROADMAP.md` + current-week changelog.
- [ ] Verify: `make examples`, `make check-links` (0 unlinked), `make pdf` (clean log);
      confirm Ch.2 = building, Ch.3 = introduction in the ToC.

## Review

Done and verified.

- New **Chapter 2 "Compiling and Running Mathilda"** written to
  `chapters/02-building.tex` (label `ch:building`); ToC confirms it sits at
  position 2 (About = 1, Introduction = 3, …, About the Author = 13, Appendix A).
- Renumbered chapters 2→3 … 12→13 with `git mv`; renamed the introduction's
  `examples/` and `figures/` dirs to `03-introduction/` and repointed its 151
  in-file paths. All cross-refs are label-based, so they auto-resolved.
- Fixed Chapter 1's forward-reference; it now points to Chapter 3 (verified in PDF).
- 6 build-verified transcripts in `examples/02-building/`
  (version, first-session, session-state ×, suppression, help). Discovered the
  build-time pipe does NOT populate `In[n]`/`Out[n]`/`%` history — those are
  interactive-REPL-only — so `%`/`Out[n]` are taught in prose, not faked in a
  transcript (saved as a memory).
- Added a reusable `shell` listing environment to `mathilda.sty` for the
  hand-authored terminal commands (git/apt/brew/make/`--help`); the `-file`
  demo output was captured from the real binary.
- `ROADMAP.md` renumbered; changelog entry added under
  `docs/spec/changelog/2026-08-24.md`.
- **Verification:** `make examples` clean; PDF builds (63 pp) with no undefined
  refs, no LaTeX errors, no multiply-defined labels, and no new overfull boxes
  from this chapter.

**Pre-existing issues (NOT introduced here), left as-is:**
- `make check-links` fails on `\B{Reduce}` (7 uses in the introduction chapter,
  from an earlier commit) — `Reduce` has no page in `builtins.json` yet.
- `make figures` needs a GUI session to render the introduction's plots; it
  fails headlessly. The book PDF still builds (plots degrade to placeholders).
