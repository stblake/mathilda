# Book Campaign 0 — Foundations + Chapter 1  ✅ DONE

Plan: `/Users/user/.claude/plans/in-book-i-would-cheeky-raven.md`

## Toolchain
- [x] `book/tools/gen_links.py` — builtins.json → generated/builtinlinks.tex (833 links)
- [x] `book/tools/build_examples.py` — examples/**/*.m → transcripts (reuses run_session)
- [x] `book/tools/check_links.py` — fail on any \B{} with no reference page
- [x] `book/mathilda.sty` — REPL transcript style, `\B{}` link macro, callout boxes
- [x] `book/Makefile` — links, examples, pdf, check-links, clean/distclean
- [x] `book/.gitignore`

## Book source
- [x] `book/TheMathildaBook.tex` — memoir master, parts + \include
- [x] `book/frontmatter/` — titlepage, copyright (GPLv3), preface, colophon
- [x] `book/chapters/01-about.tex` — Chapter 1 pilot (fully written & verified)
- [x] 17 stub chapters for every later section/appendix
- [x] `book/examples/01-about/*.m` — 3 verified example sessions
- [x] `book/references.bib`

## Docs
- [x] `book/CONTEXT.md` — principles (verified-example promise, hyperlinks, voice, stats policy)
- [x] `book/ROADMAP.md` — every section as a future campaign + status table

## Verify
- [x] gen_links spot-check: Integrate/Sin/$Version(%24)/FLINT`Det URLs correct
- [x] make examples → 3 transcripts; byte-identical to `verify_tutorial.py transcript`
- [x] make pdf (clean-slate distclean→pdf) → exit 0, 32 pp., 0 undefined/unlinked warnings
- [x] 12/12 Chapter-1 builtin hyperlinks embedded in PDF (decompressed & confirmed)
- [x] make check-links → OK (15 uses, 12 distinct, all resolve)
- [x] changelog note (docs/spec/changelog/2026-08-17.md)

## Review

Delivered the full book foundation and a verified pilot Chapter 1. Key decisions
(confirmed with the user): LaTeX/memoir; monospace REPL output; build-time-generated
examples (inputs-only source, outputs injected by running `./Mathilda` — no output
can be hand-written or drift); pilot includes a fully-written Chapter 1.

Mechanism highlights:
- Examples: one `.m` per Mathilda session under `examples/`; `build_examples.py`
  reuses `site/verify_tutorial.py::run_session` so book and site stay consistent.
- Hyperlinks: `\B{Name}` resolves via a table generated from `builtins.json`
  (the single source of truth); `check_links.py` gates any dangling `\B{}`.
- One real bug found & fixed during build: `%24` percent-encoding of `$`-symbol URLs
  contains `%`, which LaTeX reads as a comment — now escaped to `\%24` for `\href`.
- Name gotchas caught by check-links: `SVD`→`SingularValueDecomposition`; no
  `Graphics` page (used `Show`).

Not in scope (future campaigns, see ROADMAP.md): all chapters beyond 1; graphics
example capture (plots return image/plot messages, not text — needs a build_examples
extension); Data I/O last (mostly unimplemented).
