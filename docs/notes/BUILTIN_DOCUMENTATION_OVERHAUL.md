# Builtin Documentation Overhaul

A multi-campaign effort to bring **every** Mathilda builtin's documentation page
up to the detail level of [`EllipticF`](https://stblake.github.io/mathilda/documentation/special-functions/EllipticF/)
— the gold standard — with **every example run through the latest Mathilda
binary**.

This file is the living tracker: the parity bar, the authoring recipe, the
verification workflow, and the per-category campaign checklist with a progress
table. It is contributor-facing prose, so edits to it alone do **not** bump
`$VersionNumber` and are **not** tagged.

---

## 1. What "EllipticF-grade" means

The public site (`site/`) is MkDocs Material. Every builtin already gets an
auto-generated page at `site/docs/documentation/<category>/<Name>.md` — **never
hand-edit those**; they are a build product of `site/generate.py` (`make docs`).
Richness comes from four *inputs* the generator fuses in:

| Layer | File / source | Renders as |
|-------|---------------|-----------|
| **Description** | the C `symtab_set_docstring()` string | `## Description` |
| **Examples** | `In[]:=`/`Out[]=` blocks in `docs/spec/builtins/<cat>.md` | `## Examples` (verified) |
| **Worked examples + Notes** | `site/overlays/<Name>.md` | `### Worked examples` (lifted into Examples) + `## Notes & additional examples` |
| **Implementation note** | `site/impl/<Name>.md` | leads `## Implementation notes` |

A page is **EllipticF-grade** when all three hand-authored layers exist and the
build is clean:

1. **`site/impl/<Name>.md`** — source-grounded. Front matter `source:` (the real
   module path, overriding the docstring hub) and `references:` (a list of
   citations). Body of three bold-led paragraphs:
   `**Algorithm.**` · `**Data structures.**` · `**Complexity / limits.**`.
   Models: `site/impl/EllipticF.md`, `site/impl/Apart.md`.
   *Non-algorithmic symbols* (colors, option tokens, type tokens, constants) use
   the same three headings, reframed: *Algorithm* → what the symbol is and where
   it is defined in source; *Data structures* → its internal `Expr` form;
   *Complexity / limits* → evaluation/usage behaviour (often "inert symbol").
2. **`site/overlays/<Name>.md`** — `### Worked examples` + `### Notes`.
   Models: `site/overlays/EllipticF.md`, `site/overlays/Zeta.md`.
3. **Examples in `docs/spec/builtins/<cat>.md`** — a handful of basic cases plus
   application cases, auto-verified by `make docs`.

Optional, added only where real: `## Performance` (needs `site/perf.json` or a
benchmark row) and `## Options & behaviour` (only if the function takes options).

---

## 2. Authoring recipe (per function)

1. Find the source module (the page's `Source:` link / the `symtab_set_docstring`
   call / an impl `source:` override). **Read the C or `.m` source** — never guess
   the algorithm.
2. Write `site/impl/<Name>.md` from the source: the real algorithm, the data
   structures it uses, its complexity and limits, with citations in front matter.
3. Write `site/overlays/<Name>.md`:
   - `### Worked examples`: several `In[1]:=` inputs. **You write only the
     input** — the generator runs it through the current binary and fills the
     `Out[]`. Each input may carry a trailing `(* ... *)` note (lifted to a
     sentence above the cell).
   - `### Notes`: prose on semantics, gotchas, domain, relations to other heads.
4. Ensure `docs/spec/builtins/<cat>.md` has basic + application `In/Out` examples
   for the function (add them if the section is thin).
5. Regenerate and verify (§3). Spot-check the rendered page.
6. Tick the box in §4; add a note to the current week's
   `docs/spec/changelog/<Monday>.md`.

### Authoring traps (all enforced or bite silently)

- **Overlay `In[]` must be single-line.** The binary is driven over a one-line
  NDJSON protocol; a wrapped input truncates and the definition silently never
  binds.
- **A `(* ... *)` note's first word is capitalised** by the renderer — do not
  start a note with an identifier (`n`, `DK`, …); lead with an ordinary word.
- **Prose between overlay code blocks is dropped.** Per-example commentary goes
  *inside* the `(* ... *)` note; general prose goes in `### Notes`.
- **Never hand-edit `site/docs/documentation/**`.** Edit the inputs, regenerate.
- Prefer stable outputs in examples (avoid gratuitous tiny-magnitude or
  scientific-notation results where an exact form exists).

### Scaling with subagents

Dispatch per-function or per-small-batch authoring to subagents: each reads the
module source and drafts `impl/<Name>.md` + `overlay/<Name>.md` for its slice;
the main session runs the generator + gate and reviews. One function group per
subagent.

---

## 3. Verification workflow

```bash
make                         # build the latest ./Mathilda (the generator needs it)
make docs                    # regenerate + verify EVERY spec and overlay example
                             #   against the binary (overlay Out[] is binary-supplied)
make docs-build              # strict mkdocs build (fails on warnings)
make check-docs-examples     # backstop: committed pages still match the binary
python3 site/coverage_report.py          # per-category parity table (progress metric)
python3 site/coverage_report.py <slug>   # one category, listing the gaps
python3 site/coverage_report.py --missing # every function missing a layer
make docs-serve              # local preview at :8000
```

`site/generate.py` runs every overlay worked-example input through the binary and
supplies the output, so an overlay output can never silently rot. `make
check-docs-examples` (`site/verify_docs_examples.py`) re-extracts every published
`In/Out` pair and re-runs it, catching a tree that went stale because the binary
changed and nobody re-ran `make docs`.

Commit scope: stage only `site/**`, `docs/spec/**`, and this tracker. Docs-only
commits do **not** bump `$VersionNumber` and are **not** tagged.

---

## 4. Campaign roadmap

Order: **math-core first**, then non-core partials, then the bare blocks. Counts
below are the baseline snapshot (regenerate with `coverage_report.py`). Tick a
category when `coverage_report.py` reports it at 100% graded.

### Campaign 0 — infrastructure
- [x] Generator auto-verifies overlay worked-examples (binary supplies `Out[]`)
- [x] `site/coverage_report.py` (progress metric)
- [x] `site/verify_docs_examples.py` + `make check-docs-examples` (local gate)
- [ ] CI wiring deferred — a hard binary-diff gate needs a CI job with **all**
      optional deps (gmp, mpfr, readline, ecm, lapack, fftw, pcre2, raylib) so
      outputs reproduce byte-for-byte; the current Linux job omits several and
      would false-fail on graphics/fourier/lapack/string examples. Run the gate
      locally pre-commit until that job exists.
- [x] This tracker
- [x] Pilot: `elementary-functions` to 100%

### Math core — all 100%
- [x] C1 elementary-functions
- [x] C2 arithmetic
- [x] C3 calculus
- [x] C4 number-theory
- [x] C5 special-functions
- [x] C6 linear-algebra
- [x] C7 comparisons · simplification · power-series · solutions-of-equations · mathematical-constants

### Non-core partials — all 100%
- [x] expression-information · data-structures · functional-programming ·
      structural-manipulation · control-flow · assignment-and-rules ·
      scoping-constructs · pattern-matching · string-operations · statistics ·
      numerical-calculus · flint · file-io · lists-and-iteration · time-and-date ·
      random-number-generation

### Bare blocks — all 100%
- [x] graphs · other-advanced · image-processing · hypergraphs · graphics ·
      machine-learning · fourier-transforms · packed-arrays · geometry · bitwise

---

## 5. Progress

Regenerate this table with `python3 site/coverage_report.py`.

- Baseline @ v0.266 (2026-10-04): **369 / 1086 pages EllipticF-grade (33%)**.
- **Complete @ 2026-10-04: 1086 / 1086 (100%)** — every builtin page now has a
  worked-example + Notes overlay, a source-grounded implementation note, and
  verified examples. All 38 categories at 100%. `make docs` verifies 9150
  examples; `make check-docs-examples` is green on all 1086 pages; strict
  `mkdocs build` is clean.

The authoritative, always-current breakdown is `coverage_report.py`'s output,
not this snapshot.
