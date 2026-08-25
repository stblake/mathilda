# The Mathilda Book — Principles and Conventions

This file is the constitution for the book in `book/`. Read it before writing or
editing any chapter. It exists so that every contributor (human or model) writes
to the same standard, and so that the book's load-bearing promises — above all,
that every example is real — cannot quietly erode.

---

## 1. What this book is

The book has **two purposes at once**, and every chapter must serve both:

1. A **user's manual** for Mathilda — practical, example-driven, teaching a reader
   how to actually use the system.
2. A **textbook on computational mathematics / computer algebra** — showing *how
   advanced mathematics is done inside a CAS*: the algorithms, the reasoning, and
   the trade-offs beneath the surface, not just the syntax.

If a passage teaches only syntax, it is half-written. Wherever it is natural, show
the reader *what the system is doing and why* — the mathematics it embodies and
the engineering that makes it fast. Callout boxes exist for exactly this (§6).

### Audience

Numerate readers: students, researchers, working scientists, curious amateurs.
Assume comfort with undergraduate mathematics; do **not** assume prior experience
with a computer algebra system, nor with Mathematica. Part II chapters may assume
the vocabulary of Part I but should otherwise be independently readable.

---

## 2. Positioning (how the book relates to the other docs)

Mathilda already ships a large documentation corpus. The book does **not** duplicate
it — it complements it:

- **Per-builtin reference pages** (`site/docs/documentation/<category>/<Name>/`) are
  the exhaustive truth about each function's options and corner cases. **The book
  links to them** (§5) rather than restating them.
- **Site tutorials** (`site/docs/tutorials/*.md`) are short, guided paths. The book
  is the **cohesive, long-form narrative** — a book you can read cover to cover —
  and goes deeper into the mathematics and the internals than a tutorial does.
- **Design docs** (`docs/design/`) and **SPEC.md** are the engineering record; the
  book draws on them for its "Under the hood" material but is written for a reader,
  not a maintainer.

When in doubt: the reference *enumerates*, the tutorials *onboard*, the book
*explains and teaches*.

---

## 3. The verified-example promise (non-negotiable)

> **Nothing printed as `Out[...]` in this book is written by a human. Every example
> is produced by running its input through the real Mathilda binary at build time.**

This is enforced mechanically, not by discipline:

- Example **inputs** live in `book/examples/<chapter>/<name>.m` — one input
  expression per line; blank lines and `#` comments are ignored. **There is no place
  in the book source to type an output.**
- `book/tools/build_examples.py` runs each `.m` file through `./Mathilda` (via the
  NDJSON pipe, reusing `site/verify_tutorial.py::run_session` — the same driver that
  verifies the site tutorials) and writes the verified transcript to
  `book/generated/<chapter>/<name>.tex`.
- The prose includes it with `\mtranscript{<chapter>/<name>}`, which reads the
  generated file **verbatim** (via `listings`), so no output can drift or be
  mistyped, and OutputForm text needs no LaTeX escaping.
- **One `.m` file = one Mathilda session.** State carries across its lines, so an
  example may assign on one line and use the value on the next. Group inputs into a
  file precisely to control session scope; start a new file for a fresh session.

Consequences to write around, not against:

- **Output form is Mathilda's, not a textbook's.** It prints `1/6 Pi^2`, may reorder
  terms, writes `E^x` for `e^x`, and an indefinite integral may differ from a table
  by a constant. All correct; teach the reader to read it. Never silently "prettify"
  an `Out[]` — if the math form aids understanding, typeset it in the *prose* (§6)
  and let the transcript stand as-is.
- **`;`-suppressed / `Null` lines** show `In[]` with no `Out[]` (setup lines).
- **Graphics** (`Plot`, `Graphics`, `Show`) return `plot`/`image` messages over the
  pipe, not text. The Graphics chapter's campaign will define how figures are
  captured and embedded; do not expect a text transcript for a plot.

To add or change an example: edit/add the `.m` file, run `make examples` (or
`make pdf`), and read the generated transcript. If you don't like the output, change
the **input**, never the output.

---

## 4. Format: LaTeX

The book is **LaTeX** (`memoir` class) → PDF. Rationale: print-quality math,
real cross-references, a real index, and the ability to typeset derivations
alongside the monospace transcripts. Source of truth for structure is
`TheMathildaBook.tex`.

- **Transcripts** are shown as **monospace REPL text** — exactly what the REPL
  prints — never re-typeset as math. (The pipe also returns a `latex` field; we
  deliberately do not use it for `Out[]`, to preserve fidelity to the real session.)
- **Code beside its explanation (the default layout).** Worked examples pair each
  `In[k]/Out[k]` with its own note, side by side, so it is always clear which words
  explain which line. Use the `codepairs` environment with one `\pair{<snippet>}{<note>}`
  per input, where `<snippet>` is `<chapter>/<name>/<k>` (the per-pair files
  `build_examples.py` emits alongside the full transcript). The columns wrap and
  **break across pages** (built on `paracol`), the note column is 2pt smaller than the
  body (`\explanationfont`), and the first code line aligns with the first line of its
  note. Keep notes to a sentence or two. `\codecomment{<path>}{<text>}` is the older
  block-level variant (one note beside a whole transcript) — prefer `codepairs`.
- **Mathematics in prose** uses ordinary LaTeX math (`\( \)`, `\[ \]`, `amsmath`).
- Body font Latin Modern; transcripts in Inconsolata.
- Do not introduce a dependency on packages missing from a base TeX Live install
  without checking (`mdframed` and `enumitem` are **absent** here; `tcolorbox`,
  `framed`, `biblatex`, `memoir`, `microtype` are present).

---

## 5. Hyperlink every builtin

The **first time** a chapter mentions a builtin in discussion, write it with `\B{}`:

```latex
\B{Integrate} computes an antiderivative...
```

`\B{Name}` typesets the name as code and hyperlinks it to its reference page. The
link table (`generated/builtinlinks.tex`) is generated by `book/tools/gen_links.py`
from `site/docs/assets/builtins.json` — the same index the site ships — so a name is
linkable **iff** it has a reference page.

- Use the **exact builtin name** as it appears in the index. Some intuitive names
  are not the registered ones: it is `SingularValueDecomposition`, not `SVD`; there
  is no `Graphics` page (use `Plot`/`Show`). `make check-links` (`check_links.py`)
  fails the build-check on any `\B{}` that resolves to nothing — run it.
- For `$`-system symbols, write `\B{\$Version}` (the backslash keeps `$` out of math
  mode; the generated key is escaped to match).
- Inside the transcripts themselves, names are plain monospace (not linked) — the
  transcript is a faithful copy of the session, not annotated.
- Don't over-link: link on first substantive mention in a chapter, not every
  occurrence, and not in headings (the macro is fragile in moving arguments).

### The Index

The book has a back-of-book Index (`\printindex`, standard `makeindex`). Two rules:

- **Builtins index themselves.** `\B{}` auto-emits `\index{Name@\mcode{Name}}`, so a
  builtin is indexed **in code font** (the `\mcode` monospace used everywhere else)
  wherever you reference it — no manual `\index` for builtins. Because the house
  style is to `\B` only the first substantive mention per chapter, this yields one
  clean entry per builtin per chapter. `\usagebox` additionally marks the builtin's
  definition card as the **bold** page (`\bidxmain`). If you ever need the link
  without an index entry (a boxed title, a fragile spot), use `\Bnoidx{}`.
- **Concept terms are indexed by hand.** For the ideas a reader looks up — *interval
  arithmetic*, *numeric contagion*, *arbitrary precision*, an attribute, a named
  algorithm — add `\index{topic!subtopic}` at the **defining** mention only, in the
  hierarchical `topic!subtopic` style (e.g. `\index{interval arithmetic!outward
  rounding}`, `\index{attribute!Flat}`). Index the primary discussion, not every
  passing use.

"Regenerate the Index" = add/refresh the concept `\index{}` entries for whatever new
material you wrote, then rebuild (`make pdf` reruns `makeindex`). The macros live in
`mathilda.sty` under *index helpers*.

---

## 6. Voice, tone, and structure

- **Second person, warm but precise** ("type the `In[...]` lines yourself and you
  will see the same `Out[...]`"). This matches the site tutorials.
- **Honest about limitations.** If Mathilda is slower than, or lacks a feature of,
  another system, say so plainly. The book earns trust by not overclaiming.
- **Comparative where it illuminates** — Mathematica, NumPy, SymPy, PARI/GP — but the
  book is about Mathilda; comparisons serve the explanation, they are not a scorecard.
- **Show the mathematics and the machinery.** Prefer explaining an algorithm to
  merely invoking it.

Callout environments (defined in `mathilda.sty`) — use them to separate the two
voices of the book from the main thread:

| Environment      | For |
|------------------|-----|
| `theory`         | The mathematics behind what just happened (a theorem, a method). |
| `underhood`      | How Mathilda implements it (data structures, algorithm, C internals). |
| `performance`    | Cost, scaling, when it is fast/slow, packed-array or `Compile[]` notes. |
| `pitfall`        | A trap, a surprising output form, a common misreading. |

Theorem-like environments (`theorem`, `definition`, `example`, `remark`, …) are
available for formal statements.

- Punctuation: US register; `--` for parenthetical dashes; no emoji in the prose.
- Every chapter opens with a short orienting paragraph and closes by pointing forward.

---

## 7. Reproducibility and versioning

- A build's transcripts reflect **one** Mathilda version, and the title and copyright
  pages name it. That version is **not hard-coded**: `book/tools/gen_version.py` reads
  it from the real `./Mathilda --version` at build time and `\def`s `\mathildaversion`
  (in `generated/version.tex`, loaded by `mathilda.sty`), so it can never drift from
  the binary that actually produced the transcripts. To change it, rebuild Mathilda —
  the next `make pdf` picks up the new version automatically. Still note a version
  change in the changelog. (`$Version` reports the full library set.)
- Because examples are re-run from source, building against a newer Mathilda simply
  shows what that version prints — the book stays honest automatically. If an output
  changes in a way that breaks the surrounding prose, fix the **prose**.

### Single source of truth for statistics

System statistics (lines of code, module counts, builtin counts) **drift** across
`README.md`, `SPEC.md`, and the site, and must not be copied casually. Policy:

- Cite the **reference index** for the number of documented symbols: 833 entries in
  `site/docs/assets/builtins.json` at the time of writing (this counts `$`-symbols
  and `FLINT`-context routines as well as ordinary builtins).
- Use README's user-facing framing — **"780+ built-in functions across 36
  categories"** — when a round figure is wanted in prose.
- For LoC / module counts, **re-derive** at writing time rather than quoting a stale
  number; if a precise figure matters, state how it was counted.

---

## 8. Building and checking

From `book/`:

| Command            | Effect |
|--------------------|--------|
| `make pdf`         | Regenerate links + transcripts, then build `TheMathildaBook.pdf`. |
| `make links`       | Regenerate `generated/builtinlinks.tex`. |
| `make examples`    | Run every `examples/**/*.m` through `./Mathilda`. |
| `make check-links` | Fail if any `\B{}` names a builtin with no reference page. |
| `make clean` / `distclean` | Remove build artifacts (and, for distclean, `generated/`). |

`generated/` is **git-ignored** — it is build output. The PDF build therefore
requires the Mathilda binary (`make examples` builds it if needed) and a TeX
toolchain (latexmk, pdflatex, biber, makeindex).

---

## 9. Contributing a chapter (checklist)

1. Read this file and skim the target chapter's entry in `ROADMAP.md`.
2. Write the prose in `chapters/NN-slug.tex`. Link builtins with `\B{}` on first use.
3. Put every worked example's **inputs** in `examples/NN-slug/<name>.m`; include the
   transcript with `\mtranscript{NN-slug/<name>}`. Never write an `Out[]` by hand.
4. `make examples` and read the generated transcripts; adjust the *inputs* until the
   session tells the story you want.
5. `make check-links` (0 unlinked) and `make pdf` (clean log: no undefined refs, no
   `No reference link` warnings).
6. Add any cited works to `references.bib`.
7. Update `ROADMAP.md` status and add a note to the current week's changelog under
   `docs/spec/changelog/`.
