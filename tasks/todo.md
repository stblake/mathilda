# Documentation-center function regrouping

Goal: fix how the documentation center groups builtins. Driven entirely by
`site/generate.py` (a function gets a category only if its name appears in an
`## H2` heading of `docs/spec/builtins/*.md`). Thematic headings and
undocumented builtins fall through to "Other & Advanced" (104 functions).

User decisions:
- Mechanism: curated override map in generate.py + physical split of the
  arithmetic-and-algebra spec file.
- Relational operators get a new **Comparisons** category.
- Full review: also move mis-filed functions (polynomial/algebra out of
  Structural Manipulation, etc.).

## Tasks

- [x] 1. Created `docs/spec/builtins/comparisons.md` (Comparisons category).
- [x] 2. Split `arithmetic-and-algebra.md` → `arithmetic.md` (41 secs) +
      `algebra.md` (19 secs); deleted combined file.
- [x] 2b. (added) Split mathematical constants out of `special-functions.md`
      into `mathematical-constants.md` (9 secs); Special Functions keeps Gamma,
      Pochhammer, Hypergeometric* (6 secs).
- [x] 3. generate.py: NAME_RE filter excludes internal `Head`Helper`` symbols;
      added CATEGORY_OVERRIDES (~70 entries); spec-reference link follows the
      documenting file via doc_slug.
- [x] 4. Updated `Mathilda_spec.md` index (arithmetic, algebra, comparisons,
      mathematical-constants rows; broadened pattern/special notes).
- [x] 5. Regenerated: 403 pages, 0 in "Other & Advanced". Removed stale
      `arithmetic-and-algebra/` and `other-advanced/` dirs.
- [x] 6. Changelog entry added to `docs/spec/changelog/2026-06-08.md`.
- [x] 7. Fixed live tutorial link (`05-algebra.md`) to the split categories.

## Review

Root cause: the site only categorised a function when its name appeared
literally in an `## H2` heading of a `docs/spec/builtins/*.md` file. Functions
under thematic headings (`## Trig Functions`) or undocumented fell into a
104-strong "Other & Advanced" bucket; ~25 polynomial functions were mis-filed
under Structural Manipulation.

Final tally (403 functions, 23 categories, **0 in Other**):
Arithmetic 58, Expression Information 45, Algebra 44, Linear Algebra 39,
Structural Manipulation 31, Elementary Functions 21, Functional Programming 21,
Assignment and Rules 19, Control Flow 17, Pattern Matching 15, Calculus 12,
Scoping Constructs 10, Statistics 10, Comparisons 9, Mathematical Constants 9,
File I/O 7, Simplification 7, Random 6, Special Functions 6, Lists 5, String 5,
Time/Date 4, Power Series 3.

Mechanism notes: grouping corrections live in one reviewable place
(`CATEGORY_OVERRIDES` in `site/generate.py`); the doc bodies stay where they are
documented and the per-function "Specification" link follows the documenting
file, so display category and source-of-truth can differ without broken links.
Internal context-qualified helpers (`Solve`SolveLinearSystem``, `Sum`Gosper``,
…) are no longer emitted as public pages.
