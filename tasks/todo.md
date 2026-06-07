# Task: Add examples to every builtin in the documentation center

## Problem
154 builtin pages in the docs site show "_No verified examples yet for this
function._" with no examples anywhere (e.g. MapAll). 15 more have overlay-only
examples (acceptable, like Map). The user wants every builtin to have a couple
of simple examples (+ a harder one where appropriate).

## Root cause
`site/generate.py` mines `In[]:=` example blocks from `docs/spec/builtins/*.md`,
but only maps a function to a spec section when the function name appears as an
H2 heading token. Functions grouped under shared headings (e.g. `## Trig
Functions` covering Sin/Cos/Tan/...) get no section -> no mined examples. Per-
function examples cannot be expressed in a grouped spec section.

## Approach
Use the established **overlay** mechanism: `site/overlays/<Name>.md`. The
generator merges overlay worked-examples into a "Notes & additional examples"
section (this is exactly how `Map` shows its examples today). One file per
function => parallel subagents never conflict. Every example is verified by
running it through `./Mathilda` and capturing the real output.

## Plan
- [ ] Build the definitive list of 154 functions missing all examples (done: .tmp_noex.txt)
- [ ] Dispatch subagents (batched by category) to create verified overlays
- [ ] Re-run detection; mop up any stragglers
- [ ] `make docs` to regenerate + verify pages render examples
- [ ] Confirm "_No verified examples_" count dropped to ~0
- [ ] Changelog entry under docs/spec/changelog/<Monday>.md
- [ ] Review section below

## Constraints
- No edits to src/ or src/external/. No edits to generate.py.
- Every overlay example must be run through ./Mathilda; outputs must be the
  real captured output (no fabricated outputs).
- Skip examples that error / return Null / are unevaluated unless that IS the
  point being illustrated.

## Review
Done. Created 154 hand-curated overlays under `site/overlays/` (one per builtin
that had no examples), each with 2-4 worked examples verified against `./Mathilda`
plus a terse note. Dispatched 12 category-batched subagents in parallel; each
verified every example by piping inputs through the binary and capturing the real
`Out[]` (no fabricated outputs). Overlay files are one-per-function so the parallel
agents never conflicted.

Results after `make docs`:
- Overlays: 60 -> 214. Verified examples reported by the generator: 1020.
- Pages with NO examples anywhere: 154 -> 0 (every public builtin now shows examples).
- Strict MkDocs build passes clean.
- No edits to src/, src/external/, or generate.py. New overlays carry no status
  front matter (no over-claiming) — status stays auto-derived.
- Changelog entry added to docs/spec/changelog/2026-06-01.md.

Known cosmetic limitation (pre-existing, matches the 60 prior overlays incl. Map):
overlay examples render under "## Notes & additional examples"; the main
"## Examples" section still reads "No verified examples yet" for these because the
generator only auto-mines the main section from per-function spec H2 headings.
Promoting overlay examples into the main section would require a generate.py
change (out of scope; left as a possible follow-up).
