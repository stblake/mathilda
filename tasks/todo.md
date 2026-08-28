# TODO — Statistics in the Mathilda book

## Phase 0 — Review (done)
- [x] Build binary; confirm all stats builtins work
- [x] Confirm all candidate builtins are `\B{}`-linkable
- [x] Write `tasks/stats-book-review.md`

## Phase 1 — Ch 3 tour section (done)
- [x] `book/examples/03-introduction/statistics.m`
- [x] Insert `\section{Statistics}` → §3.11 (Linear algebra correctly renumbers to 3.12)

## Phase 2 — Ch 4 deep-dive section (§4.8) (done)
- [x] 8 example files under `book/examples/statistics/`
- [x] `book/chapters/math/statistics.tex` (§4.8, 8 subsections + Where this connects)
- [x] `\input{chapters/math/statistics}` in `book/chapters/04-mathematics.tex`

## Phase 3 — Build, verify, index, changelog (done)
- [x] `make examples` — all transcripts verified; prose cites printed outputs
- [x] `make check-links` — 0 unlinked *from my content* (2 pre-existing misses noted)
- [x] hand `\index{}` concept entries (Anscombe, Bessel's correction, kurtosis, …)
- [x] `make pdf` — clean build, 139 pp, no errors/undefined refs; visually spot-checked
- [x] `book/references.bib` — Anscombe (1973) added + `\cite`d (biber resolves)
- [x] changelog `docs/spec/changelog/2026-08-24.md`; ROADMAP §4.8 row + scope note

## Review / Results

**Done, document-only as agreed.** Two verified sections added:
- **§3.11 Statistics** (Ch 3 tour): exam-scores dataset — Mean/Median/StandardDeviation/
  Quartiles/Correlation/MovingAverage; forward-refs §4.8.
- **§4.8 Statistics** (Ch 4 deep dive, `chapters/math/statistics.tex`): 8 subsections on
  classical datasets — Newcomb (median robustness), Michelson (exact variance/Bessel),
  Cavendish (five-number), moments/skewness/kurtosis, **Anscombe's quartet**, moving
  statistics, frequencies, and the distribution bridge (LLN). Reference cards for
  Variance & Correlation; Theory/Pitfall/Performance/Under-the-hood callouts.

**Phase-0 review** (`tasks/stats-book-review.md`): all 15 `src/stats/` builtins + adjacent
surface exercised against the binary — nothing broken; all 37 candidates `\B{}`-linkable.

**Follow-up done — reference pages added.** `make check-links` had failed on
`ChineseRemainder` and `CylindricalDecomposition` (referenced by book prose §3.6/§4.2
but absent from `site/docs/assets/builtins.json`). Fixed by authoring
`site/overlays/<Name>.md` + `site/impl/<Name>.md` for both and regenerating via
`site/generate.py`. Reverted the regen's unrelated drift (nondeterministic
`MemoryInUse`/`Histogram`/plot-json churn, a version bump, a stale-source `Reduce`
update, and a `Prime.md` test-link mis-association) to keep a focused diff: only
`builtins.json` + 3 `index.md` (clean additions) + the 6 new files. `mkdocs build
--strict` passes (exit 0); `check-links` now green; book PDF has 0 "No reference link"
warnings.

**Build artifacts** (`book/generated/`, PDF) are git-ignored — no tracked build output.
