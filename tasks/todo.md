# Task: Book §4.5 "Numerical Calculus" — DONE

Plan: `~/.claude/plans/let-s-plan-writing-4-5-fizzy-lollipop.md`
Scope: Everything (all N-functions + optimization family), NDSolve figures via `\plotcell`.

## What shipped
- **Prose**: `book/chapters/math/numerical-calculus.tex` — 12 subsections (printed pp. 63–82):
  why-numeric, the shared goals contract, NIntegrate, NSum/NProduct, ND, NLimit/NSeries,
  NResidue, NDSolve, NRoots, NSolve/FindRoot, optimization (FindMinimum/FindMaximum +
  NMinimize/NMaximize), closing. Every function's Method sub-methods + options documented
  in detail (verbatim `\usagebox` card each) with classical examples and algorithm callouts.
- **Examples**: 16 `.m` files in `book/examples/numerical-calculus/` → 60 verified transcripts.
- **Figures**: `vanderpol.m` (phase-plane limit cycle) + `heat-profile.m` (diffusion), both
  headless PDF via `\plotcell`.
- **Bib**: 7 entries added to `references.bib` (QUADPACK, Takahasi–Mori, Hairer, Aberth,
  Storn–Price, Kirkpatrick, Jones/DIRECT).
- **Housekeeping**: ROADMAP §4.5 → Verified; changelog note (2026-08-24 week).

## Verification (all green)
- [x] `make examples` — 60 transcripts, every value confirmed correct
- [x] `make figures` — both figures render (viewed PDFs: limit cycle + decaying profiles)
- [x] `make usage` — 14 usage cards generated
- [x] `make check-links` — 0 unlinked `\B{}`
- [x] `make pdf` — clean log: no undefined refs/citations, no missing placeholders (110 pp.)
- [x] Concept `\index{}` entries captured by makeindex
- [x] All 11 numbered subsections in TOC; content confirmed via PDF render + pdftotext

## Review notes / decisions
- Pipeline captures results only (stderr `::accgl` messages don't appear) → transcripts clean;
  the Wynn-on-Basel *pitfall* shows as a visibly wrong value (1.62533 vs 1.64493).
- `NSeries` needs `Chop` (spec includes negative-power noise terms); `NResidue` on an essential
  singularity needs `Radius -> 1` (default 1/100 overflows). Both taught explicitly.
- Wilkinson roots come out clean at machine precision (Newton-polished) → framed as a strength.
- No `src/` changes (book-only). Doc/source discrepancies (NDSolve MaxSteps, NMinimize not
  HoldAll, FindMinimum auto-start 0) reflected code-true in prose; a separate docs-sync fix
  remains optional/out-of-scope.
