# Task: Write §4.2 "Algebra" of The Mathilda Book

Plan: `/Users/user/.claude/plans/let-s-plan-the-writing-curious-kernighan.md`
Template: `book/chapters/math/arithmetic.tex` (§4.1, Verified).

## Steps

- [x] Create `book/examples/algebra/*.m` (11 files, inputs only)
- [x] `make examples` — generate transcripts, read & adjust inputs
- [x] Write `book/chapters/math/algebra.tex` (10 subsections, callouts, \index)
- [x] Add missing citations to `book/references.bib` (added `clo`)
- [x] `make check-links` — 0 unlinked \B{} (331 uses resolve)
- [x] `make pdf` — clean log; Index + layout verified (78 pages)
- [x] Update `book/ROADMAP.md` (§4.2 → Verified) + changelog note

## Subsections (each = codepairs of verified examples + Theory/Under-the-hood callouts)

1. Polynomials & their anatomy — PolynomialQ, Variables, Exponent, Coefficient(List), Collect, HornerForm
2. Expanding & collecting — Expand, ExpandAll, ExpandNumerator/Denominator (FLINT fmpq_mpoly)
3. Factoring over Q — Factor, FactorList, FactorSquareFree, FactorTerms (Yun; Berlekamp–Zassenhaus+Hensel)
4. Algebraic number fields — Factor+Extension, IrreduciblePolynomialQ, MinimalPolynomial, RootReduce, Root, ToRadicals (Trager)
5. GCD, resultants, elimination — Polynomial* family, Resultant/Discriminant/Subresultants, Eliminate (subresultant PRS)
6. Rational functions — Numerator/Denominator, Together, Cancel, Apart
7. Gröbner bases — GroebnerBasis, MonomialOrder, Modulus (Buchberger + Gröbner walk)
8. Solving polynomial equations — Solve single: quadratic→Cardano(Cubics)→Root/Abel–Ruffini
9. Solving polynomial systems — Solve systems: lex GB → triangular decomposition
10. Reduce & the reals — Solve-vs-Reduce, Domain, McCallum CAD

## Do NOT \B{} (not implemented): PolynomialReduce, CoefficientRules, MonomialList, FindInstance, Roots, AlgebraicNumber

## Review

**Done.** §4.2 Algebra written (`book/chapters/math/algebra.tex`, ~10 subsections,
printed pages 42–50), replacing the stub. Documentation only — no C source changes.

- 11 example files under `book/examples/algebra/`; every `Out[]` build-verified by
  `make examples`. Two long outputs (Cardano radicals; the six-permutation system
  `{x+y+z==6, xy+yz+zx==11, xyz==6}` with its triangular Gröbner basis) use
  full-width `\mtranscript`; the rest use `codepairs`.
- Callouts carry the algorithms per ROADMAP: Yun square-free, Berlekamp–Zassenhaus +
  Hensel, Trager, subresultant PRS, Buchberger + Gröbner walk, Cardano/Abel–Ruffini,
  McCallum CAD. FLINT fast-path and the packed `fmpq_mpoly`/`Overflow[]` guard noted.
- Avoided §3 duplication (went deeper on extensions, GCD/resultants, Gröbner orders,
  Root-vs-radical, systems, and Reduce/CAD). `Simplify` deferred per user decision.
- Guards NOT tripped: did not `\B{}` the non-builtins (`PolynomialReduce`,
  `CoefficientRules`, `MonomialList`, `FindInstance`, `Roots`, `AlgebraicNumber`);
  used `Extension -> I` for the Gaussian example after confirming
  `GaussianIntegers -> True` does not factor `x^2+1`.
- `make check-links` clean; `make pdf` clean (no undefined refs/citations/links);
  Index concept entries verified in `.idx`. ROADMAP → Verified; changelog note added.
