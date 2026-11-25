# Integration shootout: 261-compact-integrals-x4 (2026-11-24)

Re-tabulated from the run of 2026-11-24: 1492 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

58 cases, timeout 30 s per run, wall 1.3 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | A algebraic (14) | M mixed (radical on top) (44) | all |
|---|---|---|---|
| Mathilda:Automatic | 100% (14) | 100% (44) | 100% (58) |
| Mathilda:ParallelMixedTower | 50% (7) | 57% (25) | 55% (32) |
| Mathilda:RischTranscendental | · | · | - (0) |
| Mathilda:Trager | 100% (14) | · | 100% (14) |
| Mathilda:Bronstein | · | 98% (43) | 98% (43) |
| FriCAS | 79% (11) | 7% (3) | 24% (14) |
| Axiom | 79% (11) | 7% (3) | 24% (14) |
| Maxima | 50% (7) | 5% (2) | 16% (9) |
| sympy | 21% (3) | 0% (0) | 5% (3) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 58 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.46 |
| Mathilda:ParallelMixedTower | 32 | 0 | 0 | 3 | 23 | 0 | 0 | 0 | 0.26 |
| Mathilda:RischTranscendental | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | - |
| Mathilda:Trager | 14 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.84 |
| Mathilda:Bronstein | 43 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0.45 |
| FriCAS | 14 | 0 | 0 | 0 | 10 | 34 | 0 | 0 | 0.17 |
| Axiom | 14 | 0 | 0 | 0 | 5 | 39 | 0 | 0 | 0.17 |
| Maxima | 9 | 0 | 0 | 49 | 0 | 0 | 0 | 0 | 0.04 |
| sympy | 0 | 3 | 0 | 44 | 7 | 0 | 4 | 0 | 2.14 |

## Cross-check in Mathilda (answers that failed their own check)

Own check: `bad` = a sample point failed; `na` = no point could be evaluated. Verdict of PossibleZeroQ[D[F, x] - f] in Mathilda.

| run | own | n | zero (verified) | nonzero (wrong) | undecided | untranslatable |
|---|---|---|---|---|---|---|
| Axiom | bad | 1 | 1 | 0 | 0 | 0 |
| sympy | bad | 4 | 0 | 0 | 4 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| A | cA-mathilda | 14 | 14 | 11 | 11 | 7 | 3 |
| M | cM-mathilda | 41 | 41 | 0 | 0 | 0 | 0 |
| M | cM-fricas | 3 | 3 | 3 | 3 | 2 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 44

## Wrong answers from other systems (their bugs)

