# Integration shootout: 257-compact-integrals (2026-11-22)

Re-tabulated from the run of 2026-11-22: 982 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

18 cases, timeout 30 s per run, wall 6.3 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | T transcendental (2) | A algebraic (11) | M mixed (radical on top) (5) | all |
|---|---|---|---|---|
| Mathilda:Automatic | 100% (2) | 100% (11) | 100% (5) | 100% (18) |
| Mathilda:ParallelMixedTower | 0% (0) | 82% (9) | 60% (3) | 67% (12) |
| Mathilda:RischTranscendental | 100% (2) | · | · | 100% (2) |
| Mathilda:Trager | · | 100% (11) | · | 100% (11) |
| Mathilda:Bronstein | · | · | 80% (4) | 80% (4) |
| FriCAS | 100% (2) | 100% (11) | 20% (1) | 78% (14) |
| Maxima | 100% (2) | 82% (9) | 0% (0) | 61% (11) |
| sympy | 0% (0) | 64% (7) | 0% (0) | 39% (7) |
| Wolfram | 100% (2) | 100% (11) | 80% (4) | 94% (17) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 18 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.03 |
| Mathilda:ParallelMixedTower | 12 | 0 | 0 | 4 | 2 | 0 | 0 | 0 | 0.03 |
| Mathilda:RischTranscendental | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.08 |
| Mathilda:Trager | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.32 |
| Mathilda:Bronstein | 4 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0.11 |
| FriCAS | 14 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0.09 |
| Maxima | 11 | 0 | 0 | 7 | 0 | 0 | 0 | 0 | 0.02 |
| sympy | 0 | 7 | 0 | 5 | 4 | 0 | 1 | 1 | 0.73 |
| Wolfram | 16 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0.30 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 2 | 2 | 2 | 2 | 0 | 2 |
| A | cA-mathilda | 11 | 11 | 11 | 9 | 7 | 11 |
| M | cM-mathilda | 4 | 4 | 0 | 0 | 0 | 4 |
| M | cM-fricas | 1 | 1 | 1 | 0 | 0 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 1

## Wrong answers from other systems (their bugs)

- sympy `cA-3-52`: `x^3*(-1 + x^3)^(2/3)`
