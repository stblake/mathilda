# Integration shootout: 261-compact-integrals-x4 (2026-11-24)

Re-tabulated from the run of 2026-11-24: 1492 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

58 cases, timeout 30 s per run, wall 1.3 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

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
| sympy | 50% (7) | 0% (0) | 12% (7) |

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
| sympy | 0 | 7 | 0 | 44 | 7 | 0 | 0 | 0 | 0.99 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 58 | 58 | 0 | 0 | 0 | 0 | 58 | 0 |
| Mathilda:ParallelMixedTower | 32 | 32 | 0 | 0 | 0 | 0 | 32 | 0 |
| Mathilda:Trager | 14 | 14 | 0 | 0 | 0 | 0 | 14 | 0 |
| Mathilda:Bronstein | 43 | 43 | 0 | 0 | 0 | 0 | 43 | 0 |
| FriCAS | 14 | 14 | 0 | 0 | 0 | 0 | 14 | 0 |
| Axiom | 14 | 14 | 0 | 0 | 0 | 0 | 13 | 1 |
| Maxima | 9 | 9 | 0 | 0 | 0 | 0 | 9 | 0 |
| sympy | 7 | 7 | 0 | 0 | 0 | 0 | 3 | 4 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| A | cA-mathilda | 14 | 14 | 11 | 11 | 7 | 7 |
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

