# Integration shootout: 257-compact-integrals (2026-11-26)

Re-tabulated from the run of 2026-11-26: 2 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

11 cases, timeout 30 s per run, wall 0.0 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | A algebraic (6) | M mixed (radical on top) (5) | all |
|---|---|---|---|
| Mathilda:Automatic | 100% (6) | 100% (5) | 100% (11) |
| Mathilda:ParallelMixedTower | 100% (6) | 60% (3) | 82% (9) |
| Mathilda:RischTranscendental | · | · | - (0) |
| Mathilda:Trager | 100% (6) | · | 100% (6) |
| Mathilda:Bronstein | · | 80% (4) | 80% (4) |
| FriCAS | 100% (6) | 20% (1) | 64% (7) |
| Maxima | 100% (6) | 0% (0) | 55% (6) |
| sympy | 100% (6) | 0% (0) | 55% (6) |
| Wolfram | 100% (6) | 80% (4) | 91% (10) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:ParallelMixedTower | 9 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0.03 |
| Mathilda:RischTranscendental | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | - |
| Mathilda:Trager | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.27 |
| Mathilda:Bronstein | 4 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0.09 |
| FriCAS | 7 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0.09 |
| Maxima | 6 | 0 | 0 | 5 | 0 | 0 | 0 | 0 | 0.01 |
| sympy | 0 | 6 | 0 | 4 | 1 | 0 | 0 | 0 | 0.71 |
| Wolfram | 9 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0.31 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 11 | 11 | 0 | 0 | 0 | 0 | 11 | 0 |
| Mathilda:ParallelMixedTower | 9 | 9 | 0 | 0 | 0 | 0 | 9 | 0 |
| Mathilda:Trager | 6 | 6 | 0 | 0 | 0 | 0 | 6 | 0 |
| Mathilda:Bronstein | 4 | 4 | 0 | 0 | 0 | 0 | 4 | 0 |
| FriCAS | 7 | 7 | 0 | 0 | 0 | 0 | 7 | 0 |
| Maxima | 6 | 6 | 0 | 0 | 0 | 0 | 6 | 0 |
| sympy | 6 | 6 | 0 | 0 | 0 | 0 | 5 | 1 |
| Wolfram | 10 | 10 | 0 | 0 | 0 | 0 | 10 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| A | cA-mathilda | 6 | 6 | 6 | 6 | 6 | 6 |
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

