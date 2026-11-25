# Integration shootout: 264-mixed-integrals (2026-11-26)

500 cases, timeout 30 s per run, wall 49.7 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | M mixed (radical on top) (500) | all |
|---|---|---|
| Mathilda:Automatic | 100% (500) | 100% (500) |
| Mathilda:ParallelMixedTower | 0% (0) | 0% (0) |
| Mathilda:RischTranscendental | · | - (0) |
| Mathilda:Trager | · | - (0) |
| Mathilda:Bronstein | 100% (500) | 100% (500) |
| FriCAS | 0% (0) | 0% (0) |
| Axiom | 0% (0) | 0% (0) |
| Maxima | 0% (0) | 0% (0) |
| sympy | 0% (0) | 0% (0) |
| Wolfram | 0% (1) | 0% (1) |
| Reduce | 0% (0) | 0% (0) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 500 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.03 |
| Mathilda:ParallelMixedTower | 0 | 0 | 0 | 495 | 5 | 0 | 0 | 0 | - |
| Mathilda:RischTranscendental | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | - |
| Mathilda:Trager | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | - |
| Mathilda:Bronstein | 500 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.04 |
| FriCAS | 0 | 0 | 0 | 0 | 0 | 500 | 0 | 0 | - |
| Axiom | 0 | 0 | 0 | 0 | 0 | 500 | 0 | 0 | - |
| Maxima | 0 | 0 | 0 | 500 | 0 | 0 | 0 | 0 | - |
| sympy | 0 | 0 | 0 | 299 | 201 | 0 | 0 | 0 | - |
| Wolfram | 1 | 0 | 0 | 499 | 0 | 0 | 0 | 0 | 0.46 |
| Reduce | 0 | 0 | 0 | 446 | 54 | 0 | 0 | 0 | - |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 500 | 500 | 0 | 0 | 0 | 0 | 500 | 0 |
| Mathilda:Bronstein | 500 | 500 | 0 | 0 | 0 | 0 | 500 | 0 |
| Wolfram | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy | Wolfram | Reduce |
|---|---|---|---|---|---|---|---|---|---|
| M | cM-mathilda | 500 | 500 | 0 | 0 | 0 | 0 | 1 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 499

## Wrong answers from other systems (their bugs)

