# Integration shootout: 261-compact-integrals-x4 (2026-11-26)

Re-tabulated from the run of 2026-11-26: 3 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

43 cases, timeout 30 s per run, wall 0.0 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | A algebraic (8) | M mixed (radical on top) (35) | all |
|---|---|---|---|
| Mathilda:Automatic | 100% (8) | 100% (35) | 100% (43) |
| Mathilda:ParallelMixedTower | 88% (7) | 63% (22) | 67% (29) |
| Mathilda:Trager | 100% (8) | · | 100% (8) |
| Mathilda:Bronstein | · | 97% (34) | 97% (34) |
| FriCAS | 88% (7) | 3% (1) | 19% (8) |
| Axiom | 88% (7) | 3% (1) | 19% (8) |
| Maxima | 88% (7) | 3% (1) | 19% (8) |
| sympy | 88% (7) | 0% (0) | 16% (7) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 43 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.18 |
| Mathilda:ParallelMixedTower | 29 | 0 | 0 | 3 | 11 | 0 | 0 | 0 | 0.11 |
| Mathilda:Trager | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.28 |
| Mathilda:Bronstein | 34 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0.20 |
| FriCAS | 8 | 0 | 0 | 0 | 2 | 33 | 0 | 0 | 0.17 |
| Axiom | 8 | 0 | 0 | 0 | 2 | 33 | 0 | 0 | 0.20 |
| Maxima | 8 | 0 | 0 | 35 | 0 | 0 | 0 | 0 | 0.04 |
| sympy | 0 | 7 | 0 | 33 | 3 | 0 | 0 | 0 | 0.99 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 43 | 43 | 0 | 0 | 0 | 0 | 43 | 0 |
| Mathilda:ParallelMixedTower | 29 | 29 | 0 | 0 | 0 | 0 | 29 | 0 |
| Mathilda:Trager | 8 | 8 | 0 | 0 | 0 | 0 | 8 | 0 |
| Mathilda:Bronstein | 34 | 34 | 0 | 0 | 0 | 0 | 34 | 0 |
| FriCAS | 8 | 8 | 0 | 0 | 0 | 0 | 8 | 0 |
| Axiom | 8 | 8 | 0 | 0 | 0 | 0 | 8 | 0 |
| Maxima | 8 | 8 | 0 | 0 | 0 | 0 | 8 | 0 |
| sympy | 7 | 7 | 0 | 0 | 0 | 0 | 3 | 4 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| A | cA-mathilda | 8 | 8 | 7 | 7 | 7 | 7 |
| M | cM-mathilda | 34 | 34 | 0 | 0 | 0 | 0 |
| M | cM-fricas | 1 | 1 | 1 | 1 | 1 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 35

## Wrong answers from other systems (their bugs)

