# Integration shootout: 260-radical-bottom (2026-11-26)

Re-tabulated from the run of 2026-11-26: 2 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

269 cases, timeout 30 s per run, wall 7.6 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | B radical at the bottom (269) | all |
|---|---|---|
| Mathilda:Automatic | 46% (124) | 46% (124) |
| Mathilda:ParallelMixedTower | 26% (71) | 26% (71) |
| Mathilda:Bronstein | 41% (111) | 41% (111) |
| FriCAS | 48% (129) | 48% (129) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 111 | 13 | 141 | 4 | 0 | 0 | 0 | 0 | 0.17 |
| Mathilda:ParallelMixedTower | 71 | 0 | 23 | 139 | 36 | 0 | 0 | 0 | 1.73 |
| Mathilda:Bronstein | 111 | 0 | 154 | 4 | 0 | 0 | 0 | 0 | 0.16 |
| FriCAS | 109 | 20 | 0 | 131 | 9 | 0 | 0 | 0 | 0.17 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 124 | 124 | 0 | 0 | 0 | 0 | 124 | 0 |
| Mathilda:ParallelMixedTower | 71 | 71 | 0 | 0 | 0 | 0 | 71 | 0 |
| Mathilda:Bronstein | 111 | 111 | 0 | 0 | 0 | 0 | 111 | 0 |
| FriCAS | 129 | 129 | 0 | 0 | 0 | 0 | 129 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| B | one-elem | 28 | 28 | 28 |
| B | one-rand | 37 | 9 | 5 |
| B | mid-elem | 29 | 29 | 29 |
| B | mid-rand | 31 | 0 | 0 |
| B | two-elem | 50 | 50 | 48 |
| B | two-rand | 49 | 0 | 0 |
| B | pert-nonelem | 45 | 8 | 19 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 2

## Wrong answers from other systems (their bugs)

