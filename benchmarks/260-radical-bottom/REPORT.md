# Integration shootout: 260-radical-bottom (2026-11-22)

Re-tabulated from the run of 2026-11-22: 149 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

271 cases, timeout 30 s per run, wall 8.6 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | B radical at the bottom (271) | all |
|---|---|---|
| Mathilda:Automatic | 45% (121) | 45% (121) |
| Mathilda:ParallelMixedTower | 27% (73) | 27% (73) |
| Mathilda:Bronstein | 41% (110) | 41% (110) |
| FriCAS | 48% (131) | 48% (131) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 108 | 13 | 143 | 7 | 0 | 0 | 0 | 0 | 0.20 |
| Mathilda:ParallelMixedTower | 73 | 0 | 23 | 143 | 32 | 0 | 0 | 0 | 0.74 |
| Mathilda:Bronstein | 110 | 0 | 156 | 5 | 0 | 0 | 0 | 0 | 0.08 |
| FriCAS | 110 | 21 | 0 | 132 | 7 | 1 | 0 | 0 | 0.13 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 121 | 121 | 0 | 0 | 0 | 0 | 121 | 0 |
| Mathilda:ParallelMixedTower | 73 | 73 | 0 | 0 | 0 | 0 | 73 | 0 |
| Mathilda:Bronstein | 110 | 110 | 0 | 0 | 0 | 0 | 110 | 0 |
| FriCAS | 131 | 131 | 0 | 0 | 0 | 0 | 131 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| B | one-elem | 28 | 28 | 28 |
| B | one-rand | 37 | 9 | 5 |
| B | mid-elem | 29 | 29 | 29 |
| B | mid-rand | 31 | 0 | 0 |
| B | two-elem | 50 | 49 | 49 |
| B | two-rand | 49 | 0 | 0 |
| B | pert-nonelem | 47 | 8 | 20 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (2)

- `two-elem-6` AUTO-GAP Mathilda:Bronstein: `(-3*E^(x + x^2)*x)/(2 + x^2)^(3/2) + (3*E^(x + x^2)*(1 + 2*x))/Sqrt[2 + x^2] - (Sqrt[2 + x^2]*(1 + x^3)^(1/3))/(1 + x) + (1 + x^2/(1 + x^3)^(2/3))/(2*(x + (1 + x^3)^(1/3))) - (x^2*Sqrt[2 + x^2]*Log[1 + x])/(1 + x^3)^(2/3`
- `two-elem-49` AUTO-GAP Mathilda:Bronstein: `(E^(x + x^2)*x^3)/(1 + x^3)^(1/3) + (E^(x + x^2)*(1 + x^3)^(2/3))/2 + (E^(x + x^2)*x*(1 + 2*x)*(1 + x^3)^(2/3))/2 + (E^x^2*x^3*(1 + x^3)^(1/3))/Sqrt[3 + x^4] + E^(2*x)*x*Sqrt[3 + x^4] + (E^x^2*x^2*Sqrt[3 + x^4])/(2*(1 + `

## Only one system closes it

- Mathilda: 8
- FriCAS: 1

## Wrong answers from other systems (their bugs)

