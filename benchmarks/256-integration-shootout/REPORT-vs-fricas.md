# Integration shootout: 256-integration-shootout (2026-11-21)

1266 cases, timeout 30 s per run, wall 9.9 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | R rational (72) | T transcendental (217) | A algebraic (400) | M mixed (radical on top) (408) | B radical at the bottom (169) | all |
|---|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (72) | 98% (212) | 58% (230) | 84% (343) | 45% (76) | 74% (933) |
| Mathilda:ParallelMixedTower | 0% (0) | 64% (138) | 52% (210) | 79% (323) | 40% (68) | 58% (739) |
| Mathilda:RischTranscendental | 100% (72) | 67% (145) | · | · | · | 75% (217) |
| Mathilda:Trager | · | · | 53% (211) | · | · | 53% (211) |
| Mathilda:Bronstein | · | · | · | 83% (340) | 41% (70) | 71% (410) |
| FriCAS | 99% (71) | 94% (203) | 51% (204) | 40% (163) | 43% (72) | 56% (713) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 832 | 101 | 278 | 46 | 8 | 0 | 1 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 739 | 0 | 234 | 235 | 58 | 0 | 0 | 0 | 0.07 |
| Mathilda:RischTranscendental | 210 | 7 | 67 | 5 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:Trager | 211 | 0 | 189 | 0 | 0 | 0 | 0 | 0 | 0.04 |
| Mathilda:Bronstein | 410 | 0 | 117 | 49 | 1 | 0 | 0 | 0 | 0.09 |
| FriCAS | 644 | 69 | 0 | 239 | 52 | 261 | 0 | 1 | 0.07 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| R | rherm | 36 | 36 | 36 |
| R | rlog | 36 | 36 | 35 |
| T | tlog | 33 | 33 | 33 |
| T | texp | 33 | 33 | 32 |
| T | tmixed | 15 | 15 | 14 |
| T | tnested | 13 | 13 | 13 |
| T | ttrig | 18 | 18 | 18 |
| T | tconstk | 26 | 26 | 26 |
| T | tlog-perturbed | 35 | 34 | 30 |
| T | texp-perturbed | 36 | 33 | 31 |
| A | agen | 36 | 36 | 35 |
| B | bbottom | 20 | 20 | 20 |
| B | bbottom-perturbed | 33 | 5 | 1 |
| M | alg | 34 | 34 | 3 |
| M | logxy | 31 | 31 | 0 |
| M | asinh | 35 | 35 | 34 |
| M | genus1 | 35 | 35 | 34 |
| M | algres | 35 | 35 | 34 |
| M | ratb | 30 | 30 | 19 |
| M | constk | 32 | 32 | 0 |
| M | tower2 | 35 | 35 | 3 |
| M | lccontent | 35 | 34 | 0 |
| M | expcube | 36 | 36 | 34 |
| M | mrandom | 33 | 0 | 0 |
| M | mperturb | 34 | 4 | 2 |
| B | bgenus | 24 | 24 | 24 |
| B | bcurve | 26 | 26 | 26 |
| B | bgenus-perturbed | 34 | 0 | 0 |
| B | bcurve-perturbed | 31 | 0 | 0 |
| T | classic-nonelem | 8 | 7 | 6 |
| M | classic-nonelem | 3 | 2 | 0 |
| A | classic-nonelem | 2 | 2 | 0 |
| B | classic-nonelem | 1 | 1 | 1 |
| A | trager-binomial-xk | 45 | 23 | 20 |
| A | trager-third-kind | 30 | 18 | 14 |
| A | trager-rational-times-radical | 50 | 19 | 16 |
| A | trager-chebyshev | 40 | 24 | 20 |
| A | trager-goursat | 20 | 12 | 11 |
| A | trager-pseudo-elliptic-recip | 10 | 10 | 10 |
| A | trager-even-substitution | 6 | 6 | 6 |
| A | trager-genus-hyper | 20 | 6 | 4 |
| A | trager-genus-super | 10 | 1 | 1 |
| A | trager-rational-radicand | 20 | 14 | 12 |
| A | trager-torsion-heavy | 18 | 0 | 0 |
| A | trager-deriv-alg | 18 | 18 | 18 |
| A | trager-deriv-log | 23 | 23 | 23 |
| A | trager-random | 52 | 18 | 14 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 188

## Wrong answers from other systems (their bugs)

- FriCAS `trager-random-048`: `(2 + 3*x - 3*x^2)/(Sqrt[1 + 2*x]*(2 - x + 3*x^2))`
