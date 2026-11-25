# Integration shootout: 256-integration-shootout (2026-11-22)

Re-tabulated from the run of 2026-11-22: 489 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

777 cases, timeout 30 s per run, wall 12.1 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (152) | A algebraic (223) | M mixed (radical on top) (266) | B radical at the bottom (136) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 97% (147) | 26% (58) | 68% (180) | 32% (44) | 55% (429) |
| Mathilda:ParallelMixedTower | 48% (73) | 19% (43) | 70% (186) | 31% (42) | 44% (344) |
| Mathilda:RischTranscendental | 53% (80) | · | · | · | 53% (80) |
| Mathilda:Trager | · | 19% (43) | · | · | 19% (43) |
| Mathilda:Bronstein | · | · | 67% (177) | 32% (44) | 55% (221) |
| FriCAS | 93% (141) | 18% (41) | 12% (31) | 32% (44) | 33% (257) |
| Maxima | 59% (90) | 6% (14) | 30% (81) | 22% (30) | 28% (215) |
| sympy | 72% (109) | 28% (62) | 2% (6) | 8% (11) | 24% (188) |
| Wolfram | 99% (151) | 65% (145) | 41% (110) | 36% (49) | 59% (455) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 338 | 91 | 320 | 4 | 2 | 0 | 22 | 0 | 0.06 |
| Mathilda:ParallelMixedTower | 344 | 0 | 224 | 152 | 57 | 0 | 0 | 0 | 0.06 |
| Mathilda:RischTranscendental | 73 | 7 | 72 | 0 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:Trager | 43 | 0 | 180 | 0 | 0 | 0 | 0 | 0 | 0.11 |
| Mathilda:Bronstein | 221 | 0 | 158 | 1 | 1 | 0 | 21 | 0 | 0.07 |
| FriCAS | 189 | 68 | 0 | 233 | 52 | 235 | 0 | 0 | 0.07 |
| Maxima | 186 | 29 | 0 | 539 | 5 | 0 | 12 | 6 | 0.07 |
| sympy | 89 | 99 | 0 | 373 | 205 | 1 | 10 | 0 | 0.25 |
| Wolfram | 256 | 199 | 0 | 262 | 13 | 0 | 43 | 4 | 0.32 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 451 | 429 | 0 | 0 | 0 | 22 | 429 | 0 |
| Mathilda:ParallelMixedTower | 344 | 344 | 0 | 0 | 0 | 0 | 344 | 0 |
| Mathilda:RischTranscendental | 80 | 80 | 0 | 0 | 0 | 0 | 80 | 0 |
| Mathilda:Trager | 43 | 43 | 0 | 0 | 0 | 0 | 43 | 0 |
| Mathilda:Bronstein | 242 | 221 | 0 | 0 | 0 | 21 | 221 | 0 |
| FriCAS | 257 | 257 | 0 | 0 | 0 | 0 | 257 | 0 |
| Maxima | 233 | 215 | 6 | 2 | 10 | 0 | 215 | 5 |
| sympy | 198 | 188 | 0 | 0 | 10 | 0 | 162 | 26 |
| Wolfram | 502 | 455 | 4 | 35 | 0 | 8 | 444 | 11 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| T | tlog | 19 | 19 | 19 | 19 | 19 | 19 |
| T | texp | 20 | 20 | 20 | 13 | 20 | 20 |
| T | tmixed | 10 | 10 | 10 | 10 | 10 | 10 |
| T | tnested | 5 | 5 | 5 | 2 | 5 | 5 |
| T | ttrig | 7 | 7 | 7 | 7 | 7 | 7 |
| T | tconstk | 12 | 12 | 12 | 12 | 12 | 12 |
| T | tlog-perturbed | 35 | 34 | 31 | 13 | 16 | 35 |
| T | texp-perturbed | 36 | 33 | 31 | 9 | 15 | 36 |
| A | agen | 3 | 3 | 3 | 2 | 0 | 3 |
| B | bbottom-perturbed | 27 | 0 | 0 | 0 | 0 | 8 |
| M | alg | 34 | 34 | 3 | 8 | 0 | 33 |
| M | logxy | 31 | 31 | 0 | 12 | 0 | 0 |
| M | asinh | 6 | 6 | 6 | 6 | 0 | 6 |
| M | ratb | 22 | 22 | 11 | 9 | 0 | 22 |
| M | constk | 32 | 32 | 0 | 0 | 0 | 0 |
| M | tower2 | 35 | 35 | 3 | 35 | 0 | 24 |
| M | lccontent | 30 | 29 | 0 | 5 | 0 | 14 |
| M | expcube | 4 | 4 | 4 | 4 | 4 | 4 |
| M | mrandom | 33 | 0 | 0 | 0 | 0 | 0 |
| M | mperturb | 32 | 2 | 0 | 0 | 0 | 1 |
| B | bgenus | 24 | 24 | 24 | 24 | 0 | 21 |
| B | bcurve | 20 | 20 | 20 | 6 | 11 | 20 |
| B | bgenus-perturbed | 34 | 0 | 0 | 0 | 0 | 0 |
| B | bcurve-perturbed | 31 | 0 | 0 | 0 | 0 | 0 |
| T | classic-nonelem | 8 | 7 | 6 | 5 | 5 | 7 |
| M | classic-nonelem | 3 | 2 | 0 | 2 | 2 | 2 |
| A | classic-nonelem | 2 | 2 | 0 | 0 | 2 | 2 |
| A | trager-binomial-xk | 25 | 5 | 4 | 4 | 24 | 25 |
| A | trager-third-kind | 16 | 4 | 0 | 0 | 0 | 15 |
| A | trager-rational-times-radical | 33 | 2 | 0 | 0 | 5 | 21 |
| A | trager-chebyshev | 21 | 6 | 2 | 2 | 21 | 21 |
| A | trager-goursat | 15 | 7 | 6 | 2 | 3 | 11 |
| A | trager-pseudo-elliptic-recip | 4 | 4 | 4 | 0 | 0 | 4 |
| A | trager-genus-hyper | 12 | 0 | 0 | 0 | 0 | 0 |
| A | trager-genus-super | 8 | 0 | 0 | 0 | 4 | 4 |
| A | trager-rational-radicand | 8 | 2 | 0 | 0 | 0 | 7 |
| A | trager-torsion-heavy | 18 | 0 | 0 | 0 | 2 | 5 |
| A | trager-deriv-alg | 11 | 11 | 11 | 4 | 1 | 11 |
| A | trager-deriv-log | 11 | 11 | 11 | 0 | 0 | 8 |
| A | trager-random | 36 | 1 | 0 | 0 | 0 | 8 |
| M | algres | 3 | 3 | 3 | 0 | 0 | 3 |
| M | genus1 | 1 | 1 | 1 | 0 | 0 | 1 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 67

## Wrong answers from other systems (their bugs)

- Maxima `lccontent-1-7`: `(4*x*Log[x] + x*Log[x]^2 + 2*x*(3 + Log[x]))/(2*x*Sqrt[x*(3 + Log[x])])`
- Wolfram `bcurvep-1-10`: `(6*E^(x*Sqrt[3 + x^3])*x + 5*E^(x*Sqrt[3 + x^3])*x^4 - 2*E^(x*Sqrt[3 + x^3])*Sqrt[3 + x^3] + 2*x^2*Sqrt[3 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[3 + x^3])`
- Maxima `trager-deriv-alg-03`: `(-22 - 22*x + 10*x^2 + 10*x^3)/(4*(2 - x - 2*x^2 + x^3)^(7/4))`
- Maxima `trager-deriv-alg-04`: `(12 + 16*x + 4*x^2)/(4*(-3 - 2*x + x^2)^(7/4))`
- Maxima `trager-deriv-alg-11`: `(-39 - 66*x + 27*x^2 + 36*x^3 - 3*x^4)/(5*(3 - 2*x - 2*x^2 + x^3)^(8/5))`
- Maxima `logxy-2-10`: `(2 - 2*E^x - 3*(2 + 2*E^x - 2*x)^(2/3))/(3*(2 + 2*E^x - 2*x)^(2/3) + 3*(2 + 2*E^x - 2*x) + 3*(2 + 2*E^x - 2*x)^(2/3)*x)`
- Wolfram `bgenusp-2-11`: `(-3*E^(2*x) + 6*E^(2*x)*x - E^(2*x)*x^2 + 6*E^(2*x)*x^3 + 3*x^2*(1 + x^2)^(2/3)*Log[x + Sqrt[1 + x^3]])/(3*x^2*(1 + x^2)^(2/3))`
- Wolfram `bcurvep-2-3`: `(2*E^(x*Sqrt[1 + x^3])*x + 5*E^(x*Sqrt[1 + x^3])*x^4 - 2*E^(x*Sqrt[1 + x^3])*Sqrt[1 + x^3] + 2*x^2*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[1 + x^3])`
- Maxima `lccontent-3-8`: `(3*x*Log[x] + x*Log[x]^2 + 2*x*(2 + Log[x]))/(2*x*Sqrt[x*(2 + Log[x])])`
- Wolfram `bcurvep-3-6`: `(2*E^(x*Sqrt[1 + x^3])*x + 5*E^(x*Sqrt[1 + x^3])*x^4 - 2*E^(x*Sqrt[1 + x^3])*Sqrt[1 + x^3] + 2*x^2*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[1 + x^3])`
