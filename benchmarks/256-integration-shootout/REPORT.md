# Integration shootout: 256-integration-shootout (2026-11-26)

Re-tabulated from the run of 2026-11-26: 22 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

699 cases, timeout 30 s per run, wall 16.3 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (125) | A algebraic (216) | M mixed (radical on top) (256) | B radical at the bottom (102) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 96% (120) | 25% (54) | 74% (189) | 25% (25) | 56% (388) |
| Mathilda:ParallelMixedTower | 39% (49) | 19% (41) | 69% (176) | 24% (24) | 41% (290) |
| Mathilda:RischTranscendental | 45% (56) | · | · | · | 45% (56) |
| Mathilda:Trager | · | 19% (41) | · | · | 19% (41) |
| Mathilda:Bronstein | · | · | 73% (186) | 25% (25) | 59% (211) |
| FriCAS | 91% (114) | 18% (39) | 11% (29) | 25% (25) | 30% (207) |
| Maxima | 54% (68) | 6% (13) | 28% (71) | 18% (18) | 24% (170) |
| sympy | 67% (84) | 27% (59) | 2% (6) | 6% (6) | 22% (155) |
| Wolfram | 99% (124) | 60% (130) | 40% (103) | 29% (30) | 55% (387) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 302 | 86 | 308 | 1 | 0 | 0 | 2 | 0 | 0.08 |
| Mathilda:ParallelMixedTower | 290 | 0 | 216 | 136 | 57 | 0 | 0 | 0 | 0.08 |
| Mathilda:RischTranscendental | 49 | 7 | 69 | 0 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:Trager | 41 | 0 | 175 | 0 | 0 | 0 | 0 | 0 | 0.10 |
| Mathilda:Bronstein | 211 | 0 | 143 | 1 | 1 | 0 | 2 | 0 | 0.08 |
| FriCAS | 142 | 65 | 0 | 214 | 51 | 226 | 0 | 1 | 0.07 |
| Maxima | 142 | 28 | 0 | 506 | 5 | 0 | 12 | 6 | 0.07 |
| sympy | 60 | 95 | 0 | 348 | 184 | 1 | 11 | 0 | 0.27 |
| Wolfram | 207 | 180 | 0 | 243 | 16 | 0 | 50 | 3 | 0.35 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 390 | 388 | 0 | 2 | 0 | 0 | 388 | 0 |
| Mathilda:ParallelMixedTower | 290 | 290 | 0 | 0 | 0 | 0 | 290 | 0 |
| Mathilda:RischTranscendental | 56 | 56 | 0 | 0 | 0 | 0 | 56 | 0 |
| Mathilda:Trager | 41 | 41 | 0 | 0 | 0 | 0 | 41 | 0 |
| Mathilda:Bronstein | 213 | 211 | 0 | 2 | 0 | 0 | 211 | 0 |
| FriCAS | 208 | 207 | 1 | 0 | 0 | 0 | 207 | 1 |
| Maxima | 188 | 170 | 6 | 2 | 10 | 0 | 170 | 5 |
| sympy | 166 | 155 | 0 | 1 | 10 | 0 | 130 | 25 |
| Wolfram | 440 | 387 | 3 | 47 | 0 | 3 | 381 | 6 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| T | tlog | 15 | 15 | 15 | 15 | 15 | 15 |
| T | texp | 17 | 17 | 17 | 11 | 17 | 17 |
| T | tmixed | 4 | 4 | 4 | 4 | 4 | 4 |
| T | tnested | 2 | 2 | 2 | 1 | 2 | 2 |
| T | ttrig | 3 | 3 | 3 | 3 | 3 | 3 |
| T | tconstk | 8 | 8 | 8 | 8 | 8 | 8 |
| T | tlog-perturbed | 34 | 33 | 30 | 13 | 16 | 34 |
| T | texp-perturbed | 34 | 31 | 29 | 8 | 14 | 34 |
| A | agen | 3 | 3 | 3 | 2 | 0 | 3 |
| B | bbottom-perturbed | 21 | 0 | 0 | 0 | 0 | 7 |
| M | alg | 34 | 34 | 3 | 8 | 0 | 33 |
| M | logxy | 31 | 31 | 0 | 12 | 0 | 0 |
| M | asinh | 5 | 5 | 5 | 5 | 0 | 5 |
| M | ratb | 22 | 22 | 11 | 9 | 0 | 22 |
| M | constk | 32 | 32 | 0 | 0 | 0 | 0 |
| M | tower2 | 27 | 27 | 3 | 27 | 0 | 19 |
| M | lccontent | 29 | 28 | 0 | 4 | 0 | 13 |
| M | expcube | 4 | 4 | 4 | 4 | 4 | 4 |
| M | mrandom | 33 | 0 | 0 | 0 | 0 | 0 |
| M | mperturb | 32 | 2 | 0 | 0 | 0 | 1 |
| B | bgenus | 15 | 15 | 15 | 15 | 0 | 13 |
| B | bcurve | 10 | 10 | 10 | 3 | 6 | 10 |
| B | bgenus-perturbed | 31 | 0 | 0 | 0 | 0 | 0 |
| B | bcurve-perturbed | 25 | 0 | 0 | 0 | 0 | 0 |
| T | classic-nonelem | 8 | 7 | 6 | 5 | 5 | 7 |
| M | classic-nonelem | 3 | 2 | 0 | 2 | 2 | 2 |
| A | classic-nonelem | 2 | 2 | 0 | 0 | 2 | 2 |
| A | trager-binomial-xk | 21 | 5 | 4 | 4 | 19 | 21 |
| A | trager-third-kind | 16 | 4 | 0 | 0 | 0 | 15 |
| A | trager-rational-times-radical | 33 | 2 | 0 | 0 | 5 | 15 |
| A | trager-chebyshev | 22 | 5 | 2 | 2 | 22 | 22 |
| A | trager-goursat | 15 | 7 | 6 | 2 | 3 | 10 |
| A | trager-pseudo-elliptic-recip | 4 | 4 | 4 | 0 | 0 | 4 |
| A | trager-genus-hyper | 12 | 0 | 0 | 0 | 0 | 0 |
| A | trager-genus-super | 9 | 0 | 0 | 0 | 5 | 5 |
| A | trager-rational-radicand | 8 | 2 | 0 | 0 | 0 | 6 |
| A | trager-torsion-heavy | 18 | 0 | 0 | 0 | 2 | 3 |
| A | trager-deriv-alg | 10 | 10 | 10 | 3 | 1 | 10 |
| A | trager-deriv-log | 10 | 10 | 10 | 0 | 0 | 8 |
| A | trager-random | 33 | 0 | 0 | 0 | 0 | 6 |
| M | algres | 3 | 3 | 2 | 0 | 0 | 3 |
| M | genus1 | 1 | 1 | 1 | 0 | 0 | 1 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 68

## Wrong answers from other systems (their bugs)

- Maxima `lccontent-1-7`: `(4*x*Log[x] + x*Log[x]^2 + 2*x*(3 + Log[x]))/(x*Sqrt[x*(3 + Log[x])])`
- Wolfram `bcurvep-1-10`: `(6*E^(x*Sqrt[3 + x^3])*x + 5*E^(x*Sqrt[3 + x^3])*x^4 - 2*E^(x*Sqrt[3 + x^3])*Sqrt[3 + x^3] + 2*x^2*Sqrt[3 + x^3]*Log[x + Sqrt[1 + x^3]])/(x^2*Sqrt[3 + x^3])`
- Maxima `trager-deriv-alg-03`: `(-11 - 11*x + 5*x^2 + 5*x^3)/(2 - x - 2*x^2 + x^3)^(7/4)`
- Maxima `trager-deriv-alg-04`: `(3 + 4*x + x^2)/(-3 - 2*x + x^2)^(7/4)`
- Maxima `trager-deriv-alg-11`: `(-13 - 22*x + 9*x^2 + 12*x^3 - x^4)/(3 - 2*x - 2*x^2 + x^3)^(8/5)`
- Maxima `logxy-2-10`: `(2 - 2*E^x - 3*(2 + 2*E^x - 2*x)^(2/3))/(2 + 2*E^x + (2 + 2*E^x - 2*x)^(2/3) - 2*x + (2 + 2*E^x - 2*x)^(2/3)*x)`
- Wolfram `bgenusp-2-11`: `(-3*E^(2*x) + 6*E^(2*x)*x - E^(2*x)*x^2 + 6*E^(2*x)*x^3 + 3*x^2*(1 + x^2)^(2/3)*Log[x + Sqrt[1 + x^3]])/(x^2*(1 + x^2)^(2/3))`
- Wolfram `bcurvep-2-3`: `(2*E^(x*Sqrt[1 + x^3])*x + 5*E^(x*Sqrt[1 + x^3])*x^4 - 2*E^(x*Sqrt[1 + x^3])*Sqrt[1 + x^3] + 2*x^2*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(x^2*Sqrt[1 + x^3])`
- FriCAS `algres-3-9`: `(2*E^x^2*x)/Sqrt[2 + E^(2*x^2)] - (2*Sqrt[2 + E^(2*x^2)]*x)/E^x^2 + ((-1 + E^x^2 + Sqrt[2 + E^(2*x^2)])*((2*E^(2*x^2)*x + 2*E^x^2*Sqrt[2 + E^(2*x^2)]*x)/(2 + E^`
- Maxima `lccontent-3-8`: `(3*x*Log[x] + x*Log[x]^2 + 2*x*(2 + Log[x]))/(x*Sqrt[x*(2 + Log[x])])`
