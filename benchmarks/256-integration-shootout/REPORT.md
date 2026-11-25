# Integration shootout: 256-integration-shootout (2026-11-22)

Re-tabulated from the run of 2026-11-22: 489 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

777 cases, timeout 30 s per run, wall 12.1 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | T transcendental (152) | A algebraic (223) | M mixed (radical on top) (266) | B radical at the bottom (136) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 97% (147) | 26% (59) | 76% (201) | 32% (44) | 58% (451) |
| Mathilda:ParallelMixedTower | 48% (73) | 19% (43) | 70% (186) | 31% (42) | 44% (344) |
| Mathilda:RischTranscendental | 53% (80) | · | · | · | 53% (80) |
| Mathilda:Trager | · | 19% (43) | · | · | 19% (43) |
| Mathilda:Bronstein | · | · | 74% (198) | 32% (44) | 60% (242) |
| FriCAS | 91% (139) | 18% (41) | 12% (31) | 32% (44) | 33% (255) |
| Maxima | 67% (102) | 7% (16) | 32% (84) | 22% (30) | 30% (232) |
| sympy | 72% (109) | 17% (37) | 2% (6) | 8% (11) | 21% (163) |
| Wolfram | 99% (151) | 71% (158) | 41% (110) | 56% (76) | 64% (495) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 359 | 92 | 320 | 4 | 2 | 0 | 0 | 0 | 0.07 |
| Mathilda:ParallelMixedTower | 344 | 0 | 224 | 152 | 57 | 0 | 0 | 0 | 0.06 |
| Mathilda:RischTranscendental | 73 | 7 | 72 | 0 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:Trager | 43 | 0 | 180 | 0 | 0 | 0 | 0 | 0 | 0.11 |
| Mathilda:Bronstein | 242 | 0 | 158 | 1 | 1 | 0 | 0 | 0 | 0.08 |
| FriCAS | 188 | 67 | 0 | 226 | 52 | 244 | 0 | 0 | 0.07 |
| Maxima | 191 | 41 | 0 | 539 | 5 | 0 | 0 | 1 | 0.07 |
| sympy | 89 | 74 | 0 | 373 | 205 | 1 | 6 | 29 | 0.22 |
| Wolfram | 256 | 239 | 0 | 262 | 13 | 0 | 1 | 6 | 0.34 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| T | tlog | 19 | 19 | 19 | 19 | 19 | 19 |
| T | texp | 20 | 20 | 20 | 13 | 20 | 20 |
| T | tmixed | 10 | 10 | 9 | 10 | 10 | 10 |
| T | tnested | 5 | 5 | 5 | 2 | 5 | 5 |
| T | ttrig | 7 | 7 | 7 | 7 | 7 | 7 |
| T | tconstk | 12 | 12 | 12 | 12 | 12 | 12 |
| T | tlog-perturbed | 35 | 34 | 30 | 18 | 16 | 35 |
| T | texp-perturbed | 36 | 33 | 31 | 14 | 15 | 36 |
| A | agen | 3 | 3 | 3 | 2 | 0 | 3 |
| B | bbottom-perturbed | 27 | 0 | 0 | 0 | 0 | 13 |
| M | alg | 34 | 34 | 3 | 8 | 0 | 33 |
| M | logxy | 31 | 31 | 0 | 13 | 0 | 0 |
| M | asinh | 6 | 6 | 6 | 6 | 0 | 6 |
| M | ratb | 22 | 22 | 11 | 9 | 0 | 22 |
| M | constk | 32 | 32 | 0 | 0 | 0 | 0 |
| M | tower2 | 35 | 35 | 3 | 35 | 0 | 24 |
| M | lccontent | 30 | 29 | 0 | 7 | 0 | 14 |
| M | expcube | 4 | 4 | 4 | 4 | 4 | 4 |
| M | mrandom | 33 | 0 | 0 | 0 | 0 | 0 |
| M | mperturb | 32 | 2 | 0 | 0 | 0 | 1 |
| B | bgenus | 24 | 24 | 24 | 24 | 0 | 21 |
| B | bcurve | 20 | 20 | 20 | 6 | 11 | 20 |
| B | bgenus-perturbed | 34 | 0 | 0 | 0 | 0 | 11 |
| B | bcurve-perturbed | 31 | 0 | 0 | 0 | 0 | 11 |
| T | classic-nonelem | 8 | 7 | 6 | 7 | 5 | 7 |
| M | classic-nonelem | 3 | 2 | 0 | 2 | 2 | 2 |
| A | classic-nonelem | 2 | 2 | 0 | 0 | 1 | 2 |
| A | trager-binomial-xk | 25 | 5 | 4 | 4 | 14 | 25 |
| A | trager-third-kind | 16 | 4 | 0 | 0 | 0 | 16 |
| A | trager-rational-times-radical | 33 | 2 | 0 | 0 | 5 | 21 |
| A | trager-chebyshev | 21 | 6 | 2 | 2 | 9 | 21 |
| A | trager-goursat | 15 | 7 | 6 | 2 | 2 | 10 |
| A | trager-pseudo-elliptic-recip | 4 | 4 | 4 | 0 | 0 | 4 |
| A | trager-genus-hyper | 12 | 0 | 0 | 0 | 0 | 3 |
| A | trager-genus-super | 8 | 0 | 0 | 0 | 4 | 5 |
| A | trager-rational-radicand | 8 | 2 | 0 | 0 | 0 | 6 |
| A | trager-torsion-heavy | 18 | 0 | 0 | 0 | 2 | 5 |
| A | trager-deriv-alg | 11 | 11 | 11 | 6 | 0 | 11 |
| A | trager-deriv-log | 11 | 11 | 11 | 0 | 0 | 8 |
| A | trager-random | 36 | 2 | 0 | 0 | 0 | 18 |
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

- sympy `tlogp-1-6`: `(x^2*Log[x] + 2*x*Log[1 + x] - Log[1 + x]^2 - x*Log[1 + x]^2)/(x^2*(1 + x))`
- sympy `texpp-1-8`: `(E^x*x + Log[x] + x*Log[x])/(1 + x)^2`
- sympy `texpp-1-11`: `(1 + E^x + x - E^x*x^2 + Log[x] + 2*E^x*Log[x] + E^(2*x)*Log[x])/((1 + E^x)^2*(1 + x))`
- Wolfram `bcurvep-1-10`: `(6*E^(x*Sqrt[3 + x^3])*x + 5*E^(x*Sqrt[3 + x^3])*x^4 - 2*E^(x*Sqrt[3 + x^3])*Sqrt[3 + x^3] + 2*x^2*Sqrt[3 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[3 + x^3])`
- sympy `classic-5`: `Log[x]/(1 + x)`
- sympy `trager-binom-008`: `x^3/(-1 + x^3)^(1/4)`
- sympy `trager-binom-009`: `x/(-1 + x^4)^(1/3)`
- sympy `trager-binom-019`: `1/(-1 + x^2)^(1/3)`
- sympy `trager-binom-020`: `x^3/(-1 + x^6)^(1/4)`
- sympy `trager-binom-024`: `1/(-1 + x^6)^(1/3)`
- sympy `trager-binom-025`: `x^3/Sqrt[-3 + x^3]`
- sympy `trager-binom-026`: `1/Sqrt[-1 + x^6]`
- sympy `trager-binom-038`: `x/(-3 + x^4)^(1/4)`
- sympy `trager-binom-039`: `x^3/Sqrt[-3 + x^6]`
- sympy `trager-binom-040`: `1/(-3 + x^6)^(1/3)`
- sympy `trager-binom-044`: `1/Sqrt[-3 + x^6]`
- sympy `trager-cheb-002`: `x^2*Sqrt[-1 + 2*x^4]`
- sympy `trager-cheb-004`: `x/(-1 + 2*x^3)^(3/4)`
- sympy `trager-cheb-005`: `(-1 + 2*x^4)^(1/4)`
- sympy `trager-cheb-014`: `x^3*(-1 + 2*x^3)^(2/3)`
- sympy `trager-cheb-022`: `x*(-1 + 2*x^3)^(3/2)`
- sympy `trager-cheb-030`: `x^4*(-1 + 2*x^2)^(1/4)`
- sympy `trager-cheb-040`: `x^4/(-1 + 2*x^4)^(1/3)`
- Wolfram `trager-goursat-007`: `x/((-2 + x^3)*(1 + x^3)^(1/3))`
- Wolfram `trager-goursat-009`: `(-2 + x^3)/((1 + x^3)^(2/3)*(4 + x^3))`
- sympy `trager-goursat-011`: `x^2/(-1 + x^4)^(3/4)`
- Maxima `trager-deriv-alg-03`: `(-22 - 22*x + 10*x^2 + 10*x^3)/(4*(2 - x - 2*x^2 + x^3)^(7/4))`
- sympy `trager-deriv-alg-07`: `(6 - 24*x + 42*x^4 - 48*x^5)/(6*(-1 + 3*x^4)^(11/6))`
- sympy `tlogp-2-6`: `(Log[x] + 6*x*Log[1 + x] + 3*Log[1 + x]^2 + 3*x*Log[1 + x]^2)/(1 + x)`
- sympy `texpp-2-5`: `(2*E^(2*x)*x + 4*E^(2*x)*x^2 + 2*E^(2*x)*x^3 + Log[x])/(1 + x)`
- Wolfram `bgenusp-2-11`: `(-3*E^(2*x) + 6*E^(2*x)*x - E^(2*x)*x^2 + 6*E^(2*x)*x^3 + 3*x^2*(1 + x^2)^(2/3)*Log[x + Sqrt[1 + x^3]])/(3*x^2*(1 + x^2)^(2/3))`
- Wolfram `bcurvep-2-3`: `(2*E^(x*Sqrt[1 + x^3])*x + 5*E^(x*Sqrt[1 + x^3])*x^4 - 2*E^(x*Sqrt[1 + x^3])*Sqrt[1 + x^3] + 2*x^2*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[1 + x^3])`
- sympy `tlogp-3-12`: `(3*Log[x] + 2*x*Log[x] + 6*Log[3 + 2*x]^2 + 12*x*Log[3 + 2*x]^2 + 6*x^2*Log[3 + 2*x]^2 + 3*Log[3 + 2*x]^3 + 5*x*Log[3 + 2*x]^3 + 2*x^2*Log[3 + 2*x]^3)/((1 + x)*`
- sympy `texpp-3-1`: `(3*E^(x/2) + 3*E^(x/2)*x + 2*Log[x] + 6*E^(x/2)*Log[x])/(2*(1 + 3*E^(x/2))*(1 + x))`
- sympy `texpp-3-5`: `(-3 - 4*x - x^2 + 4*E^x*Log[x] + 4*E^x*x*Log[x] + E^x*x^2*Log[x])/(E^x*(1 + x)*(2 + x)^2)`
- Wolfram `bcurvep-3-6`: `(2*E^(x*Sqrt[1 + x^3])*x + 5*E^(x*Sqrt[1 + x^3])*x^4 - 2*E^(x*Sqrt[1 + x^3])*Sqrt[1 + x^3] + 2*x^2*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*x^2*Sqrt[1 + x^3])`
