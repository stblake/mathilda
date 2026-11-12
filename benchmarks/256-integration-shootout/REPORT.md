# Integration shootout (2026-11-12)

1151 cases, timeout 30 s per run, wall 12.4 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | R rational (72) | T transcendental (217) | A algebraic (400) | M mixed (radical on top) (408) | B radical at the bottom (54) | all |
|---|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (72) | 97% (210) | 57% (228) | 84% (342) | 48% (26) | 76% (878) |
| Mathilda:ParallelMixedTower | 0% (0) | 64% (138) | 52% (208) | 79% (322) | 37% (20) | 60% (688) |
| Mathilda:RischTranscendental | 100% (72) | 66% (143) | · | · | · | 74% (215) |
| Mathilda:Trager | · | · | 53% (211) | · | · | 53% (211) |
| Mathilda:Bronstein | · | · | · | 83% (337) | 11% (6) | 74% (343) |
| FriCAS | 99% (71) | 94% (204) | 51% (204) | 40% (165) | 41% (22) | 58% (666) |
| Maxima | 100% (72) | 77% (167) | 36% (145) | 42% (170) | 43% (23) | 50% (577) |
| sympy | 100% (72) | 80% (174) | 32% (130) | 13% (54) | 24% (13) | 38% (443) |
| Wolfram | 100% (72) | 100% (216) | 83% (331) | 61% (250) | 72% (39) | 79% (908) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 778 | 100 | 208 | 55 | 10 | 0 | 0 | 0 | 0.07 |
| Mathilda:ParallelMixedTower | 688 | 0 | 244 | 167 | 52 | 0 | 0 | 0 | 0.05 |
| Mathilda:RischTranscendental | 208 | 7 | 63 | 11 | 0 | 0 | 0 | 0 | 0.00 |
| Mathilda:Trager | 211 | 0 | 189 | 0 | 0 | 0 | 0 | 0 | 0.03 |
| Mathilda:Bronstein | 343 | 0 | 52 | 66 | 1 | 0 | 0 | 0 | 0.12 |
| FriCAS | 597 | 69 | 0 | 185 | 43 | 256 | 0 | 1 | 0.07 |
| Maxima | 527 | 50 | 0 | 570 | 2 | 0 | 0 | 2 | 0.02 |
| sympy | 348 | 95 | 0 | 427 | 240 | 1 | 8 | 32 | 0.12 |
| Wolfram | 657 | 251 | 0 | 227 | 13 | 0 | 1 | 2 | 0.15 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Maxima | sympy | Wolfram |
|---|---|---|---|---|---|---|---|
| R | rherm | 36 | 36 | 36 | 36 | 36 | 36 |
| R | rlog | 36 | 36 | 35 | 36 | 36 | 36 |
| T | tlog | 33 | 33 | 33 | 33 | 33 | 33 |
| T | texp | 33 | 33 | 32 | 26 | 33 | 33 |
| T | tmixed | 15 | 15 | 15 | 15 | 15 | 15 |
| T | tnested | 13 | 13 | 13 | 10 | 13 | 13 |
| T | ttrig | 18 | 18 | 18 | 18 | 18 | 18 |
| T | tconstk | 26 | 26 | 26 | 26 | 26 | 26 |
| T | tlog-perturbed | 35 | 32 | 30 | 18 | 16 | 35 |
| T | texp-perturbed | 36 | 33 | 31 | 14 | 15 | 36 |
| A | agen | 36 | 36 | 35 | 33 | 20 | 35 |
| B | bbottom | 20 | 20 | 20 | 17 | 12 | 20 |
| B | bbottom-perturbed | 33 | 5 | 1 | 5 | 0 | 18 |
| M | alg | 34 | 34 | 3 | 8 | 0 | 33 |
| M | logxy | 31 | 31 | 0 | 13 | 0 | 0 |
| M | asinh | 35 | 35 | 35 | 28 | 5 | 35 |
| M | genus1 | 35 | 35 | 35 | 19 | 8 | 35 |
| M | algres | 35 | 35 | 34 | 0 | 0 | 35 |
| M | ratb | 30 | 30 | 19 | 17 | 2 | 30 |
| M | constk | 32 | 32 | 0 | 0 | 0 | 0 |
| M | tower2 | 35 | 35 | 3 | 35 | 0 | 24 |
| M | lccontent | 35 | 33 | 0 | 12 | 1 | 17 |
| M | expcube | 36 | 36 | 34 | 36 | 35 | 36 |
| M | mrandom | 33 | 0 | 0 | 0 | 0 | 0 |
| M | mperturb | 34 | 4 | 2 | 0 | 1 | 3 |
| T | classic-nonelem | 8 | 7 | 6 | 7 | 5 | 7 |
| M | classic-nonelem | 3 | 2 | 0 | 2 | 2 | 2 |
| A | classic-nonelem | 2 | 2 | 0 | 0 | 1 | 2 |
| B | classic-nonelem | 1 | 1 | 1 | 1 | 1 | 1 |
| A | trager-binomial-xk | 45 | 22 | 20 | 21 | 32 | 45 |
| A | trager-third-kind | 30 | 18 | 14 | 13 | 4 | 30 |
| A | trager-rational-times-radical | 50 | 19 | 16 | 17 | 19 | 38 |
| A | trager-chebyshev | 40 | 23 | 20 | 20 | 25 | 40 |
| A | trager-goursat | 20 | 12 | 11 | 4 | 4 | 15 |
| A | trager-pseudo-elliptic-recip | 10 | 10 | 10 | 0 | 0 | 10 |
| A | trager-even-substitution | 6 | 6 | 6 | 0 | 0 | 6 |
| A | trager-genus-hyper | 20 | 6 | 4 | 4 | 4 | 11 |
| A | trager-genus-super | 10 | 1 | 1 | 1 | 6 | 7 |
| A | trager-rational-radicand | 20 | 14 | 12 | 5 | 0 | 18 |
| A | trager-torsion-heavy | 18 | 0 | 0 | 0 | 2 | 5 |
| A | trager-deriv-alg | 18 | 18 | 18 | 13 | 3 | 18 |
| A | trager-deriv-log | 23 | 23 | 23 | 0 | 0 | 17 |
| A | trager-random | 52 | 18 | 14 | 14 | 10 | 34 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 67
- Maxima: 1

## Wrong answers from other systems (their bugs)

- sympy `tlogp-1-6`: `(x^2*Log[x] + 2*x*Log[1 + x] - Log[1 + x]^2 - x*Log[1 + x]^2)/(x^2*(1 + x))`
- sympy `texpp-1-8`: `(E^x*x + Log[x] + x*Log[x])/(1 + x)^2`
- sympy `texpp-1-11`: `(1 + E^x + x - E^x*x^2 + Log[x] + 2*E^x*Log[x] + E^(2*x)*Log[x])/((1 + E^x)^2*(1 + x))`
- sympy `classic-5`: `Log[x]/(1 + x)`
- sympy `trager-binom-008`: `x^3/(-1 + x^3)^(1/4)`
- sympy `trager-binom-009`: `x/(-1 + x^4)^(1/3)`
- sympy `trager-binom-019`: `1/(-1 + x^2)^(1/3)`
- sympy `trager-binom-020`: `x^3/(-1 + x^6)^(1/4)`
- sympy `trager-binom-024`: `1/(-1 + x^6)^(1/3)`
- sympy `trager-binom-025`: `x^3/Sqrt[-3 + x^3]`
- sympy `trager-binom-026`: `1/Sqrt[-1 + x^6]`
- sympy `trager-binom-027`: `x/(-1 + x^4)^(1/4)`
- sympy `trager-binom-038`: `x/(-3 + x^4)^(1/4)`
- sympy `trager-binom-039`: `x^3/Sqrt[-3 + x^6]`
- sympy `trager-binom-040`: `1/(-3 + x^6)^(1/3)`
- sympy `trager-binom-043`: `x/Sqrt[-1 + x^6]`
- sympy `trager-binom-044`: `1/Sqrt[-3 + x^6]`
- sympy `trager-cheb-002`: `x^2*Sqrt[-1 + 2*x^4]`
- sympy `trager-cheb-004`: `x/(-1 + 2*x^3)^(3/4)`
- sympy `trager-cheb-005`: `(-1 + 2*x^4)^(1/4)`
- sympy `trager-cheb-014`: `x^3*(-1 + 2*x^3)^(2/3)`
- sympy `trager-cheb-022`: `x*(-1 + 2*x^3)^(3/2)`
- sympy `trager-cheb-027`: `1/(x*(3 - x^2)^(1/3))`
- sympy `trager-cheb-030`: `x^4*(-1 + 2*x^2)^(1/4)`
- sympy `trager-cheb-040`: `x^4/(-1 + 2*x^4)^(1/3)`
- Wolfram `trager-goursat-007`: `x/((-2 + x^3)*(1 + x^3)^(1/3))`
- Wolfram `trager-goursat-009`: `(-2 + x^3)/((1 + x^3)^(2/3)*(4 + x^3))`
- sympy `trager-goursat-011`: `x^2/(-1 + x^4)^(3/4)`
- Maxima `trager-ratradicand-016`: `Sqrt[(1 + x)/(-1 + x)]/(-2 + x)`
- Maxima `trager-deriv-alg-03`: `(-22 - 22*x + 10*x^2 + 10*x^3)/(4*(2 - x - 2*x^2 + x^3)^(7/4))`
- sympy `trager-deriv-alg-07`: `(6 - 24*x + 42*x^4 - 48*x^5)/(6*(-1 + 3*x^4)^(11/6))`
- FriCAS `trager-random-048`: `(2 + 3*x - 3*x^2)/(Sqrt[1 + 2*x]*(2 - x + 3*x^2))`
- sympy `tlogp-2-6`: `(Log[x] + 6*x*Log[1 + x] + 3*Log[1 + x]^2 + 3*x*Log[1 + x]^2)/(1 + x)`
- sympy `texpp-2-5`: `(2*E^(2*x)*x + 4*E^(2*x)*x^2 + 2*E^(2*x)*x^3 + Log[x])/(1 + x)`
- sympy `tlogp-3-12`: `(3*Log[x] + 2*x*Log[x] + 6*Log[3 + 2*x]^2 + 12*x*Log[3 + 2*x]^2 + 6*x^2*Log[3 + 2*x]^2 + 3*Log[3 + 2*x]^3 + 5*x*Log[3 + 2*x]^3 + 2*x^2*Log[3 + 2*x]^3)/((1 + x)*`
- sympy `texpp-3-1`: `(3*E^(x/2) + 3*E^(x/2)*x + 2*Log[x] + 6*E^(x/2)*Log[x])/(2*(1 + 3*E^(x/2))*(1 + x))`
- sympy `texpp-3-5`: `(-3 - 4*x - x^2 + 4*E^x*Log[x] + 4*E^x*x*Log[x] + E^x*x^2*Log[x])/(E^x*(1 + x)*(2 + x)^2)`
