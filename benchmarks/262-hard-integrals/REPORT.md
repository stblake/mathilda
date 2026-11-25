# Integration shootout: 262-hard-integrals (2026-11-25)

1000 cases, timeout 30 s per run, wall 81.7 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | T transcendental (330) | A algebraic (336) | M mixed (radical on top) (216) | B radical at the bottom (118) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (330) | 94% (317) | 100% (216) | 91% (107) | 97% (970) |
| Mathilda:ParallelMixedTower | 99% (328) | 81% (271) | 0% (1) | 51% (60) | 66% (660) |
| Mathilda:RischTranscendental | 100% (330) | · | · | · | 100% (330) |
| Mathilda:Trager | · | 93% (312) | · | · | 93% (312) |
| Mathilda:Bronstein | · | · | 100% (216) | 100% (118) | 100% (334) |
| FriCAS | 98% (323) | 79% (266) | 0% (0) | 43% (51) | 64% (640) |
| Axiom | 100% (330) | 72% (243) | 0% (0) | 0% (0) | 57% (573) |
| Maxima | 100% (330) | 71% (238) | 0% (0) | 62% (73) | 64% (641) |
| sympy | 89% (293) | 35% (116) | 0% (0) | 0% (0) | 41% (409) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 970 | 0 | 0 | 30 | 0 | 0 | 0 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 660 | 0 | 0 | 235 | 105 | 0 | 0 | 0 | 0.15 |
| Mathilda:RischTranscendental | 330 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.01 |
| Mathilda:Trager | 312 | 0 | 0 | 3 | 21 | 0 | 0 | 0 | 0.76 |
| Mathilda:Bronstein | 334 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.08 |
| FriCAS | 640 | 0 | 0 | 0 | 130 | 230 | 0 | 0 | 0.10 |
| Axiom | 573 | 0 | 0 | 5 | 42 | 380 | 0 | 0 | 0.17 |
| Maxima | 641 | 0 | 0 | 359 | 0 | 0 | 0 | 0 | 0.04 |
| sympy | 293 | 116 | 0 | 235 | 224 | 3 | 114 | 15 | 2.09 |

## Cross-check in Mathilda (answers that failed their own check)

Own check: `bad` = a sample point failed; `na` = no point could be evaluated. Verdict of PossibleZeroQ[D[F, x] - f] in Mathilda.

| run | own | n | zero (verified) | nonzero (wrong) | undecided | untranslatable |
|---|---|---|---|---|---|---|
| FriCAS | bad | 1 | 1 | 0 | 0 | 0 |
| FriCAS | na | 6 | 6 | 0 | 0 | 0 |
| sympy | bad | 124 | 0 | 15 | 109 | 0 |
| sympy | na | 18 | 13 | 0 | 4 | 1 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 330 | 330 | 323 | 330 | 330 | 293 |
| A | cA-mathilda | 318 | 318 | 251 | 228 | 235 | 115 |
| A | cA-axiom | 15 | 2 | 15 | 15 | 0 | 0 |
| A | cA-maxima | 3 | 0 | 0 | 0 | 3 | 1 |
| M | cM-mathilda | 216 | 216 | 0 | 0 | 0 | 0 |
| B | cB-mathilda | 110 | 110 | 51 | 0 | 65 | 0 |
| B | cB-maxima | 8 | 8 | 0 | 0 | 8 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (14)

- `cA-201-105` AUTO-GAP Mathilda:Trager: `(-1 + x^6)/((-3 + x^6)^(1/6)*(1 + x^6))`
- `cA-201-110` AUTO-GAP Mathilda:Trager: `(-1 + x^6)/((1 + x^6)*(2 + x^6)^(1/6))`
- `cA-212-83` AUTO-GAP Mathilda:Trager: `(-2 + 2*x^3)^(2/3)/(2 + x^3)`
- `cB-201-559` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(3 + 2*x^6)^(1/6))`
- `cB-202-297` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(2 + x^6)^(5/6))`
- `cB-202-470` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(1 + 2*x^6)^(1/6))`
- `cB-203-340` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(-2 + x^6)^(5/6))`
- `cB-204-491` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(3 + 3*x^6)^(5/6))`
- `cB-206-513` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(-3 + 3*x^6)^(5/6))`
- `cB-207-392` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(3 + x^6)^(1/6))`
- `cB-210-330` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(2 + 3*x^6)^(5/6))`
- `cB-210-386` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^2*(2 + 2*x^6)^(5/6))`
- `cB-214-388` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(-3 + 2*x^6)^(1/6))`
- `cB-214-404` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(1 + 3*x^6)^(1/6))`

## Only one system closes it

- Mathilda: 232
- Maxima: 3

## Wrong answers from other systems (their bugs)

- sympy `cT-201-71`: `(3*Log[3 + x] + 4*x*Log[3 + x] + 5*x^2*Log[3 + x] - 2*x^3*Log[3 + x] + Log[2 + x^3] + 2*x*Log[2 + x^3] - 2*x^3*Log[2 + x^3] - x^4*Log[2 + x^3])/((1 + x)^2*(1 + `
- sympy `cT-202-33`: `(-3*x^2*Log[1 + x^2])/(-2 + x^3)^2`
- sympy `cT-202-36`: `x*(-(((-6 + x^3)*Log[x])/(3 + x^3)^2) - (2*Log[2 + x + x^2])/(3 + x^2)^2)`
- sympy `cT-203-8`: `x*((-2*Log[x])/(1 + x^2)^2 - (3*x*Log[(-2 + x)/(2 + x)])/(2 + x^3)^2)`
- sympy `cT-203-31`: `(x^2*Log[1 + x^2])/(2 + x^3)^2`
- sympy `cT-203-60`: `-(((-2 + x^2)*(Log[3 + x^2] + Log[-2 + x^3]))/(2 + x^2)^2)`
- sympy `cT-204-40`: `-((x*(-6 + x^3)*Log[-1 + x^3])/(3 + x^3)^2)`
- sympy `cT-205-36`: `(-2*Log[-2 + x^3])/(3 + x)^3`
- sympy `cT-205-68`: `3*x^2*(-(Log[x]/(1 + x^3)^2) - Log[2 + x^2]/(2 + x^3)^2)`
- sympy `cT-206-48`: `(-2*Log[3 + x^2])/(-1 + x)^3 - (2*Log[-2 + x^3])/(-2 + x)^3`
- sympy `cT-208-20`: `-((x*(-4 + x^3)*Log[1 + x^3])/(2 + x^3)^2)`
- sympy `cT-208-35`: `(-3*x^2*Log[(2 + x)/(-3 + x)])/(2 + x^3)^2`
- sympy `cT-209-22`: `(-3*x^2*Log[1 + x^2])/(-2 + x^3)^2 - ((-3 - 4*x - x^2 + 2*x^3)*Log[1 + x + x^2])/((2 + x)^2*(1 + x^2)^2)`
- sympy `cT-209-89`: `-((x*(-4 + x^3)*Log[-1 + x^3])/(2 + x^3)^2)`
- sympy `cT-212-76`: `-(((1 + x^2)*Log[-1 + x])/((-1 + x)^2*(1 + x)^2)) - (x*(-2 + x^3)*Log[-2 + x^3])/(1 + 2*x^3 + x^6)`
