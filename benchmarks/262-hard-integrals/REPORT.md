# Integration shootout: 262-hard-integrals (2026-11-25)

1000 cases, timeout 30 s per run, wall 81.7 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (330) | A algebraic (336) | M mixed (radical on top) (216) | B radical at the bottom (118) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (330) | 94% (317) | 100% (216) | 91% (107) | 97% (970) |
| Mathilda:ParallelMixedTower | 99% (328) | 81% (271) | 0% (1) | 51% (60) | 66% (660) |
| Mathilda:RischTranscendental | 100% (330) | · | · | · | 100% (330) |
| Mathilda:Trager | · | 93% (312) | · | · | 93% (312) |
| Mathilda:Bronstein | · | · | 100% (216) | 100% (118) | 100% (334) |
| FriCAS | 100% (330) | 80% (268) | 0% (0) | 41% (48) | 65% (646) |
| Axiom | 100% (330) | 72% (241) | 0% (0) | 0% (0) | 57% (571) |
| Maxima | 99% (327) | 71% (238) | 0% (0) | 62% (73) | 64% (638) |
| sympy | 89% (293) | 68% (228) | 0% (0) | 0% (0) | 52% (521) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 970 | 0 | 0 | 30 | 0 | 0 | 0 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 660 | 0 | 0 | 235 | 105 | 0 | 0 | 0 | 0.15 |
| Mathilda:RischTranscendental | 330 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.01 |
| Mathilda:Trager | 312 | 0 | 0 | 3 | 21 | 0 | 0 | 0 | 0.76 |
| Mathilda:Bronstein | 334 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.08 |
| FriCAS | 646 | 0 | 0 | 0 | 131 | 219 | 4 | 0 | 0.10 |
| Axiom | 571 | 0 | 0 | 5 | 45 | 377 | 0 | 2 | 0.17 |
| Maxima | 638 | 0 | 0 | 359 | 0 | 0 | 0 | 3 | 0.04 |
| sympy | 293 | 228 | 0 | 235 | 224 | 3 | 2 | 15 | 2.50 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 970 | 970 | 0 | 0 | 0 | 0 | 970 | 0 |
| Mathilda:ParallelMixedTower | 660 | 660 | 0 | 0 | 0 | 0 | 660 | 0 |
| Mathilda:RischTranscendental | 330 | 330 | 0 | 0 | 0 | 0 | 330 | 0 |
| Mathilda:Trager | 312 | 312 | 0 | 0 | 0 | 0 | 312 | 0 |
| Mathilda:Bronstein | 334 | 334 | 0 | 0 | 0 | 0 | 334 | 0 |
| FriCAS | 650 | 646 | 0 | 0 | 4 | 0 | 639 | 7 |
| Axiom | 573 | 571 | 2 | 0 | 0 | 0 | 571 | 2 |
| Maxima | 641 | 638 | 3 | 0 | 0 | 0 | 638 | 3 |
| sympy | 538 | 521 | 15 | 0 | 2 | 0 | 395 | 126 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 330 | 330 | 330 | 330 | 327 | 293 |
| A | cA-mathilda | 318 | 318 | 253 | 228 | 235 | 225 |
| A | cA-axiom | 15 | 2 | 15 | 13 | 0 | 0 |
| A | cA-maxima | 3 | 0 | 0 | 0 | 3 | 3 |
| M | cM-mathilda | 216 | 216 | 0 | 0 | 0 | 0 |
| B | cB-mathilda | 110 | 110 | 48 | 0 | 65 | 0 |
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
- Maxima `cT-207-5`: `x*((-3*x*Log[(-3 + x)/(1 + x)])/(-1 + x^3)^2 - ((-6 + x^3)*Log[1 + (2 + x)^(-1)])/(3 + x^3)^2)`
- sympy `cT-208-20`: `-((x*(-4 + x^3)*Log[1 + x^3])/(2 + x^3)^2)`
- sympy `cT-208-35`: `(-3*x^2*Log[(2 + x)/(-3 + x)])/(2 + x^3)^2`
- sympy `cT-209-22`: `(-3*x^2*Log[1 + x^2])/(-2 + x^3)^2 - ((-3 - 4*x - x^2 + 2*x^3)*Log[1 + x + x^2])/((2 + x)^2*(1 + x^2)^2)`
- sympy `cT-209-89`: `-((x*(-4 + x^3)*Log[-1 + x^3])/(2 + x^3)^2)`
- Maxima `cT-210-29`: `x*(-(((-4 + x^3)*Log[(2 + x)/(-1 + x)])/(2 + x^3)^2) - (2*Log[1 + x + x^2])/(3 + x^2)^2)`
- Maxima `cT-212-49`: `(-(((-1 + 2*x)*(-5 + 2*x + x^2)*Log[1 + (-3 + x)^(-1)])/(-3 + x)^2) - (-1 + x)*(1 + x)*Log[1 + x^3])/(1 + x^2)^2`
- sympy `cT-212-76`: `-(((1 + x^2)*Log[-1 + x])/((-1 + x)^2*(1 + x)^2)) - (x*(-2 + x^3)*Log[-2 + x^3])/(1 + 2*x^3 + x^6)`
- Axiom `cA-202-43`: `(-2 + x^6)/((2 + x^6)*(3 + x^6)^(1/6))`
- Axiom `cA-205-117`: `(-1 + x^6)/((1 + x^6)*(3 + x^6)^(1/6))`
