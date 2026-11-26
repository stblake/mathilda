# Integration shootout: 262-hard-integrals (2026-11-26)

Re-tabulated from the run of 2026-11-26: 1 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

955 cases, timeout 30 s per run, wall 28.4 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (330) | A algebraic (316) | M mixed (radical on top) (215) | B radical at the bottom (94) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (330) | 95% (299) | 100% (215) | 99% (93) | 98% (937) |
| Mathilda:ParallelMixedTower | 99% (328) | 81% (255) | 0% (1) | 61% (57) | 67% (641) |
| Mathilda:RischTranscendental | 100% (330) | · | · | · | 100% (330) |
| Mathilda:Trager | · | 93% (294) | · | · | 93% (294) |
| Mathilda:Bronstein | · | · | 100% (215) | 99% (93) | 100% (308) |
| FriCAS | 100% (330) | 80% (253) | 0% (0) | 45% (42) | 65% (625) |
| Axiom | 100% (330) | 73% (230) | 0% (0) | 0% (0) | 59% (560) |
| Maxima | 99% (327) | 69% (219) | 0% (0) | 52% (49) | 62% (595) |
| sympy | 89% (293) | 69% (217) | 0% (0) | 0% (0) | 53% (510) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 937 | 0 | 0 | 0 | 18 | 0 | 0 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 641 | 0 | 0 | 229 | 85 | 0 | 0 | 0 | 0.13 |
| Mathilda:RischTranscendental | 330 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.01 |
| Mathilda:Trager | 294 | 0 | 0 | 3 | 19 | 0 | 0 | 0 | 0.70 |
| Mathilda:Bronstein | 308 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.05 |
| FriCAS | 625 | 0 | 0 | 0 | 108 | 218 | 4 | 0 | 0.10 |
| Axiom | 560 | 0 | 0 | 5 | 45 | 343 | 0 | 2 | 0.17 |
| Maxima | 595 | 0 | 0 | 357 | 0 | 0 | 0 | 3 | 0.03 |
| sympy | 293 | 217 | 0 | 290 | 135 | 3 | 2 | 15 | 2.08 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 937 | 937 | 0 | 0 | 0 | 0 | 937 | 0 |
| Mathilda:ParallelMixedTower | 641 | 641 | 0 | 0 | 0 | 0 | 641 | 0 |
| Mathilda:RischTranscendental | 330 | 330 | 0 | 0 | 0 | 0 | 330 | 0 |
| Mathilda:Trager | 294 | 294 | 0 | 0 | 0 | 0 | 294 | 0 |
| Mathilda:Bronstein | 308 | 308 | 0 | 0 | 0 | 0 | 308 | 0 |
| FriCAS | 629 | 625 | 0 | 0 | 4 | 0 | 618 | 7 |
| Axiom | 562 | 560 | 2 | 0 | 0 | 0 | 560 | 2 |
| Maxima | 598 | 595 | 3 | 0 | 0 | 0 | 595 | 3 |
| sympy | 527 | 510 | 15 | 0 | 2 | 0 | 390 | 120 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 330 | 330 | 330 | 330 | 327 | 293 |
| A | cA-mathilda | 298 | 298 | 238 | 217 | 216 | 214 |
| A | cA-axiom | 15 | 2 | 15 | 13 | 0 | 0 |
| A | cA-maxima | 3 | 0 | 0 | 0 | 3 | 3 |
| M | cM-mathilda | 215 | 215 | 0 | 0 | 0 | 0 |
| B | cB-mathilda | 86 | 86 | 42 | 0 | 41 | 0 |
| B | cB-maxima | 8 | 8 | 0 | 0 | 8 | 0 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (2)

- `cA-201-105` AUTO-GAP Mathilda:Trager: `(-1 + x^6)/((-3 + x^6)^(1/6)*(1 + x^6))`
- `cB-214-388` AUTO-GAP Mathilda:Bronstein: `Log[x]/(x^6*(-3 + 2*x^6)^(1/6))`

## Only one system closes it

- Mathilda: 234
- Maxima: 3
- FriCAS: 2

## Wrong answers from other systems (their bugs)

- sympy `cT-201-71`: `(3*Log[3 + x] + 4*x*Log[3 + x] + 5*x^2*Log[3 + x] - 2*x^3*Log[3 + x] + Log[2 + x^3] + 2*x*Log[2 + x^3] - 2*x^3*Log[2 + x^3] - x^4*Log[2 + x^3])/((1 + x)^2*(1 + `
- sympy `cT-202-33`: `(x^2*Log[1 + x^2])/(-2 + x^3)^2`
- sympy `cT-202-36`: `x*(-(((-6 + x^3)*Log[x])/(3 + x^3)^2) - (2*Log[2 + x + x^2])/(3 + x^2)^2)`
- sympy `cT-203-8`: `x*((-2*Log[x])/(1 + x^2)^2 - (3*x*Log[(-2 + x)/(2 + x)])/(2 + x^3)^2)`
- sympy `cT-203-31`: `(x^2*Log[1 + x^2])/(2 + x^3)^2`
- sympy `cT-203-60`: `((-2 + x^2)*(Log[3 + x^2] + Log[-2 + x^3]))/(2 + x^2)^2`
- sympy `cT-204-40`: `(x*(-6 + x^3)*Log[-1 + x^3])/(3 + x^3)^2`
- sympy `cT-205-36`: `Log[-2 + x^3]/(3 + x)^3`
- sympy `cT-205-68`: `x^2*(-(Log[x]/(1 + x^3)^2) - Log[2 + x^2]/(2 + x^3)^2)`
- sympy `cT-206-48`: `-(Log[3 + x^2]/(-1 + x)^3) - Log[-2 + x^3]/(-2 + x)^3`
- Maxima `cT-207-5`: `x*((-3*x*Log[(-3 + x)/(1 + x)])/(-1 + x^3)^2 - ((-6 + x^3)*Log[1 + (2 + x)^(-1)])/(3 + x^3)^2)`
- sympy `cT-208-20`: `(x*(-4 + x^3)*Log[1 + x^3])/(2 + x^3)^2`
- sympy `cT-208-35`: `(x^2*Log[(2 + x)/(-3 + x)])/(2 + x^3)^2`
- sympy `cT-209-22`: `(-3*x^2*Log[1 + x^2])/(-2 + x^3)^2 - ((-3 - 4*x - x^2 + 2*x^3)*Log[1 + x + x^2])/((2 + x)^2*(1 + x^2)^2)`
- sympy `cT-209-89`: `(x*(-4 + x^3)*Log[-1 + x^3])/(2 + x^3)^2`
- Maxima `cT-210-29`: `x*(-(((-4 + x^3)*Log[(2 + x)/(-1 + x)])/(2 + x^3)^2) - (2*Log[1 + x + x^2])/(3 + x^2)^2)`
- Maxima `cT-212-49`: `(-(((-1 + 2*x)*(-5 + 2*x + x^2)*Log[1 + (-3 + x)^(-1)])/(-3 + x)^2) - (-1 + x)*(1 + x)*Log[1 + x^3])/(1 + x^2)^2`
- sympy `cT-212-76`: `-(((1 + x^2)*Log[-1 + x])/((-1 + x)^2*(1 + x)^2)) - (x*(-2 + x^3)*Log[-2 + x^3])/(1 + 2*x^3 + x^6)`
- Axiom `cA-202-43`: `(-2 + x^6)/((2 + x^6)*(3 + x^6)^(1/6))`
- Axiom `cA-205-117`: `(-1 + x^6)/((1 + x^6)*(3 + x^6)^(1/6))`
