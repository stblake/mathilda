# Integration shootout: 263-hard-integrals-ii (2026-11-25)

1850 cases, timeout 30 s per run, wall 70.8 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (500) | A algebraic (550) | M mixed (radical on top) (400) | B radical at the bottom (400) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (500) | 96% (526) | 100% (400) | 100% (400) | 99% (1826) |
| Mathilda:ParallelMixedTower | 99% (497) | 83% (454) | 0% (0) | 97% (388) | 72% (1339) |
| Mathilda:RischTranscendental | 100% (500) | · | · | · | 100% (500) |
| Mathilda:Trager | · | 94% (517) | · | · | 94% (517) |
| Mathilda:Bronstein | · | · | 100% (400) | 100% (400) | 100% (800) |
| FriCAS | 100% (500) | 78% (427) | 0% (0) | 97% (389) | 71% (1316) |
| Axiom | 100% (500) | 61% (333) | 0% (0) | 92% (370) | 65% (1203) |
| Maxima | 97% (486) | 73% (403) | 0% (0) | 60% (238) | 61% (1127) |
| sympy | 89% (447) | 70% (385) | 0% (0) | 29% (116) | 51% (948) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1826 | 0 | 0 | 24 | 0 | 0 | 0 | 0 | 0.13 |
| Mathilda:ParallelMixedTower | 1339 | 0 | 0 | 406 | 105 | 0 | 0 | 0 | 0.49 |
| Mathilda:RischTranscendental | 500 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.08 |
| Mathilda:Trager | 517 | 0 | 0 | 4 | 29 | 0 | 0 | 0 | 2.01 |
| Mathilda:Bronstein | 800 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.11 |
| FriCAS | 1316 | 0 | 0 | 11 | 103 | 420 | 0 | 0 | 0.18 |
| Axiom | 1203 | 0 | 0 | 34 | 89 | 520 | 0 | 4 | 0.15 |
| Maxima | 1127 | 0 | 0 | 697 | 0 | 7 | 0 | 19 | 0.03 |
| sympy | 557 | 391 | 0 | 540 | 320 | 0 | 14 | 28 | 1.96 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1826 | 1826 | 0 | 0 | 0 | 0 | 1826 | 0 |
| Mathilda:ParallelMixedTower | 1339 | 1339 | 0 | 0 | 0 | 0 | 1339 | 0 |
| Mathilda:RischTranscendental | 500 | 500 | 0 | 0 | 0 | 0 | 500 | 0 |
| Mathilda:Trager | 517 | 517 | 0 | 0 | 0 | 0 | 517 | 0 |
| Mathilda:Bronstein | 800 | 800 | 0 | 0 | 0 | 0 | 800 | 0 |
| FriCAS | 1316 | 1316 | 0 | 0 | 0 | 0 | 1309 | 7 |
| Axiom | 1207 | 1203 | 4 | 0 | 0 | 0 | 1203 | 0 |
| Maxima | 1146 | 1127 | 19 | 0 | 0 | 0 | 1127 | 15 |
| sympy | 990 | 948 | 28 | 0 | 14 | 0 | 741 | 207 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 494 | 494 | 494 | 494 | 480 | 442 |
| T | cT-fricas | 6 | 6 | 6 | 6 | 6 | 5 |
| A | cA-mathilda | 515 | 515 | 400 | 310 | 391 | 373 |
| A | cA-fricas | 24 | 23 | 24 | 22 | 2 | 2 |
| A | cA-maxima | 10 | 10 | 2 | 0 | 10 | 10 |
| A | cA-axiom | 1 | 1 | 1 | 1 | 0 | 0 |
| M | cM-mathilda | 400 | 400 | 0 | 0 | 0 | 0 |
| B | cB-mathilda | 400 | 400 | 389 | 370 | 238 | 116 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (23)

- `cA-312-12` AUTO-GAP Mathilda:Trager: `(-2 + 2*x^6)^(5/6)/(1 + x^6)`
- `cA-314-18` AUTO-GAP Mathilda:Trager: `(x^3*(-2 + 3*x^6)^(1/3))/(3 + x^6)`
- `cA-314-86` AUTO-GAP Mathilda:Trager: `(x^6*(-3 + 3*x^3)^(2/3))/(3 + x^3)`
- `cA-316-17` AUTO-GAP Mathilda:Trager: `x^6/((2 + x^3)*(1 + 3*x^3)^(1/3))`
- `cA-316-39` AUTO-GAP Mathilda:Trager: `(x^4*(-1 + 2*x^3)^(4/3))/(1 + x^3)`
- `cA-316-80` AUTO-GAP Mathilda:Trager: `(-1 + x)*x^4*(-1 + 2*x^6)^(1/6)`
- `cA-317-41` AUTO-GAP Mathilda:Trager: `(2*(1 + x)*(3 + 2*x^3)^(4/3))/x^2`
- `cA-317-104` AUTO-GAP Mathilda:Trager: `(x^4*(1 + x))/(-1 + 3*x^6)^(5/6)`
- `cA-319-94` AUTO-GAP Mathilda:Trager: `x^3/((3 + x^3)*(3 + 3*x^3)^(1/3))`
- `cA-324-122` AUTO-GAP Mathilda:Trager: `(x^6*(-2 + x^3)^(2/3))/(2 + x^3)`
- `cA-325-54` AUTO-GAP Mathilda:Trager: `((-2 + x)*x^4)/(-1 + 3*x^6)^(5/6)`
- `cA-326-26` AUTO-GAP Mathilda:Trager: `(-2 + x)/(x^2*(2 + x^6)^(5/6))`
- `cA-326-54` AUTO-GAP Mathilda:ParallelMixedTower: `1/((1 + 2*x^3)^(4/3)*(3 + x^6))`
- `cA-327-90` AUTO-GAP Mathilda:Trager: `(x^7*(1 + 2*x^3)^(4/3))/(1 + x^3)`
- `cA-329-86` AUTO-GAP Mathilda:Trager: `(-2 + x)*x^4*(-3 + 3*x^6)^(1/6)`
- `cA-334-59` AUTO-GAP Mathilda:Trager: `(x^4*(1 + 2*x))/(2 + x^6)^(5/6)`
- `cA-335-75` AUTO-GAP Mathilda:Trager: `(x^7*(2 + 2*x^3)^(4/3))/(2 + x^3)`
- `cA-336-100` AUTO-GAP Mathilda:Trager: `((1 + x)*(-2 + 2*x^6)^(1/6))/x^2`
- `cA-338-86` AUTO-GAP Mathilda:Trager: `((-1 + x)*x^4)/(-2 + x^6)^(5/6)`
- `cA-339-33` AUTO-GAP Mathilda:Trager: `(x^3*(-2 + 2*x^3)^(2/3))/(2 + x^3)`
- `cA-339-42` AUTO-GAP Mathilda:Trager: `(x^3*(2 + 2*x^3)^(2/3))/(2 + x^3)`
- `cA-340-68` AUTO-GAP Mathilda:Trager: `(-3 + 3*x^3)^(2/3)/(3 + x^3)`
- `cA-343-26` AUTO-GAP Mathilda:Trager: `((-1 + x)*x^4)/(-1 + 3*x^6)^(5/6)`

## Only one system closes it

- Mathilda: 430

## Wrong answers from other systems (their bugs)

- sympy `cT-301-6`: `-(((-3 + x^2)*Log[x])/(3 + x^2)^2) - (x*(-4 + x^3)*Log[1 + x + x^2])/(2 + x^3)^2`
- sympy `cT-301-9`: `-((x*(-6 + x^3)*Log[-1 + x^2])/(3 + x^3)^2)`
- sympy `cT-301-51`: `(-2*Log[3 + x])/(3 + x)^3 - (x*(-4 + x^3)*Log[-2 + x^3])/(2 + x^3)^2`
- sympy `cT-302-70`: `-(Log[-3 + x]/(-3 + x)^2) - ((-10 - 12*x - 3*x^2 + 2*x^3)*Log[2 + x^3])/((3 + x)^2*(2 + x^2)^2)`
- sympy `cT-303-49`: `x*((-2*Log[x])/(3 + x^2)^2 - ((-4 + x^3)*Log[3 + x^2])/(2 + x^3)^2)`
- Maxima `cT-304-3`: `(E^x*(-2 + x))/(-1 + x)^2 - (2*x*Log[(3 + x)/(-2 + x)])/(2 + x^2)^2`
- sympy `cT-305-23`: `-((x*(-4 + x^3)*Log[-2 + x^3])/(2 + x^3)^2)`
- sympy `cT-306-56`: `((10 + 12*x + 3*x^2 - 2*x^3)*Log[2 + x^3])/((3 + x)^2*(2 + x^2)^2)`
- sympy `cT-307-24`: `(E^x*(-4 + x))/(-3 + x)^2 - (2*(-2 + x)*(1 + x + x^2)*Log[2 + x^3])/((1 + x)^2*(2 + x^2)^2)`
- sympy `cT-308-28`: `-(((-3 + x^2)*Log[3 + x^2])/(3 + x^2)^2) - Log[-2 + x^3]/(-3 + x)^2`
- sympy `cT-308-53`: `-(Log[-2 + x^3]/(2 + x)^2)`
- sympy `cT-309-13`: `((1 - 4*x - 5*x^2 - 2*x^3)*Log[-2 + x^3])/((2 + x)^2*(1 + x^2)^2)`
- Maxima `cT-309-46`: `x*((E^x*(4 + 2*x - x^3 + x^4))/(2 + x^3)^2 - (2*Log[-1 + x^2])/(1 + x^2)^2)`
- sympy `cT-310-9`: `(x*(2 - x^3)*Log[2 + x^3])/(1 + x^3)^2`
- Maxima `cT-310-37`: `x*(-(((-2 + x^3)*Log[-2 + x^2])/(1 + 2*x^3 + x^6)) - ((-4 + x^3)*Log[-1 + x^2])/(2 + x^3)^2)`
- sympy `cT-311-44`: `-(((-2 + 2*x + x^2)*Log[1 + (-2 + x)^(-1)])/(3 + x + x^2)^2) - (x*(-4 + x^3)*Log[-2 + x])/(2 + x^3)^2`
- sympy `cT-312-50`: `-((x*(-4 + x^3)*Log[-2 + x])/(2 + x^3)^2)`
- sympy `cT-313-18`: `(x*(2 - x^3)*Log[-2 + x^3])/(1 + x^3)^2`
- sympy `cT-313-26`: `(-3*x^2*Log[-2 + x])/(2 + x^3)^2 - (2*Log[3 + x + x^2])/(-2 + x)^3`
- sympy `cT-313-28`: `(-2*(1 - x + x^2 + x^3)*Log[-2 + x^3])/((-1 + x)^2*(1 + x^2)^2)`
- sympy `cT-313-43`: `-(((1 + x^2)*Log[x])/((-1 + x)^2*(1 + x)^2)) - (3*x^2*Log[(-2 + x)/(2 + x)])/(2 + x^3)^2`
- sympy `cT-314-22`: `(-2*x*Log[(3 + x)/(-2 + x)])/(2 + x^2)^2 - Log[-2 + x^3]/(-2 + x)^2`
- Maxima `cT-314-46`: `(E^x*(-4 + x))/(-2 + x)^3 - (2*(-4 - 4*x - 2*x^2 + x^3)*Log[-1 + x^2])/((2 + x)^2*(2 + x^2)^2)`
- Maxima `cT-315-42`: `(-E^x - 3*E^x*x^2 + E^x*x^3 - 3*x^2*Log[(-1 + x)/(-3 + x)])/(-1 + x^3)^2`
- sympy `cT-316-18`: `(E^x*(1 + x))/(3 + x)^3 - (2*(4 - 4*x + 2*x^2 + x^3)*Log[2 + x^3])/((-2 + x)^2*(2 + x^2)^2)`
- sympy `cT-316-56`: `(E^x*(-4 + x))/(-2 + x)^3 - (x*(-6 + x^3)*Log[-1 + x^2])/(3 + x^3)^2`
- Maxima `cT-317-18`: `x*((-2*Log[-2 + x])/(3 + x^2)^2 - ((-2 + x^3)*Log[(2 + x)/(-1 + x)])/(1 + 2*x^3 + x^6))`
- sympy `cT-317-45`: `-(((1 + 4*x + x^2)*Log[x])/(1 + x + x^2)^2) + ((1 + x)^(-2) + 2/(1 + x^2)^2 - (1 + x^2)^(-1))*Log[-2 + x^3]`
- sympy `cT-318-34`: `(-2*Log[3 + x + x^2])/(2 + x)^3 - ((-2 + 2*x + x^2)*Log[2 + x^3])/(3 + x + x^2)^2`
- sympy `cT-318-40`: `(-3*x^2*Log[1 + x^2])/(2 + x^3)^2`
- sympy `cT-319-17`: `(-(-1 + x)^(-2) - 2/(1 + x^2)^2 + (1 + x^2)^(-1))*Log[x] - ((1 + x^2)*Log[-2 + x^3])/((-1 + x)^2*(1 + x)^2)`
- sympy `cT-320-38`: `x*((-3*x*Log[x])/(2 + x^3)^2 - (2*Log[-2 + x^3])/(2 + x^2)^2)`
- sympy `cT-320-62`: `-((x*(-6 + x^3)*Log[x])/(3 + x^3)^2) - (2*Log[-1 + x^2])/(3 + x)^3`
- sympy `cT-322-10`: `(-2*Log[x])/(-3 + x)^3 - ((-2 + x^2)*Log[-2 + x^3])/(2 + x^2)^2`
- Axiom `cA-313-80`: `(x^7*(3 + x^3)^(4/3))/(1 + x^6)`
- Axiom `cA-313-88`: `(x^7*(-3 + x^3)^(4/3))/(1 + x^6)`
- Axiom `cA-332-118`: `(x^4*(1 + x^3)^(4/3))/(2 + x^6)`
- Axiom `cA-337-73`: `x^7/((2 + x^3)^(2/3)*(2 + x^6))`
- Maxima `cB-304-10`: `(-3*(2 + x + x^3) + x*Log[2 + x^(-1)] + 2*x^2*Log[2 + x^(-1)] + 3*x^3*Log[2 + x^(-1)] + 6*x^4*Log[2 + x^(-1)])/((3*x + 6*x^2)*(2 + x + x^3)^(2/3))`
- Maxima `cB-304-35`: `(-1 - x^4 + x^4*Log[2 + x^(-1)] + 2*x^5*Log[2 + x^(-1)])/(x*(1 + 2*x)*(1 + x^4)^(3/4))`
- Maxima `cB-305-23`: `(E^x*(-1 + 2*x)*(-2 + x^2))/(2*Sqrt[1 - 2*x^2 + x^3])`
- Maxima `cB-306-12`: `(-2 - x^4 + x^4*Log[2 + x^(-1)] + 2*x^5*Log[2 + x^(-1)])/(x*(1 + 2*x)*(2 + x^4)^(3/4))`
- Maxima `cB-306-16`: `(E^Sqrt[-1 + x^4]*(2*x^3 + 2*x^4 - Sqrt[-1 + x^4]))/((1 + x)^2*Sqrt[-1 + x^4])`
- Maxima `cB-317-11`: `(-4 - 2*x^3 + 3*x^3*Log[1 + x^(-1)] + 3*x^4*Log[1 + x^(-1)])/(2*x*(1 + x)*Sqrt[2 + x^3])`
- Maxima `cB-317-37`: `(Sqrt[1 - 2*x^2 + x^3]*(-2 + 4*x^2 - 2*x^3 + 8*x*Log[2 + x] - 2*x^2*Log[2 + x] - 3*x^3*Log[2 + x]))/(2*(2 + x - 4*x^2 + x^4))`
- Maxima `cB-317-42`: `(Sqrt[1 - 2*x^2 + x^3]*(4 - 6*x^2 + 4*x^3 - 4*x^4 + 2*x^5 + 4*Log[x] - 18*x^2*Log[x] + 10*x^3*Log[x] + x^5*Log[x]))/(2*(-1 + x)*(2 + x^2)^2*(-1 - x + x^2))`
- Maxima `cB-318-8`: `(-4 - 4*x^2 - 2*x^3 + 4*x^2*Log[2 + x^(-1)] + 11*x^3*Log[2 + x^(-1)] + 6*x^4*Log[2 + x^(-1)])/(2*x*(1 + 2*x)*Sqrt[2 + 2*x^2 + x^3])`
- Maxima `cB-319-12`: `(Sqrt[1 - 2*x^2 + x^3]*(2 - 4*x^2 + 2*x^3 - 8*x*Log[2 + x] + 2*x^2*Log[2 + x] + 3*x^3*Log[2 + x]))/(2*(2 + x - 4*x^2 + x^4))`
- Maxima `cB-321-13`: `-(Sqrt[2 + x + x^4]/(x + 2*x^2)) + ((Sqrt[2 + x + x^4] + 4*x^3*Sqrt[2 + x + x^4])*Log[2 + x^(-1)])/(4 + 2*x + 2*x^4)`
- Maxima `cB-323-2`: `(E^(2*x)*(-5 - 2*x + 4*x^3 + 7*x^4 + 2*x^5))/Sqrt[-1 + x^4]`
- Maxima `cB-329-7`: `(E^x*x*(-2 - x + 4*x^4 + x^5))/Sqrt[-1 + x^4]`
