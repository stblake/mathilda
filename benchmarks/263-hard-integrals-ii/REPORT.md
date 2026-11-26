# Integration shootout: 263-hard-integrals-ii (2026-11-26)

Re-tabulated from the run of 2026-11-26: 6 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

1764 cases, timeout 30 s per run, wall 58.0 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (496) | A algebraic (473) | M mixed (radical on top) (399) | B radical at the bottom (396) | all |
|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (495) | 100% (471) | 100% (399) | 100% (396) | 100% (1761) |
| Mathilda:ParallelMixedTower | 99% (493) | 85% (402) | 0% (0) | 98% (390) | 73% (1285) |
| Mathilda:RischTranscendental | 100% (494) | · | · | · | 100% (494) |
| Mathilda:Trager | · | 94% (446) | · | · | 94% (446) |
| Mathilda:Bronstein | · | · | 100% (399) | 100% (396) | 100% (795) |
| FriCAS | 100% (496) | 79% (374) | 0% (0) | 97% (385) | 71% (1255) |
| Axiom | 100% (496) | 64% (301) | 0% (0) | 92% (366) | 66% (1163) |
| Maxima | 97% (482) | 71% (334) | 0% (0) | 59% (234) | 60% (1050) |
| sympy | 89% (443) | 69% (328) | 0% (0) | 28% (112) | 50% (883) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1761 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 1285 | 0 | 0 | 406 | 73 | 0 | 0 | 0 | 0.16 |
| Mathilda:RischTranscendental | 494 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0.02 |
| Mathilda:Trager | 446 | 0 | 0 | 4 | 23 | 0 | 0 | 0 | 0.64 |
| Mathilda:Bronstein | 795 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.03 |
| FriCAS | 1255 | 0 | 0 | 11 | 80 | 418 | 0 | 0 | 0.14 |
| Axiom | 1163 | 0 | 0 | 23 | 88 | 486 | 0 | 4 | 0.12 |
| Maxima | 1050 | 0 | 0 | 688 | 0 | 7 | 0 | 19 | 0.03 |
| sympy | 551 | 332 | 0 | 538 | 308 | 0 | 7 | 28 | 1.44 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1761 | 1761 | 0 | 0 | 0 | 0 | 1761 | 0 |
| Mathilda:ParallelMixedTower | 1285 | 1285 | 0 | 0 | 0 | 0 | 1285 | 0 |
| Mathilda:RischTranscendental | 494 | 494 | 0 | 0 | 0 | 0 | 494 | 0 |
| Mathilda:Trager | 446 | 446 | 0 | 0 | 0 | 0 | 446 | 0 |
| Mathilda:Bronstein | 795 | 795 | 0 | 0 | 0 | 0 | 795 | 0 |
| FriCAS | 1255 | 1255 | 0 | 0 | 0 | 0 | 1248 | 7 |
| Axiom | 1167 | 1163 | 4 | 0 | 0 | 0 | 1163 | 0 |
| Maxima | 1069 | 1050 | 19 | 0 | 0 | 0 | 1050 | 15 |
| sympy | 918 | 883 | 28 | 1 | 6 | 0 | 706 | 177 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| T | cT-mathilda | 490 | 490 | 490 | 490 | 476 | 438 |
| T | cT-fricas | 6 | 6 | 6 | 6 | 6 | 5 |
| A | cA-mathilda | 444 | 444 | 348 | 279 | 327 | 321 |
| A | cA-fricas | 23 | 21 | 23 | 21 | 2 | 2 |
| A | cA-maxima | 5 | 5 | 2 | 0 | 5 | 5 |
| A | cA-axiom | 1 | 1 | 1 | 1 | 0 | 0 |
| M | cM-mathilda | 399 | 399 | 0 | 0 | 0 | 0 |
| B | cB-mathilda | 396 | 396 | 385 | 366 | 234 | 112 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (1)

- `cT-312-33` AUTO-GAP Mathilda:ParallelMixedTower: `x*(-(((-2 + x^3)*Log[-2 + x^3])/(1 + 2*x^3 + x^6)) - (2*Log[2 + x^3])/(3 + x^2)^2)`

## Only one system closes it

- Mathilda: 428

## Wrong answers from other systems (their bugs)

- sympy `cT-301-6`: `-(((-3 + x^2)*Log[x])/(3 + x^2)^2) - (x*(-4 + x^3)*Log[1 + x + x^2])/(2 + x^3)^2`
- sympy `cT-301-9`: `(x*(-6 + x^3)*Log[-1 + x^2])/(3 + x^3)^2`
- sympy `cT-301-51`: `(-2*Log[3 + x])/(3 + x)^3 - (x*(-4 + x^3)*Log[-2 + x^3])/(2 + x^3)^2`
- sympy `cT-302-70`: `-(Log[-3 + x]/(-3 + x)^2) - ((-10 - 12*x - 3*x^2 + 2*x^3)*Log[2 + x^3])/((3 + x)^2*(2 + x^2)^2)`
- sympy `cT-303-49`: `x*((-2*Log[x])/(3 + x^2)^2 - ((-4 + x^3)*Log[3 + x^2])/(2 + x^3)^2)`
- Maxima `cT-304-3`: `(E^x*(-2 + x))/(-1 + x)^2 - (2*x*Log[(3 + x)/(-2 + x)])/(2 + x^2)^2`
- sympy `cT-305-23`: `(x*(-4 + x^3)*Log[-2 + x^3])/(2 + x^3)^2`
- sympy `cT-306-56`: `((10 + 12*x + 3*x^2 - 2*x^3)*Log[2 + x^3])/((3 + x)^2*(2 + x^2)^2)`
- sympy `cT-307-24`: `(E^x*(-4 + x))/(-3 + x)^2 - (2*(-2 + x)*(1 + x + x^2)*Log[2 + x^3])/((1 + x)^2*(2 + x^2)^2)`
- sympy `cT-308-28`: `-(((-3 + x^2)*Log[3 + x^2])/(3 + x^2)^2) - Log[-2 + x^3]/(-3 + x)^2`
- sympy `cT-308-53`: `Log[-2 + x^3]/(2 + x)^2`
- sympy `cT-309-13`: `((1 - 4*x - 5*x^2 - 2*x^3)*Log[-2 + x^3])/((2 + x)^2*(1 + x^2)^2)`
- Maxima `cT-309-46`: `x*((E^x*(4 + 2*x - x^3 + x^4))/(2 + x^3)^2 - (2*Log[-1 + x^2])/(1 + x^2)^2)`
- sympy `cT-310-9`: `(x*(2 - x^3)*Log[2 + x^3])/(1 + x^3)^2`
- Maxima `cT-310-37`: `x*(-(((-2 + x^3)*Log[-2 + x^2])/(1 + 2*x^3 + x^6)) - ((-4 + x^3)*Log[-1 + x^2])/(2 + x^3)^2)`
- sympy `cT-311-44`: `-(((-2 + 2*x + x^2)*Log[1 + (-2 + x)^(-1)])/(3 + x + x^2)^2) - (x*(-4 + x^3)*Log[-2 + x])/(2 + x^3)^2`
- sympy `cT-312-50`: `(x*(-4 + x^3)*Log[-2 + x])/(2 + x^3)^2`
- sympy `cT-313-18`: `(x*(2 - x^3)*Log[-2 + x^3])/(1 + x^3)^2`
- sympy `cT-313-26`: `(-3*x^2*Log[-2 + x])/(2 + x^3)^2 - (2*Log[3 + x + x^2])/(-2 + x)^3`
- sympy `cT-313-28`: `((1 - x + x^2 + x^3)*Log[-2 + x^3])/((-1 + x)^2*(1 + x^2)^2)`
- sympy `cT-313-43`: `-(((1 + x^2)*Log[x])/((-1 + x)^2*(1 + x)^2)) - (3*x^2*Log[(-2 + x)/(2 + x)])/(2 + x^3)^2`
- sympy `cT-314-22`: `(-2*x*Log[(3 + x)/(-2 + x)])/(2 + x^2)^2 - Log[-2 + x^3]/(-2 + x)^2`
- Maxima `cT-314-46`: `(E^x*(-4 + x))/(-2 + x)^3 - (2*(-4 - 4*x - 2*x^2 + x^3)*Log[-1 + x^2])/((2 + x)^2*(2 + x^2)^2)`
- Maxima `cT-315-42`: `(-E^x - 3*E^x*x^2 + E^x*x^3 - 3*x^2*Log[(-1 + x)/(-3 + x)])/(-1 + x^3)^2`
- sympy `cT-316-18`: `(E^x*(1 + x))/(3 + x)^3 - (2*(4 - 4*x + 2*x^2 + x^3)*Log[2 + x^3])/((-2 + x)^2*(2 + x^2)^2)`
- sympy `cT-316-56`: `(E^x*(-4 + x))/(-2 + x)^3 - (x*(-6 + x^3)*Log[-1 + x^2])/(3 + x^3)^2`
- Maxima `cT-317-18`: `x*((-2*Log[-2 + x])/(3 + x^2)^2 - ((-2 + x^3)*Log[(2 + x)/(-1 + x)])/(1 + 2*x^3 + x^6))`
- sympy `cT-317-45`: `-(((1 + 4*x + x^2)*Log[x])/(1 + x + x^2)^2) + ((1 + x)^(-2) + 2/(1 + x^2)^2 - (1 + x^2)^(-1))*Log[-2 + x^3]`
- sympy `cT-318-34`: `(-2*Log[3 + x + x^2])/(2 + x)^3 - ((-2 + 2*x + x^2)*Log[2 + x^3])/(3 + x + x^2)^2`
- sympy `cT-318-40`: `(x^2*Log[1 + x^2])/(2 + x^3)^2`
- sympy `cT-319-17`: `(-(-1 + x)^(-2) - 2/(1 + x^2)^2 + (1 + x^2)^(-1))*Log[x] - ((1 + x^2)*Log[-2 + x^3])/((-1 + x)^2*(1 + x)^2)`
- sympy `cT-320-38`: `x*((-3*x*Log[x])/(2 + x^3)^2 - (2*Log[-2 + x^3])/(2 + x^2)^2)`
- sympy `cT-320-62`: `-((x*(-6 + x^3)*Log[x])/(3 + x^3)^2) - (2*Log[-1 + x^2])/(3 + x)^3`
- sympy `cT-322-10`: `(-2*Log[x])/(-3 + x)^3 - ((-2 + x^2)*Log[-2 + x^3])/(2 + x^2)^2`
- Axiom `cA-313-80`: `(x^7*(3 + x^3)^(4/3))/(1 + x^6)`
- Axiom `cA-313-88`: `(x^7*(-3 + x^3)^(4/3))/(1 + x^6)`
- Axiom `cA-332-118`: `(x^4*(1 + x^3)^(4/3))/(2 + x^6)`
- Axiom `cA-337-73`: `x^7/((2 + x^3)^(2/3)*(2 + x^6))`
- Maxima `cB-304-10`: `(-3*(2 + x + x^3) + x*Log[2 + x^(-1)] + 2*x^2*Log[2 + x^(-1)] + 3*x^3*Log[2 + x^(-1)] + 6*x^4*Log[2 + x^(-1)])/((x + 2*x^2)*(2 + x + x^3)^(2/3))`
- Maxima `cB-304-35`: `(-1 - x^4 + x^4*Log[2 + x^(-1)] + 2*x^5*Log[2 + x^(-1)])/(x*(1 + 2*x)*(1 + x^4)^(3/4))`
- Maxima `cB-305-23`: `(E^x*(-1 + 2*x)*(-2 + x^2))/Sqrt[1 - 2*x^2 + x^3]`
- Maxima `cB-306-12`: `(-2 - x^4 + x^4*Log[2 + x^(-1)] + 2*x^5*Log[2 + x^(-1)])/(x*(1 + 2*x)*(2 + x^4)^(3/4))`
- Maxima `cB-306-16`: `(E^Sqrt[-1 + x^4]*(2*x^3 + 2*x^4 - Sqrt[-1 + x^4]))/((1 + x)^2*Sqrt[-1 + x^4])`
- Maxima `cB-317-11`: `(-4 - 2*x^3 + 3*x^3*Log[1 + x^(-1)] + 3*x^4*Log[1 + x^(-1)])/(x*(1 + x)*Sqrt[2 + x^3])`
- Maxima `cB-317-37`: `(Sqrt[1 - 2*x^2 + x^3]*(-2 + 4*x^2 - 2*x^3 + 8*x*Log[2 + x] - 2*x^2*Log[2 + x] - 3*x^3*Log[2 + x]))/(2 + x - 4*x^2 + x^4)`
- Maxima `cB-317-42`: `(Sqrt[1 - 2*x^2 + x^3]*(4 - 6*x^2 + 4*x^3 - 4*x^4 + 2*x^5 + 4*Log[x] - 18*x^2*Log[x] + 10*x^3*Log[x] + x^5*Log[x]))/((-1 + x)*(2 + x^2)^2*(-1 - x + x^2))`
- Maxima `cB-318-8`: `(-4 - 4*x^2 - 2*x^3 + 4*x^2*Log[2 + x^(-1)] + 11*x^3*Log[2 + x^(-1)] + 6*x^4*Log[2 + x^(-1)])/(x*(1 + 2*x)*Sqrt[2 + 2*x^2 + x^3])`
- Maxima `cB-319-12`: `(Sqrt[1 - 2*x^2 + x^3]*(2 - 4*x^2 + 2*x^3 - 8*x*Log[2 + x] + 2*x^2*Log[2 + x] + 3*x^3*Log[2 + x]))/(2 + x - 4*x^2 + x^4)`
- Maxima `cB-321-13`: `-(Sqrt[2 + x + x^4]/(x + 2*x^2)) + ((Sqrt[2 + x + x^4] + 4*x^3*Sqrt[2 + x + x^4])*Log[2 + x^(-1)])/(4 + 2*x + 2*x^4)`
- Maxima `cB-323-2`: `(E^(2*x)*(-5 - 2*x + 4*x^3 + 7*x^4 + 2*x^5))/Sqrt[-1 + x^4]`
- Maxima `cB-329-7`: `(E^x*x*(-2 - x + 4*x^4 + x^5))/Sqrt[-1 + x^4]`
