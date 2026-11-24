# Integration shootout: 261-compact-integrals-x4 (2026-11-24)

1550 cases, timeout 30 s per run, wall 1.3 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | R rational (50) | T transcendental (250) | A algebraic (250) | M mixed (radical on top) (550) | B radical at the bottom (450) | all |
|---|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (50) | 100% (250) | 100% (250) | 100% (549) | 100% (450) | 100% (1549) |
| Mathilda:ParallelMixedTower | 0% (0) | 100% (250) | 96% (240) | 90% (497) | 86% (386) | 89% (1373) |
| Mathilda:RischTranscendental | 98% (49) | 100% (250) | · | · | · | 100% (299) |
| Mathilda:Trager | · | · | 100% (250) | · | · | 100% (250) |
| Mathilda:Bronstein | · | · | · | 100% (548) | 100% (450) | 100% (998) |
| FriCAS | 100% (50) | 100% (249) | 97% (243) | 87% (481) | 100% (449) | 95% (1472) |
| Axiom | 98% (49) | 96% (239) | 82% (206) | 70% (383) | 98% (442) | 85% (1319) |
| Maxima | 60% (30) | 43% (108) | 50% (124) | 39% (214) | 94% (421) | 58% (897) |
| sympy | 60% (30) | 69% (173) | 10% (26) | 4% (24) | 24% (107) | 23% (360) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1549 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.12 |
| Mathilda:ParallelMixedTower | 1373 | 0 | 0 | 61 | 116 | 0 | 0 | 0 | 0.30 |
| Mathilda:RischTranscendental | 299 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0.03 |
| Mathilda:Trager | 250 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.21 |
| Mathilda:Bronstein | 998 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0.19 |
| FriCAS | 1472 | 0 | 0 | 0 | 32 | 43 | 3 | 0 | 0.17 |
| Axiom | 1319 | 0 | 0 | 1 | 13 | 61 | 0 | 156 | 0.16 |
| Maxima | 897 | 0 | 0 | 653 | 0 | 0 | 0 | 0 | 0.04 |
| sympy | 319 | 41 | 0 | 536 | 555 | 62 | 26 | 11 | 3.60 |

## Cross-check in Mathilda (answers that failed their own check)

Own check: `bad` = a sample point failed; `na` = no point could be evaluated. Verdict of PossibleZeroQ[D[F, x] - f] in Mathilda.

| run | own | n | zero (verified) | nonzero (wrong) | undecided | untranslatable |
|---|---|---|---|---|---|---|
| FriCAS | bad | 60 | 58 | 0 | 2 | 0 |
| FriCAS | na | 40 | 39 | 0 | 1 | 0 |
| Axiom | bad | 163 | 17 | 146 | 0 | 0 |
| Axiom | na | 10 | 0 | 10 | 0 | 0 |
| sympy | bad | 37 | 0 | 11 | 4 | 22 |
| sympy | na | 8 | 8 | 0 | 0 | 0 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS | Axiom | Maxima | sympy |
|---|---|---|---|---|---|---|---|
| R | cR-mathilda | 50 | 50 | 50 | 49 | 30 | 30 |
| T | cT-mathilda | 250 | 250 | 249 | 239 | 108 | 173 |
| A | cA-mathilda | 250 | 250 | 243 | 206 | 124 | 26 |
| M | cM-mathilda | 539 | 539 | 472 | 375 | 209 | 24 |
| M | cM-fricas | 6 | 6 | 6 | 5 | 3 | 0 |
| M | cM-maxima | 2 | 2 | 0 | 0 | 2 | 0 |
| M | cM-axiom | 3 | 3 | 3 | 3 | 0 | 0 |
| B | cB-mathilda | 450 | 450 | 449 | 442 | 421 | 107 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (1)

- `cM-102-291` AUTO-GAP Mathilda:Bronstein: `(2*(1 + E^x^2)*(1 + E^(6*x^2))^(1/6)*x)/E^x^2`

## Only one system closes it

- Mathilda: 69

## Wrong answers from other systems (their bugs)

- sympy `cR-101-54`: `(-2 + x)/(1 + x^6)`
- sympy `cR-101-126`: `(2 + x^6)^(-1)`
- sympy `cR-101-401`: `(1 + 2*x)/(1 + x^6)`
- sympy `cR-101-513`: `(-2 + x + 2*x^2)/(1 + x^6)`
- sympy `cR-101-699`: `(-2 - x + 2*x^2)/(1 + x^6)`
- Axiom `cR-101-808`: `(1 + 2*x^2)/(1 + x^6)`
- sympy `cR-101-876`: `(2 + x)/(1 + x^6)`
- sympy `cR-101-1030`: `(1 + x)/(1 + x^6)`
- sympy `cR-101-1035`: `(-1 + x)/(1 + x^6)`
- sympy `cR-101-1065`: `(1 - 2*x + 2*x^2)/(1 + x^6)`
- sympy `cR-101-1108`: `(-1 - 2*x + 2*x^2)/(1 + x^6)`
- sympy `cR-101-1283`: `(2*(-1 + x))/(1 + x^6)`
- Axiom `cT-101-1752`: `2/(x*(2 - Log[x]^2 + 2*Log[x]^4))`
- Axiom `cT-102-5378`: `(2*E^x)/(2 - 2*E^(2*x) + E^(4*x))`
- Axiom `cT-103-5964`: `(2*E^(2*x))/(1 + E^x + 2*E^(2*x) - E^(3*x) + E^(4*x))`
- Axiom `cT-105-2055`: `E^x/((1 + E^x)*(2 - Log[1 + E^x]^2 + Log[1 + E^x]^4))`
- Axiom `cT-107-1063`: `E^x/((1 + E^x)*(1 - 2*Log[1 + E^x]^2 + 2*Log[1 + E^x]^4))`
- Axiom `cT-108-71`: `(2*E^x)/(1 - 2*E^(2*x) + 2*E^(4*x))`
- Axiom `cT-109-3488`: `E^(3*x)/(-2 + E^(2*x) - 2*E^(4*x))`
- Axiom `cT-109-4032`: `(2*E^(2*x))/(2 - E^(4*x) + E^(8*x))`
- Axiom `cT-111-2874`: `(-2*E^x^2*x)/(-2 + E^(2*x^2) - 2*E^(4*x^2))`
- Axiom `cT-111-4565`: `(4*E^x^2*x)/(2 - E^(2*x^2) + E^(4*x^2))`
- Axiom `cT-112-1638`: `(2*E^(2*x))/(2 - 2*E^(4*x) + E^(8*x))`
- Axiom `cA-101-189`: `(2 + x)/((1 + x)^(3/2)*(2 - x + 2*x^2))`
- Axiom `cA-101-303`: `Sqrt[1 + x]/(2 + x^2)`
- Axiom `cA-101-357`: `x^3/((2 + x^2)^(3/2)*(2 - x + 2*x^2))`
- Axiom `cA-102-2`: `Sqrt[1 + x]/(3 + x^2)`
- Axiom `cA-102-25`: `(2 + x)/((3 + x)^(3/2)*(2 + x + x^2))`
- Axiom `cA-102-160`: `((-1 + x)*(2 + x + x^2)^(3/2))/(2 - x + x^2)`
- Axiom `cA-102-296`: `x^2/((2 + x + x^2)^(3/2)*(2 + x + 2*x^2))`
- Axiom `cA-102-411`: `1/(Sqrt[2 + x]*(1 + x^2))`
- Axiom `cA-102-437`: `((-1 + x)*Sqrt[3 + x^2])/(1 + 2*x^2)`
- Axiom `cB-102-985`: `(2*x)/((1 + x^2)*(2 + x + x^2)^(3/2))`
- Axiom `cA-103-105`: `1/((3 + x)^(3/2)*(2 + x^2))`
- Axiom `cA-103-237`: `(Sqrt[1 + x]*(2 + x))/(1 + x^2 + 2*x^3)`
- Axiom `cA-103-358`: `((-1 + x)*Sqrt[2 + x])/(2 + 2*x^2)`
- Axiom `cA-103-479`: `((-2 + x)*Sqrt[2 + x + x^2])/(1 + 2*x + 2*x^2)`
- Axiom `cA-103-531`: `x/((2 + x + x^2)^(3/2)*(2 - x + x^2 + x^3))`
- Axiom `cA-104-12`: `((-1 + x)*Sqrt[2 + x + x^2])/(2 + x^2)`
- Axiom `cA-104-312`: `((2 + x)*(3 + x^2)^(3/2))/(1 + x + x^2)`
- Axiom `cA-104-544`: `1/(Sqrt[3 + x]*(2 + x^2))`
- Axiom `cA-105-83`: `((1 + x)*(1 + x + x^2)^(3/2))/(2 - x + x^2)`
- Axiom `cA-105-386`: `1/((3 + x)^(3/2)*(1 + x^2))`
- Axiom `cA-105-443`: `x/((2 + x)^(1/6)*(2 - x^2 + x^3))`
- Axiom `cA-105-457`: `x/((3 + x^2)^(3/2)*(1 + x + x^2))`
- Axiom `cA-105-532`: `((-2 + x)*Sqrt[2 + x])/(1 - 2*x + 2*x^2)`
- Axiom `cA-105-609`: `(2 + x)/((2 + x^2)^(3/2)*(1 - x + x^2))`
- Axiom `cA-106-52`: `x^2/((2 + x + x^2)^(3/2)*(2 - x^2 + x^3))`
- Axiom `cA-106-64`: `1/(Sqrt[1 + x]*(1 + x^2))`
- Axiom `cA-106-116`: `((-2 + x)*(3 + x)^(3/2))/(2*(1 + x + x^2))`
- Axiom `cA-106-184`: `x^2/((1 + x)*(3 + x^2)^(3/2)*(2 - 3*x + 2*x^2))`
- Axiom `cA-106-197`: `Sqrt[1 + x]/(2 - x + 2*x^2)`
- Axiom `cA-106-218`: `Sqrt[1 + x]/(1 + x^2)`
- Axiom `cA-106-230`: `1/(Sqrt[1 + x]*(3 + x^2))`
- Axiom `cA-107-65`: `x/((1 + x)^(3/2)*(1 - x + 2*x^2))`
- Axiom `cA-107-143`: `Sqrt[3 + x]/(2 + x^2)`
- Axiom `cA-107-244`: `((2 + x)*(3 + x)^(3/2))/(1 - x + 2*x^2)`
- Axiom `cA-107-281`: `(2 + x)/(Sqrt[1 + x]*(2 + x + 2*x^2))`
- Axiom `cA-108-62`: `((-2 + x)*Sqrt[1 + x + x^2])/(2 - 2*x + x^2)`
- Axiom `cA-108-148`: `((-1 + x)*(2 + x + x^2)^(3/2))/(2 - x + 2*x^2)`
- Axiom `cM-101-172`: `(2*(2 + Log[1 + 2*x]))/((1 + 2*x)*(1 + Log[1 + 2*x] + Log[1 + 2*x]^2)^(3/2)*(1 - 2*Log[1 + 2*x] + 2*Log[1 + 2*x]^2))`
- Axiom `cM-101-618`: `(2*(2 + Log[1 + 2*x])*(1 + Log[1 + 2*x]^2)^(3/2))/((1 + 2*x)*(2 - Log[1 + 2*x] + 2*Log[1 + 2*x]^2))`
- Axiom `cM-101-653`: `(2*E^(2*x))/(Sqrt[1 + E^(2*x)]*(2 + E^(4*x)))`
- Axiom `cM-101-676`: `(2*Sqrt[3 + Log[1 + 2*x]])/((1 + 2*x)*(3 + Log[1 + 2*x]^2))`
- Axiom `cM-101-681`: `(2*E^(2*x)*Sqrt[2 + E^(2*x)])/(2 + E^(4*x))`
- Axiom `cM-101-691`: `Sqrt[2 + Log[x]]/(x*(3 + Log[x]^2))`
- Axiom `cM-101-721`: `E^x/((1 + E^x)*(2 + Log[1 + E^x])^(1/4)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-101-835`: `(E^x*(2 + Log[1 + E^x])^(1/3))/((1 + E^x)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-101-848`: `(E^x*(1 + E^x + E^(2*x))^(3/2))/(2 + E^(2*x))`
- Axiom `cM-101-888`: `(2*E^(2*x)*(1 + E^(2*x))^(1/4))/(3 + E^(4*x))`
- Axiom `cM-102-132`: `E^x/((1 + E^x)*(1 + Log[1 + E^x])^(1/4)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-102-446`: `((-2 + Log[x])*(1 + Log[x]^2)^(3/2))/(x - 2*x*Log[x] + 2*x*Log[x]^2)`
- Axiom `cM-102-452`: `(2*(3 + Log[1 + 2*x])^(3/2))/((1 + 2*x)*(2 + Log[1 + 2*x]^2))`
- Axiom `cM-102-567`: `(E^x*(3 + E^x)^(3/2))/(1 + E^(2*x))`
- Axiom `cM-102-778`: `(1 + Log[1 + x])/(2*(1 + x)*(1 + Log[1 + x]^2)*(2 + Log[1 + x] + Log[1 + x]^2)^(3/2))`
- Axiom `cM-102-957`: `(2*E^x^2*(2 + E^x^2)^(1/4)*x)/(1 + E^(2*x^2))`
- Axiom `cM-102-978`: `-((E^x*(1 + E^(-2*x) + E^(-x))^(3/2))/(1 + 2*E^(2*x)))`
- Axiom `cM-103-596`: `(E^x*(1 + Log[1 + E^x] + Log[1 + E^x]^2)^(3/2))/((1 + E^x)*(3 + Log[1 + E^x]^2))`
- Axiom `cM-103-686`: `1/(x*Sqrt[2 + Log[x]]*(3 + Log[x]^2))`
- Axiom `cM-103-792`: `(E^x*(1 + Log[1 + E^x])*(2 + Log[1 + E^x])^(1/6))/(2*(1 + E^x)*(1 + Log[1 + E^x] + Log[1 + E^x]^2))`
- Axiom `cM-103-827`: `-((E^x*(1 + E^x))/((-2 + E^x - 2*E^(2*x))*Sqrt[2 + E^x + E^(2*x)]))`
- Axiom `cM-103-884`: `1/(x*Sqrt[1 + Log[x]]*(3 + Log[x]^2))`
- Axiom `cM-103-1053`: `(Sqrt[2 + E^(-x)]*(-1 - 2*E^x))/(1 + E^(2*x))`
- Axiom `cM-104-153`: `(2*E^x^2*x)/((3 + E^x^2)^(3/2)*(2 + E^(2*x^2)))`
- Axiom `cM-104-381`: `((-1 + Log[x])*(1 + Log[x])^(3/2))/(x*(2 + Log[x] + 2*Log[x]^2))`
- Axiom `cM-104-459`: `(2 + Log[1 + x])/((1 + x)*Sqrt[1 + Log[1 + x]]*(2 + Log[1 + x] + 2*Log[1 + x]^2))`
- Axiom `cM-104-585`: `(E^x*(3 + Log[1 + E^x])^(1/6))/((1 + E^x)*(3 + Log[1 + E^x]^2))`
- Axiom `cM-104-676`: `(E^x*(1 + Log[1 + E^x] + Log[1 + E^x]^2)^(3/2))/((1 + E^x)*(2 + Log[1 + E^x]^2))`
- Axiom `cM-104-833`: `-((E^x*(3 + E^(-x))^(1/4))/(1 + E^(2*x)))`
- Axiom `cM-104-1155`: `(2*Sqrt[3 + Log[1 + 2*x]])/((1 + 2*x)*(2 + Log[1 + 2*x]^2))`
- Axiom `cM-105-40`: `(E^x*Sqrt[2 + E^x])/(2 + E^(2*x))`
- Axiom `cM-105-67`: `1/((1 + x)*Sqrt[1 + Log[1 + x]]*(2 + Log[1 + x]^2))`
- Axiom `cM-105-69`: `(E^x*(1 + Log[1 + E^x])*Sqrt[3 + Log[1 + E^x]^2])/((1 + E^x)*(1 - Log[1 + E^x] + 2*Log[1 + E^x]^2))`
- Axiom `cM-105-101`: `((-2 + Log[x])*(2 + Log[x] + Log[x]^2)^(3/2))/(x*(2 - Log[x] + 2*Log[x]^2))`
- Axiom `cM-105-193`: `(E^x*(-2 + Log[1 + E^x])*(1 + Log[1 + E^x])^(1/4))/((1 + E^x)*(2 + Log[1 + E^x]^2))`
- Axiom `cM-105-505`: `-((E^x*Sqrt[3 + E^(-x)])/(1 + 2*E^(2*x)))`
- Axiom `cM-105-618`: `(2 + Log[1 + x])^(3/2)/((1 + x)*(1 + Log[1 + x]^2))`
- Axiom `cM-105-632`: `(2*E^(2*x)*(2 + E^(2*x))*(1 + E^(2*x) + E^(4*x))^(3/2))/(2 + E^(2*x) + 2*E^(4*x))`
- Axiom `cM-106-196`: `(2*E^x^2*(3 + E^x^2)^(3/2)*x)/(2 + E^(2*x^2))`
- Axiom `cM-106-267`: `1/(Sqrt[3 + Log[x]]*(x + x*Log[x]^2))`
- Axiom `cM-106-322`: `(E^x*(2 + Log[1 + E^x])^(7/6))/((1 + E^x)*(2 + Log[1 + E^x] + 2*Log[1 + E^x]^2))`
- Axiom `cM-106-681`: `(-2*E^(2*x)*(1 + E^(2*x))^(5/4))/(-1 + E^(2*x) - 2*E^(4*x))`
- Axiom `cM-106-699`: `(2*x*Sqrt[2 + Log[1 + x^2]])/((1 + x^2)*(2 + Log[1 + x^2]^2))`
- Axiom `cM-106-808`: `-((E^x*(1 + E^(-x))^(1/4))/(1 + 3*E^(2*x)))`
- Axiom `cM-106-951`: `-((E^x*(3 + E^(-x))^(2/3))/(1 + 2*E^(2*x)))`
- Axiom `cM-106-1019`: `((-1 + Log[x])*Sqrt[2 + Log[x]])/(x + x*Log[x] + 2*x*Log[x]^2)`
- Axiom `cM-106-1248`: `Sqrt[3 + Log[x]]/(x + x*Log[x]^2)`
- Axiom `cM-107-334`: `(E^x*(2 + E^x)^(3/2))/(1 + E^(2*x))`
- Axiom `cM-107-613`: `1/((1 + x)*(3 + Log[1 + x]^2)*(1 + Log[1 + x] + Log[1 + x]^2)^(3/2))`
- Axiom `cM-107-665`: `(E^x*Sqrt[1 + E^x])/(1 + E^(2*x))`
- Axiom `cM-107-700`: `(E^x*(2 + Log[1 + E^x])^(1/6))/((1 + E^x)*(3 + Log[1 + E^x]^2))`
- Axiom `cM-107-738`: `-((E^x*(1 + E^(-x))^(2/3))/(1 + 3*E^(2*x)))`
- Axiom `cM-107-863`: `1/(Sqrt[2 + Log[x]]*(x + x*Log[x]^2))`
- Axiom `cM-108-133`: `(2*(-1 + Log[1 + 2*x])*Sqrt[3 + Log[1 + 2*x]^2])/((1 + 2*x)*(1 + 2*Log[1 + 2*x]^2))`
- Axiom `cM-108-465`: `(-2*E^(2*x)*(-1 + E^x)*(1 + E^x)*(2 + E^(4*x))^(3/2))/(-1 + E^(2*x) - 2*E^(4*x))`
- Axiom `cM-108-502`: `(E^x*(3 + Log[1 + E^x])^(1/6))/((1 + E^x)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-108-618`: `((-1 + Log[x])*(1 + Log[x] + Log[x]^2)^(3/2))/(2*x*(1 - Log[x] + Log[x]^2))`
- Axiom `cM-108-662`: `1/(x*Sqrt[2 + Log[x]]*(2 + Log[x]^2))`
- Axiom `cM-108-674`: `Sqrt[1 + Log[1 + x]]/((1 + x)*(2 + Log[1 + x]^2))`
- Axiom `cM-108-1078`: `(E^x*(2 + Log[1 + E^x])^(2/3))/((1 + E^x)*(2 + Log[1 + E^x]^2))`
- Axiom `cM-109-144`: `E^x/(Sqrt[3 + E^(-x)]*(-1 - 2*E^(2*x)))`
- Axiom `cM-109-200`: `(Sqrt[3 + E^(-x)]*(-1 - 3*E^x))/(1 + 2*E^(2*x))`
- Axiom `cM-109-230`: `(-2*E^x^2*(-2 + E^x^2)*x)/(Sqrt[3 + E^x^2]*(-1 + 2*E^x^2 - 2*E^(2*x^2)))`
- Axiom `cM-109-333`: `(E^x*(1 + Log[1 + E^x])^(4/3))/((1 + E^x)*(1 + 2*Log[1 + E^x]^2))`
- Axiom `cM-109-683`: `(E^x*(2 + Log[1 + E^x])^(1/4))/((1 + E^x)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-109-850`: `(E^x*(2 + E^x)^(1/4))/(3 + E^(2*x))`
- Axiom `cM-109-882`: `(2*x)/((1 + x^2)*(3 + Log[1 + x^2]^2)*(1 + Log[1 + x^2] + Log[1 + x^2]^2)^(3/2))`
- Axiom `cM-109-907`: `(E^x*Sqrt[1 + E^x])/(2 + E^(2*x))`
- Axiom `cM-110-99`: `(2*E^x^2*Sqrt[3 + E^x^2]*x)/(3 + E^(2*x^2))`
- Axiom `cM-110-158`: `(E^x*(2 + Log[1 + E^x])^(3/2))/((1 + E^x)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-110-359`: `((2 + E^(2*x))^(3/2)*(-E^(2*x) + E^(4*x)))/(1 - E^(2*x) + E^(4*x))`
- Axiom `cM-110-746`: `(1 + E^x)/((2 + E^(-x))^(1/4)*(-2 + 2*E^x - E^(2*x)))`
- Axiom `cM-110-796`: `(E^x*Sqrt[3 + E^x])/(2 + E^(2*x))`
- Axiom `cM-110-1519`: `Sqrt[3 + Log[x]]/(x*(3 + Log[x]^2))`
- Axiom `cM-110-1535`: `-((Sqrt[3 + E^(-2*x)]*(1 + 2*E^x)*(1 + 3*E^(2*x)))/(2*E^(2*x) + 2*E^(4*x)))`
- Axiom `cM-111-42`: `(E^x*(-1 + E^x))/((2 + E^x + E^(2*x))^(3/2)*(2 + 2*E^x + E^(2*x)))`
- Axiom `cM-111-196`: `(E^x*Sqrt[2 + E^x])/(1 + E^(2*x))`
- Axiom `cM-111-214`: `(2*(1 + Log[1 + 2*x])^(3/2)*(2 + Log[1 + 2*x]))/((1 + 2*x)*(1 - 2*Log[1 + 2*x] + 2*Log[1 + 2*x]^2))`
- Axiom `cM-111-232`: `(2*E^x^2*(2 + E^x^2)*x)/((2 + E^x^2 + E^(2*x^2))^(3/2)*(2 + E^x^2 + 2*E^(2*x^2)))`
- Axiom `cM-111-383`: `(E^x*(1 + Log[1 + E^x])^(1/4))/((1 + E^x)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-111-407`: `-((E^x*(-2 + E^x))/(Sqrt[2 + E^x]*(-1 + E^x - E^(2*x))))`
- Axiom `cM-111-442`: `2/((1 + 2*x)*Sqrt[1 + Log[1 + 2*x]]*(1 - 2*Log[1 + 2*x] + 2*Log[1 + 2*x]^2))`
- Axiom `cM-111-643`: `(2*x*Sqrt[3 + Log[1 + x^2]])/((1 + x^2)*(3 + Log[1 + x^2]^2))`
- Axiom `cM-111-666`: `(2*E^(2*x)*(3 + E^(2*x))^(3/2))/(2 + E^(4*x))`
- Axiom `cM-111-853`: `Sqrt[1 + Log[x]]/(x + x*Log[x]^2)`
- Axiom `cM-111-969`: `(2*E^x^2*(1 + E^x^2 + E^(2*x^2))^(3/2)*x)/(3 + E^(2*x^2))`
- Axiom `cM-111-1011`: `(2*E^x^2*Sqrt[3 + E^x^2]*x)/(2 + E^(2*x^2))`
- Axiom `cM-112-346`: `(2*E^(2*x))/(Sqrt[3 + E^(2*x)]*(3 + E^(4*x)))`
- Axiom `cM-112-744`: `(2*E^x^2*x)/((3 + E^x^2)^(3/2)*(1 + E^(2*x^2)))`
- Axiom `cM-112-989`: `(E^x*(2 + Log[1 + E^x])^(1/4))/((1 + E^x)*(2 + Log[1 + E^x]^2))`
- Axiom `cM-112-1070`: `(1 + Log[x])/(x*(2 + Log[x]^2)^(3/2)*(2 + Log[x] + Log[x]^2))`
- Axiom `cM-112-1169`: `E^x/((1 + E^x)*(1 + Log[1 + E^x])^(1/4)*(2 + Log[1 + E^x]^2))`
- Axiom `cM-112-1221`: `E^x/((1 + E^x)*(3 + Log[1 + E^x])^(3/2)*(1 + Log[1 + E^x]^2))`
- Axiom `cM-113-83`: `(2*E^(2*x)*(-2 + E^(2*x)))/(Sqrt[2 + E^(2*x) + E^(4*x)]*(2 + E^(2*x) + 2*E^(4*x)))`
- Axiom `cM-113-381`: `(2*Sqrt[2 + Log[1 + 2*x]])/((1 + 2*x)*(2 + Log[1 + 2*x]^2))`
- Axiom `cM-113-414`: `(E^x*(1 + Log[1 + E^x])^(1/6)*(2 + Log[1 + E^x]))/((1 + E^x)*(-2 + Log[1 + E^x]^2))`
- Axiom `cM-113-557`: `(3 + Log[1 + x])^(3/2)/((1 + x)*(2 + Log[1 + x]^2))`
- Axiom `cM-113-617`: `1/((1 + x)*Sqrt[3 + Log[1 + x]]*(1 + Log[1 + x]^2))`
- Axiom `cM-113-628`: `(2*E^(2*x))/(Sqrt[2 + E^(2*x)]*(3 + E^(4*x)))`
- Axiom `cM-113-782`: `(2*E^x^2*x)/(Sqrt[2 + E^x^2]*(3 + E^(2*x^2)))`
- Axiom `cM-113-935`: `2/((1 + 2*x)*(3 + Log[1 + 2*x]^2)*(1 + Log[1 + 2*x] + Log[1 + 2*x]^2)^(3/2))`
- Axiom `cB-101-284`: `(1 + x)^(3/2)*Log[1 + x^2]`
- Axiom `cB-101-903`: `Log[1 + x^2]/(2 + x)^(1/3)`
- Axiom `cB-101-2055`: `(3 + x)^(1/4)*Log[1 + x^2]`
- Axiom `cB-102-1677`: `Log[1 + x^2]/(2 + x)^(1/4)`
- Axiom `cB-102-1915`: `Sqrt[1 + x]*Log[1 + x^2]`
- Axiom `cB-103-3445`: `(1 + x)^(1/4)*Log[1 + x^2]`
