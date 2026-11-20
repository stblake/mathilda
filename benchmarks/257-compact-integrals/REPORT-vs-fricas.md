# Integration shootout: 257-compact-integrals (2026-11-20)

1000 cases, timeout 30 s per run, wall 1.9 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | R rational (100) | T transcendental (150) | A algebraic (300) | M mixed (radical on top) (250) | B radical at the bottom (200) | all |
|---|---|---|---|---|---|---|
| Mathilda:Automatic | 100% (100) | 100% (150) | 100% (300) | 100% (250) | 100% (200) | 100% (1000) |
| Mathilda:RischTranscendental | 100% (100) | 97% (146) | · | · | · | 98% (246) |
| Mathilda:Trager | · | · | 99% (297) | · | · | 99% (297) |
| Mathilda:Bronstein | · | · | · | 99% (248) | 100% (200) | 100% (448) |
| FriCAS | 100% (100) | 95% (143) | 85% (255) | 84% (211) | 95% (190) | 90% (899) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.03 |
| Mathilda:RischTranscendental | 246 | 0 | 0 | 4 | 0 | 0 | 0 | 0 | 0.01 |
| Mathilda:Trager | 297 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0.11 |
| Mathilda:Bronstein | 448 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0.15 |
| FriCAS | 899 | 0 | 0 | 0 | 4 | 36 | 22 | 39 | 0.09 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| R | cR-mathilda | 100 | 100 | 100 |
| T | cT-mathilda | 150 | 150 | 143 |
| A | cA-mathilda | 300 | 300 | 255 |
| M | cM-mathilda | 238 | 238 | 199 |
| M | cM-fricas | 12 | 12 | 12 |
| B | cB-mathilda | 200 | 200 | 190 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (0)


## Only one system closes it

- Mathilda: 101

## Wrong answers from other systems (their bugs)

- FriCAS `cA-1-47`: `(-2 + x)/((1 + x^2)^(3/2)*(1 + x + 2*x^2))`
- FriCAS `cA-1-65`: `1/(Sqrt[1 + x]*(2 + x^2))`
- FriCAS `cA-1-186`: `(2 + x)/((3 + x^2)^(3/2)*(1 + 2*x + 2*x^2))`
- FriCAS `cA-1-219`: `(1 + x + x^2)^(3/2)/(3 + x^2)`
- FriCAS `cA-2-3`: `(1 + x)/((1 - x + x^2)*Sqrt[2 + x + x^2])`
- FriCAS `cA-2-8`: `x^3/(Sqrt[2 + x]*(2 - x + 2*x^2))`
- FriCAS `cA-2-72`: `1/(Sqrt[3 + x]*(3 + x^2))`
- FriCAS `cA-2-139`: `1/((2 + x^2)*(1 + x + x^2)^(3/2))`
- FriCAS `cA-2-185`: `1/(Sqrt[2 + x]*(2 + 2*x + x^2))`
- FriCAS `cA-2-228`: `Sqrt[2 + x]/(3 + x^2)`
- FriCAS `cA-3-161`: `(1 + x)/((2 + x)^(3/2)*(1 + 2*x + 2*x^2))`
- FriCAS `cA-3-194`: `Sqrt[2 + x + x^2]/(1 + x^2)`
- FriCAS `cA-3-294`: `(2 + x)/(Sqrt[2 + x + x^2]*(1 - 2*x + 2*x^2))`
- FriCAS `cA-4-51`: `Sqrt[2 + x]/(1 + x^2)`
- FriCAS `cA-4-279`: `1/((3 + x^2)*(1 + x + x^2)^(3/2))`
- FriCAS `cM-1-44`: `(2*E^x^2*Sqrt[1 + E^x^2 + E^(2*x^2)]*x)/(-2 + E^x^2)`
- FriCAS `cM-1-98`: `(E^x*(-1 + Log[1 + E^x])*(2 + Log[1 + E^x])^(3/2))/((1 + E^x)*(1 + Log[1 + E^x] + Log[1 + E^x]^2))`
- FriCAS `cM-1-126`: `(2*E^(2*x)*(2 + E^(4*x))^(3/2))/(-1 + E^(2*x))`
- FriCAS `cM-1-143`: `(1 + Log[x] + Log[x]^2)^(3/2)/(x*(3 + Log[x]^2))`
- FriCAS `cM-1-155`: `(2*E^(2*x)*(3 + E^(4*x))^(3/2))/(-1 + E^(2*x))`
- FriCAS `cM-1-259`: `-((E^x*(2 + E^x))/((-1 + E^x - E^(2*x))*(1 + E^x + E^(2*x))^(3/2)))`
- FriCAS `cM-1-286`: `(2*(1 + Log[1 + 2*x] + Log[1 + 2*x]^2)^(3/2))/((1 + 2*x)*(3 + Log[1 + 2*x]^2))`
- FriCAS `cM-1-327`: `(2*x*Sqrt[3 + Log[1 + x^2]])/((1 + x^2)*(2 + Log[1 + x^2]^2))`
- FriCAS `cM-1-348`: `(2*E^(2*x)*Sqrt[2 + E^(2*x) + E^(4*x)])/(1 + E^(2*x))`
- FriCAS `cM-2-67`: `-((E^x*(-2 + E^x))/(Sqrt[2 + E^x]*(-1 + E^x - 2*E^(2*x))))`
- FriCAS `cM-2-122`: `1/((1 + x)*(2 + Log[1 + x]^2)*(1 + Log[1 + x] + Log[1 + x]^2)^(3/2))`
- FriCAS `cM-2-157`: `(2*x*(3 + Log[1 + x^2])^(3/2))/((1 + x^2)*(2 + Log[1 + x^2]^2))`
- FriCAS `cM-2-249`: `(E^x*(-1 + Log[1 + E^x])*Sqrt[1 + Log[1 + E^x]])/((1 + E^x)*(1 + 2*Log[1 + E^x]^2))`
- FriCAS `cM-2-286`: `(E^x*(2 + E^(-x))^(3/2))/(-1 + E^x - E^(2*x))`
- FriCAS `cM-2-329`: `(2*E^x^2*Sqrt[1 + E^x^2]*x)/(2 + E^(2*x^2))`
- FriCAS `cM-2-351`: `(2*E^x^2*(2 + E^x^2)*x)/((1 + E^(2*x^2))^(3/2)*(1 + 2*E^x^2 + 2*E^(2*x^2)))`
- FriCAS `cM-2-352`: `(E^x*Sqrt[2 + Log[1 + E^x]])/((1 + E^x)*(1 - Log[1 + E^x] + 2*Log[1 + E^x]^2))`
- FriCAS `cM-3-104`: `(2*E^x^2*x)/(Sqrt[3 + E^x^2]*(3 + E^(2*x^2)))`
- FriCAS `cM-3-158`: `1/(x*Sqrt[1 + Log[x]]*(1 + Log[x]^2))`
- FriCAS `cM-3-159`: `(E^x*(1 + E^x + E^(2*x))^(3/2))/(1 + E^(2*x))`
- FriCAS `cM-3-240`: `(E^x^2*(2 + E^x^2)*x)/(Sqrt[3 + E^x^2]*(1 + E^(2*x^2)))`
- FriCAS `cM-4-42`: `E^x/((1 + E^x)*Sqrt[1 + Log[1 + E^x]]*(1 + Log[1 + E^x]^2))`
- FriCAS `cM-4-275`: `(2*E^(2*x)*(1 + E^(2*x) + E^(4*x))^(3/2))/(-2 + E^(2*x))`
- FriCAS `cM-4-340`: `(2*E^x^2*Sqrt[1 + E^x^2 + E^(2*x^2)]*x)/(1 + E^(2*x^2))`
