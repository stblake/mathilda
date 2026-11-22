# Integration shootout: 260-radical-bottom (2026-11-22)

420 cases, timeout 30 s per run, wall 8.6 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | B radical at the bottom (420) | all |
|---|---|---|
| Mathilda:Automatic | 51% (215) | 51% (215) |
| Mathilda:ParallelMixedTower | 36% (150) | 36% (150) |
| Mathilda:Bronstein | 45% (189) | 45% (189) |
| FriCAS | 52% (217) | 52% (217) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 187 | 28 | 197 | 7 | 1 | 0 | 0 | 0 | 0.12 |
| Mathilda:ParallelMixedTower | 150 | 0 | 48 | 188 | 34 | 0 | 0 | 0 | 0.25 |
| Mathilda:Bronstein | 189 | 0 | 226 | 5 | 0 | 0 | 0 | 0 | 0.04 |
| FriCAS | 180 | 37 | 0 | 180 | 7 | 14 | 0 | 2 | 0.13 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| B | one-elem | 60 | 60 | 58 |
| B | one-rand | 60 | 18 | 12 |
| B | mid-elem | 60 | 60 | 58 |
| B | mid-rand | 60 | 4 | 2 |
| B | two-elem | 60 | 59 | 54 |
| B | two-rand | 60 | 2 | 2 |
| B | pert-nonelem | 60 | 14 | 31 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (2)

- `two-elem-6` AUTO-GAP Mathilda:Bronstein: `(-3*E^(x + x^2)*x)/(2 + x^2)^(3/2) + (3*E^(x + x^2)*(1 + 2*x))/Sqrt[2 + x^2] - (Sqrt[2 + x^2]*(1 + x^3)^(1/3))/(1 + x) + (1 + x^2/(1 + x^3)^(2/3))/(2*(x + (1 + x^3)^(1/3))) - (x^2*Sqrt[2 + x^2]*Log[1 + x])/(1 + x^3)^(2/3`
- `two-elem-49` AUTO-GAP Mathilda:Bronstein: `(E^(x + x^2)*x^3)/(1 + x^3)^(1/3) + (E^(x + x^2)*(1 + x^3)^(2/3))/2 + (E^(x + x^2)*x*(1 + 2*x)*(1 + x^3)^(2/3))/2 + (E^x^2*x^3*(1 + x^3)^(1/3))/Sqrt[3 + x^4] + E^(2*x)*x*Sqrt[3 + x^4] + (E^x^2*x^2*Sqrt[3 + x^4])/(2*(1 + `

## Only one system closes it

- Mathilda: 10
- FriCAS: 1

## Wrong answers from other systems (their bugs)

- FriCAS `one-rand-29`: `E^x^2/Sqrt[x]`
- FriCAS `one-rand-40`: `E^x^2*Sqrt[x]`
