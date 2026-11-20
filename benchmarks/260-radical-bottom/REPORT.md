# Integration shootout: 260-radical-bottom (2026-11-21)

420 cases, timeout 30 s per run, wall 8.5 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | B radical at the bottom (420) | all |
|---|---|---|
| Mathilda:Automatic | 49% (204) | 49% (204) |
| Mathilda:ParallelMixedTower | 35% (148) | 35% (148) |
| Mathilda:Bronstein | 37% (155) | 37% (155) |
| FriCAS | 52% (217) | 52% (217) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 173 | 31 | 171 | 45 | 0 | 0 | 0 | 0 | 0.07 |
| Mathilda:ParallelMixedTower | 148 | 0 | 48 | 185 | 39 | 0 | 0 | 0 | 0.20 |
| Mathilda:Bronstein | 155 | 0 | 200 | 65 | 0 | 0 | 0 | 0 | 0.04 |
| FriCAS | 180 | 37 | 0 | 180 | 7 | 14 | 0 | 2 | 0.13 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| B | one-elem | 60 | 60 | 58 |
| B | one-rand | 60 | 18 | 12 |
| B | mid-elem | 60 | 55 | 58 |
| B | mid-rand | 60 | 4 | 2 |
| B | two-elem | 60 | 52 | 54 |
| B | two-rand | 60 | 2 | 2 |
| B | pert-nonelem | 60 | 15 | 31 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (2)

- `mid-elem-23` AUTO-GAP Mathilda:ParallelMixedTower: `(E^(x + Sqrt[1 + x^3])*(1 + (3*x^2)/(2*Sqrt[1 + x^3])))/2 + (x*(1 + x^3)*(1 + (3*x^2)/(2*Sqrt[1 + x^3]))*Log[1 + Sqrt[1 + x^3]])/(x + Sqrt[1 + x^3]) + (3*x^3*Sqrt[1 + x^3]*Log[x + Sqrt[1 + x^3]])/(2*(1 + Sqrt[1 + x^3])) `
- `mid-elem-29` AUTO-GAP Mathilda:ParallelMixedTower: `-E^(x*(1 + x^3)^(1/3)) + (3*E^(1 + x^3)^(1/3))/(1 + x) - E^(x*(1 + x^3)^(1/3))*Log[x] - E^(x*(1 + x^3)^(1/3))*x*(x^3/(1 + x^3)^(2/3) + (1 + x^3)^(1/3))*Log[x] + (3*E^(1 + x^3)^(1/3)*x^2*Log[1 + x])/(1 + x^3)^(2/3)`

## Only one system closes it

- FriCAS: 11
- Mathilda: 8

## Wrong answers from other systems (their bugs)

- FriCAS `one-rand-29`: `E^x^2/Sqrt[x]`
- FriCAS `one-rand-40`: `E^x^2*Sqrt[x]`
