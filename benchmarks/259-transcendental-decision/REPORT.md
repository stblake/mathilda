# Integration shootout: 259-transcendental-decision (2026-11-20)

640 cases, timeout 30 s per run, wall 9.0 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | T transcendental (640) | all |
|---|---|---|
| Mathilda:Automatic | 65% (414) | 65% (414) |
| Mathilda:ParallelMixedTower | 45% (287) | 45% (287) |
| Mathilda:RischTranscendental | 63% (405) | 63% (405) |
| FriCAS | 67% (429) | 67% (429) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 303 | 111 | 173 | 13 | 39 | 0 | 0 | 1 | 0.33 |
| Mathilda:ParallelMixedTower | 287 | 0 | 134 | 191 | 28 | 0 | 0 | 0 | 0.12 |
| Mathilda:RischTranscendental | 341 | 64 | 211 | 19 | 5 | 0 | 0 | 0 | 0.02 |
| FriCAS | 333 | 96 | 0 | 192 | 0 | 19 | 0 | 0 | 0.11 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 40 | 40 | 39 |
| T | tanh-rand | 40 | 16 | 19 |
| T | tan-elem | 40 | 39 | 38 |
| T | tan-rand | 40 | 12 | 14 |
| T | tani-elem | 40 | 40 | 39 |
| T | tani-rand | 40 | 22 | 24 |
| T | gauss-elem | 40 | 36 | 28 |
| T | gauss-rand | 40 | 27 | 25 |
| T | lim-elem | 40 | 40 | 38 |
| T | lim-rand | 40 | 2 | 4 |
| T | prim-elem | 40 | 40 | 39 |
| T | prim-rand | 40 | 18 | 21 |
| T | dep-elem | 40 | 40 | 40 |
| T | dep-rand | 40 | 28 | 18 |
| T | f12-elem | 40 | 39 | 39 |
| T | f12-rand | 40 | 3 | 4 |

## WRONG (1)

- `tani-rand-9` WRONG Mathilda:Automatic: `I*x^2*Log[-I + Tan[2*x]]`

## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (28)

- `tanh-elem-1` AUTO-GAP Mathilda:ParallelMixedTower: `-((E^Tanh[x]*x^2*Coth[x]*(Sech[x]^2/x - Tanh[x]/x^2))/(1 + Tanh[x]/x)^2) + (2*E^Tanh[x]*x*Coth[x])/(1 + Tanh[x]/x) - (E^Tanh[x]*x^2*Csch[x]^2)/(1 + Tanh[x]/x) + (E^Tanh[x]*x^2*Csch[x]*Sech[x])/(1 + Tanh[x]/x)`
- `tanh-elem-6` AUTO-GAP Mathilda:ParallelMixedTower: `-(E^Tanh[x]*Log[1 - Tanh[x]^2/x^2]) + (Log[x]*Sech[x]^2)/(2*x^2) - E^Tanh[x]*x*Log[1 - Tanh[x]^2/x^2]*Sech[x]^2 + Tanh[x]/(2*x^3) - (Log[x]*Tanh[x])/x^3 - (E^Tanh[x]*x*((-2*Sech[x]^2*Tanh[x])/x^2 + (2*Tanh[x]^2)/x^3))/(1`
- `tanh-elem-9` AUTO-GAP Mathilda:ParallelMixedTower: `(2*E^(x*Tanh[1 + x])*Sech[1 + x]^2*Tanh[1 + x])/(x*(1 - Tanh[1 + x]^2)^2) - E^(x*Tanh[1 + x])/(x^2*(1 - Tanh[1 + x]^2)) + (E^(x*Tanh[1 + x])*(x*Sech[1 + x]^2 + Tanh[1 + x]))/(x*(1 - Tanh[1 + x]^2))`
- `tanh-elem-12` AUTO-GAP Mathilda:ParallelMixedTower: `(-8*E^Tanh[2*x]*Sech[2*x]^2*Tanh[2*x])/((1 - Tanh[2*x]^2)*(2 + Tanh[2*x]^2)^2) + (8*E^Tanh[2*x]*Sech[2*x]^2*Tanh[2*x])/((1 - Tanh[2*x]^2)^2*(2 + Tanh[2*x]^2)) + (4*E^Tanh[2*x]*Sech[2*x]^2)/((1 - Tanh[2*x]^2)*(2 + Tanh[2*`
- `tanh-elem-13` AUTO-GAP Mathilda:ParallelMixedTower: `(3*Log[x + Tanh[x]]*Sech[x]^2)/(2*x) - (3*Log[x + Tanh[x]]*Tanh[x])/(2*x^2) + (3*(1 + Sech[x]^2)*Tanh[x])/(2*x*(x + Tanh[x])) - (2*x*Sech[x]^2)/((1 + Tanh[x])*(2 + Tanh[x]^2)^2) - (x*Csch[x]*Sech[x])/((1 + Tanh[x])^2*(2 `
- `tanh-elem-14` AUTO-GAP Mathilda:ParallelMixedTower: `-((E^Tanh[2*x]*Coth[2*x])/(1 + x)^2) - (2*E^Tanh[2*x]*Csch[2*x]^2)/(1 + x) + (2*E^Tanh[2*x]*Csch[2*x]*Sech[2*x])/(1 + x)`
- `tanh-elem-18` AUTO-GAP Mathilda:ParallelMixedTower: `-1 - (3*Log[x]*Sech[x]^2*Tanh[x])/((1 - Tanh[x]^2)*(2 + Tanh[x]^2)^2) + (3*Log[x]*Sech[x]^2*Tanh[x])/((1 - Tanh[x]^2)^2*(2 + Tanh[x]^2)) + 3/(2*x*(1 - Tanh[x]^2)*(2 + Tanh[x]^2))`
- `tanh-elem-22` AUTO-GAP Mathilda:ParallelMixedTower: `-((x*Log[x]*Log[x + Tanh[x]/x]*((2*Sech[x]^2*Tanh[x])/x^2 - (2*Tanh[x]^2)/x^3))/(2 + Tanh[x]^2/x^2)^2) + Log[x + Tanh[x]/x]/(2 + Tanh[x]^2/x^2) + (Log[x]*Log[x + Tanh[x]/x])/(2 + Tanh[x]^2/x^2) + (x*Log[x]*(1 + Sech[x]^2`
- `tanh-elem-25` AUTO-GAP Mathilda:ParallelMixedTower: `(-4*x*Log[x + Tanh[2*x]]*Sech[2*x]^2*Tanh[2*x]^3)/(2 + Tanh[2*x]^2)^2 + (4*x*Log[x + Tanh[2*x]]*Sech[2*x]^2*Tanh[2*x])/(2 + Tanh[2*x]^2) + (Log[x + Tanh[2*x]]*Tanh[2*x]^2)/(2 + Tanh[2*x]^2) + (x*(1 + 2*Sech[2*x]^2)*Tanh[`
- `tanh-elem-30` AUTO-GAP Mathilda:ParallelMixedTower: `(-2*x*Csch[x]*Log[x]*Sech[x])/(1 + Tanh[x])^2 - (x*Csch[x]*Log[x + Tanh[x]]*Sech[x])/(1 + Tanh[x])^2 + (2*Coth[x])/(1 + Tanh[x]) + (2*Coth[x]*Log[x])/(1 + Tanh[x]) - (2*x*Csch[x]^2*Log[x])/(1 + Tanh[x]) + (Coth[x]*Log[x `
- `tanh-elem-39` AUTO-GAP Mathilda:ParallelMixedTower: `(-2*E^(x*Tanh[x^2])*x*Sech[x^2]^2)/(2 + Tanh[x^2]^2)^2 - (E^(x*Tanh[x^2])*x*Csch[x^2]^2)/(2 + Tanh[x^2]^2) + (E^(x*Tanh[x^2])*Coth[x^2]*(2*x^2*Sech[x^2]^2 + Tanh[x^2]))/(2*(2 + Tanh[x^2]^2))`
- `tan-elem-3` AUTO-GAP Mathilda:ParallelMixedTower: `(-6*E^x*Log[x + Tan[2*x]]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 + (3*E^x*Log[x + Tan[2*x]])/(2*(1 + Tan[2*x]^2)) + (3*E^x*(1 + 2*Sec[2*x]^2))/(2*(x + Tan[2*x])*(1 + Tan[2*x]^2))`
- `tan-elem-7` AUTO-GAP Mathilda:ParallelMixedTower: `2*E^Tan[1 + x]*Log[1 + Tan[1 + x]^2] + 2*E^Tan[1 + x]*x*Log[1 + Tan[1 + x]^2]*Sec[1 + x]^2 + (4*E^Tan[1 + x]*x*Sec[1 + x]^2*Tan[1 + x])/(1 + Tan[1 + x]^2)`
- `tan-elem-9` AUTO-GAP Mathilda:ParallelMixedTower: `E^(x*Tan[x])/x - (Sec[x]^2*Tan[x]^2)/(2*(1 + Tan[x])^2) + (Sec[x]^2*Tan[x])/(1 + Tan[x]) + E^(x*Tan[x])*Log[x]*(x*Sec[x]^2 + Tan[x])`
- `tan-elem-16` AUTO-GAP Mathilda:ParallelMixedTower: `(E^(Tan[x]^2/x^2)*Log[x + Tan[x]/x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*Tan[x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(2*(1 + Tan[x]^2/x^2)^2) + (3*Sec[x]^2)/(2*(1 + Tan[x]^2/x^`
- `tan-elem-27` AUTO-GAP Mathilda:ParallelMixedTower: `(2*x*Log[Tan[2*x]]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 - (2*x*Log[Tan[2*x]]*Log[1 + Tan[2*x]^2]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 + (Log[Tan[2*x]]*Log[1 + Tan[2*x]^2])/(2*(1 + Tan[2*x]^2)) + (x*Csc[2*x]*Log[1 +`
- `tan-elem-28` AUTO-GAP Mathilda:ParallelMixedTower: `2/(1 + x) - (2*x*Log[x])/(1 + x)^2 + (2*Log[x])/(1 + x) - (3*x*Log[x + Tan[x]/x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(2*(1 + Tan[x]/x)*(1 + Tan[x]^2/x^2)^2) - (3*x*Log[x + Tan[x]/x]*(Sec[x]^2/x - Tan[x]/x^2))/(`
- `tani-elem-10` AUTO-GAP Mathilda:RischTranscendental: `Log[x]/(x*(I + x)) + Log[I + x]/x^2 - (Log[x]*Log[I + x])/x^2 - ((2*I)*Sec[x^2]^2)/(-I + Tan[x^2])^2 + ((2*I)*Log[-I + Tan[x^2]]*Sec[x^2]^2)/(-I + Tan[x^2])^2 + (I*Log[-I + Tan[x^2]])/(x^2*(-I + Tan[x^2]))`
- `tani-elem-16` AUTO-GAP Mathilda:ParallelMixedTower: `-(E^(I*x)*Log[-I + Tan[x^2]]) - I*E^(I*x)*x*Log[-I + Tan[x^2]] - (2*E^(I*x)*x^2*Sec[x^2]^2)/(-I + Tan[x^2])`
- `tani-elem-21` AUTO-GAP Mathilda:ParallelMixedTower: `(Log[I + x]*Sec[1 + x]^2)/(2*(-I + Tan[1 + x])^2) - (Log[I + x]*Log[-I + Tan[1 + x]]*Sec[1 + x]^2)/(2*(-I + Tan[1 + x])^2) + Log[-I + Tan[1 + x]]/(2*(I + x)*(-I + Tan[1 + x]))`
- `tani-elem-35` AUTO-GAP Mathilda:ParallelMixedTower: `(1 + I)*E^(2*Tan[1 + x]) + (2 + 2*I)*E^(2*Tan[1 + x])*x*Sec[1 + x]^2 - (E^Tan[1 + x]*Log[I + x]*Sec[1 + x]^2)/(2*(1 + Tan[1 + x])^2) + E^Tan[1 + x]/(2*(I + x)*(1 + Tan[1 + x])) + (E^Tan[1 + x]*Log[I + x]*Sec[1 + x]^2)/(2`
- `gauss-elem-2` AUTO-GAP Mathilda:RischTranscendental: `1/((1 + E^((-I)*(1 + x)))*(1 + E^(I*(2 + x)))) - (I*E^(I*(2 + x))*x)/((1 + E^((-I)*(1 + x)))*(1 + E^(I*(2 + x)))^2) + (I*E^((-I)*(1 + x))*x)/((1 + E^((-I)*(1 + x)))^2*(1 + E^(I*(2 + x)))) + ((2*I)*x*Sec[2 + x]^2)/(1 + E^`
- `gauss-elem-13` AUTO-GAP Mathilda:RischTranscendental: `(E^(-I)*Sech[1 + x]^2)/(-I + Tan[1 + x]) - (E^(-I)*Sec[1 + x]^2*Tanh[1 + x])/(-I + Tan[1 + x])^2`
- `lim-elem-9` AUTO-GAP Mathilda:ParallelMixedTower: `(-2*(1 + E^x)*Log[E^x + x]*Log[1 + Tanh[x]^2])/(E^x + x) - (2*Log[E^x + x]^2*Sech[x]^2*Tanh[x])/(1 + Tanh[x]^2)`
- `lim-elem-12` AUTO-GAP Mathilda:ParallelMixedTower: `(4*(1 + E^x)*Log[1 + (1 + E^x)^(-2)]*Log[E^x + x])/((1 + (1 + E^x)^(-2))*(E^x + x)) - (4*E^x*Log[E^x + x]^2)/((1 + E^x)^3*(1 + (1 + E^x)^(-2))^2) + (4*E^x*Log[1 + (1 + E^x)^(-2)]*Log[E^x + x]^2)/((1 + E^x)^3*(1 + (1 + E^`
- `lim-elem-23` AUTO-GAP Mathilda:RischTranscendental: `(2*(1 + E^x)*Log[E^x + x])/((E^x + x)*(1 + Coth[x])) - (x*Csch[x]^2*Log[E^x + x])/((1 + Coth[x])*(1 + Coth[x]^2)) + (Csch[x]^2*Log[E^x + x]^2)/(1 + Coth[x])^2 + ((1 + E^x)*x*Log[1 + Coth[x]])/((E^x + x)*(1 + Coth[x]^2)) `
- `lim-elem-28` AUTO-GAP Mathilda:ParallelMixedTower: `-((Log[1 + x^2]^2*Log[1 + Tan[x]]*Sec[x]^2*Tan[x])/(1 + Tan[x]^2)^2) - (x*Log[1 + x^2]^2*Log[1 + Tan[x]]*Sec[x]^2*Tan[x])/(1 + Tan[x]^2)^2 + (2*x*Log[1 + x^2]*Log[1 + Tan[x]])/((1 + x^2)*(1 + Tan[x]^2)) + (2*x^2*Log[1 + `
- `lim-elem-40` AUTO-GAP Mathilda:ParallelMixedTower: `Log[1 + Tanh[2*x]^2]/x^2 + Sech[2*x]^2 - (4*Sech[2*x]^2*Tanh[2*x])/(x*(1 + Tanh[2*x]^2))`

## Only one system closes it

- Mathilda: 18
- FriCAS: 4

## Wrong answers from other systems (their bugs)

