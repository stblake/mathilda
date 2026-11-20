# Integration shootout: 259-transcendental-decision (2026-11-20)

640 cases, timeout 30 s per run, wall 5.6 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | T transcendental (640) | all |
|---|---|---|
| Mathilda:Automatic | 65% (415) | 65% (415) |
| Mathilda:RischTranscendental | 63% (406) | 63% (406) |
| FriCAS | 67% (428) | 67% (428) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 304 | 111 | 173 | 13 | 38 | 0 | 0 | 1 | 0.36 |
| Mathilda:RischTranscendental | 342 | 64 | 211 | 19 | 4 | 0 | 0 | 0 | 0.02 |
| FriCAS | 333 | 95 | 0 | 189 | 0 | 23 | 0 | 0 | 0.10 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 40 | 40 | 39 |
| T | tanh-rand | 40 | 16 | 20 |
| T | tan-elem | 40 | 39 | 40 |
| T | tan-rand | 40 | 12 | 14 |
| T | tani-elem | 40 | 40 | 39 |
| T | tani-rand | 40 | 22 | 23 |
| T | gauss-elem | 40 | 36 | 27 |
| T | gauss-rand | 40 | 27 | 24 |
| T | lim-elem | 40 | 40 | 40 |
| T | lim-rand | 40 | 2 | 4 |
| T | prim-elem | 40 | 40 | 36 |
| T | prim-rand | 40 | 18 | 21 |
| T | dep-elem | 40 | 40 | 39 |
| T | dep-rand | 40 | 28 | 19 |
| T | f12-elem | 40 | 39 | 38 |
| T | f12-rand | 40 | 3 | 5 |

## WRONG (1)

- `tani-rand-9` WRONG Mathilda:Automatic: `I*x^2*Log[-I + Tan[2*x]]`

## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (27)

- `tanh-elem-6` AUTO-GAP Mathilda:RischTranscendental: `-(E^Tanh[x]*Log[1 - Tanh[x]^2/x^2]) + (Log[x]*Sech[x]^2)/(2*x^2) - E^Tanh[x]*x*Log[1 - Tanh[x]^2/x^2]*Sech[x]^2 + Tanh[x]/(2*x^3) - (Log[x]*Tanh[x])/x^3 - (E^Tanh[x]*x*((-2*Sech[x]^2*Tanh[x])/x^2 + (2*Tanh[x]^2)/x^3))/(1`
- `tanh-elem-9` AUTO-GAP Mathilda:RischTranscendental: `(2*E^(x*Tanh[1 + x])*Sech[1 + x]^2*Tanh[1 + x])/(x*(1 - Tanh[1 + x]^2)^2) - E^(x*Tanh[1 + x])/(x^2*(1 - Tanh[1 + x]^2)) + (E^(x*Tanh[1 + x])*(x*Sech[1 + x]^2 + Tanh[1 + x]))/(x*(1 - Tanh[1 + x]^2))`
- `tanh-elem-12` AUTO-GAP Mathilda:RischTranscendental: `(-8*E^Tanh[2*x]*Sech[2*x]^2*Tanh[2*x])/((1 - Tanh[2*x]^2)*(2 + Tanh[2*x]^2)^2) + (8*E^Tanh[2*x]*Sech[2*x]^2*Tanh[2*x])/((1 - Tanh[2*x]^2)^2*(2 + Tanh[2*x]^2)) + (4*E^Tanh[2*x]*Sech[2*x]^2)/((1 - Tanh[2*x]^2)*(2 + Tanh[2*`
- `tanh-elem-13` AUTO-GAP Mathilda:RischTranscendental: `(3*Log[x + Tanh[x]]*Sech[x]^2)/(2*x) - (3*Log[x + Tanh[x]]*Tanh[x])/(2*x^2) + (3*(1 + Sech[x]^2)*Tanh[x])/(2*x*(x + Tanh[x])) - (2*x*Sech[x]^2)/((1 + Tanh[x])*(2 + Tanh[x]^2)^2) - (x*Csch[x]*Sech[x])/((1 + Tanh[x])^2*(2 `
- `tanh-elem-14` AUTO-GAP Mathilda:RischTranscendental: `-((E^Tanh[2*x]*Coth[2*x])/(1 + x)^2) - (2*E^Tanh[2*x]*Csch[2*x]^2)/(1 + x) + (2*E^Tanh[2*x]*Csch[2*x]*Sech[2*x])/(1 + x)`
- `tanh-elem-18` AUTO-GAP Mathilda:RischTranscendental: `-1 - (3*Log[x]*Sech[x]^2*Tanh[x])/((1 - Tanh[x]^2)*(2 + Tanh[x]^2)^2) + (3*Log[x]*Sech[x]^2*Tanh[x])/((1 - Tanh[x]^2)^2*(2 + Tanh[x]^2)) + 3/(2*x*(1 - Tanh[x]^2)*(2 + Tanh[x]^2))`
- `tanh-elem-22` AUTO-GAP Mathilda:RischTranscendental: `-((x*Log[x]*Log[x + Tanh[x]/x]*((2*Sech[x]^2*Tanh[x])/x^2 - (2*Tanh[x]^2)/x^3))/(2 + Tanh[x]^2/x^2)^2) + Log[x + Tanh[x]/x]/(2 + Tanh[x]^2/x^2) + (Log[x]*Log[x + Tanh[x]/x])/(2 + Tanh[x]^2/x^2) + (x*Log[x]*(1 + Sech[x]^2`
- `tanh-elem-25` AUTO-GAP Mathilda:RischTranscendental: `(-4*x*Log[x + Tanh[2*x]]*Sech[2*x]^2*Tanh[2*x]^3)/(2 + Tanh[2*x]^2)^2 + (4*x*Log[x + Tanh[2*x]]*Sech[2*x]^2*Tanh[2*x])/(2 + Tanh[2*x]^2) + (Log[x + Tanh[2*x]]*Tanh[2*x]^2)/(2 + Tanh[2*x]^2) + (x*(1 + 2*Sech[2*x]^2)*Tanh[`
- `tanh-elem-30` AUTO-GAP Mathilda:RischTranscendental: `(-2*x*Csch[x]*Log[x]*Sech[x])/(1 + Tanh[x])^2 - (x*Csch[x]*Log[x + Tanh[x]]*Sech[x])/(1 + Tanh[x])^2 + (2*Coth[x])/(1 + Tanh[x]) + (2*Coth[x]*Log[x])/(1 + Tanh[x]) - (2*x*Csch[x]^2*Log[x])/(1 + Tanh[x]) + (Coth[x]*Log[x `
- `tanh-elem-39` AUTO-GAP Mathilda:RischTranscendental: `(-2*E^(x*Tanh[x^2])*x*Sech[x^2]^2)/(2 + Tanh[x^2]^2)^2 - (E^(x*Tanh[x^2])*x*Csch[x^2]^2)/(2 + Tanh[x^2]^2) + (E^(x*Tanh[x^2])*Coth[x^2]*(2*x^2*Sech[x^2]^2 + Tanh[x^2]))/(2*(2 + Tanh[x^2]^2))`
- `tan-elem-3` AUTO-GAP Mathilda:RischTranscendental: `(-6*E^x*Log[x + Tan[2*x]]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 + (3*E^x*Log[x + Tan[2*x]])/(2*(1 + Tan[2*x]^2)) + (3*E^x*(1 + 2*Sec[2*x]^2))/(2*(x + Tan[2*x])*(1 + Tan[2*x]^2))`
- `tan-elem-7` AUTO-GAP Mathilda:RischTranscendental: `2*E^Tan[1 + x]*Log[1 + Tan[1 + x]^2] + 2*E^Tan[1 + x]*x*Log[1 + Tan[1 + x]^2]*Sec[1 + x]^2 + (4*E^Tan[1 + x]*x*Sec[1 + x]^2*Tan[1 + x])/(1 + Tan[1 + x]^2)`
- `tan-elem-9` AUTO-GAP Mathilda:RischTranscendental: `E^(x*Tan[x])/x - (Sec[x]^2*Tan[x]^2)/(2*(1 + Tan[x])^2) + (Sec[x]^2*Tan[x])/(1 + Tan[x]) + E^(x*Tan[x])*Log[x]*(x*Sec[x]^2 + Tan[x])`
- `tan-elem-16` AUTO-GAP Mathilda:RischTranscendental: `(E^(Tan[x]^2/x^2)*Log[x + Tan[x]/x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*Tan[x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(2*(1 + Tan[x]^2/x^2)^2) + (3*Sec[x]^2)/(2*(1 + Tan[x]^2/x^`
- `tan-elem-27` AUTO-GAP Mathilda:RischTranscendental: `(2*x*Log[Tan[2*x]]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 - (2*x*Log[Tan[2*x]]*Log[1 + Tan[2*x]^2]*Sec[2*x]^2*Tan[2*x])/(1 + Tan[2*x]^2)^2 + (Log[Tan[2*x]]*Log[1 + Tan[2*x]^2])/(2*(1 + Tan[2*x]^2)) + (x*Csc[2*x]*Log[1 +`
- `tan-elem-28` AUTO-GAP Mathilda:RischTranscendental: `2/(1 + x) - (2*x*Log[x])/(1 + x)^2 + (2*Log[x])/(1 + x) - (3*x*Log[x + Tan[x]/x]*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(2*(1 + Tan[x]/x)*(1 + Tan[x]^2/x^2)^2) - (3*x*Log[x + Tan[x]/x]*(Sec[x]^2/x - Tan[x]/x^2))/(`
- `tan-elem-32` AUTO-GAP Mathilda:RischTranscendental: `(x^3*Cot[x]^2*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*x^2*Cot[x]^2)/(1 + Tan[x]^2/x^2) + (2*x^3*Cot[x]*Csc[x]^2)/(1 + Tan[x]^2/x^2)`
- `tani-elem-10` AUTO-GAP Mathilda:RischTranscendental: `Log[x]/(x*(I + x)) + Log[I + x]/x^2 - (Log[x]*Log[I + x])/x^2 - ((2*I)*Sec[x^2]^2)/(-I + Tan[x^2])^2 + ((2*I)*Log[-I + Tan[x^2]]*Sec[x^2]^2)/(-I + Tan[x^2])^2 + (I*Log[-I + Tan[x^2]])/(x^2*(-I + Tan[x^2]))`
- `tani-elem-16` AUTO-GAP Mathilda:RischTranscendental: `-(E^(I*x)*Log[-I + Tan[x^2]]) - I*E^(I*x)*x*Log[-I + Tan[x^2]] - (2*E^(I*x)*x^2*Sec[x^2]^2)/(-I + Tan[x^2])`
- `tani-elem-35` AUTO-GAP Mathilda:RischTranscendental: `(1 + I)*E^(2*Tan[1 + x]) + (2 + 2*I)*E^(2*Tan[1 + x])*x*Sec[1 + x]^2 - (E^Tan[1 + x]*Log[I + x]*Sec[1 + x]^2)/(2*(1 + Tan[1 + x])^2) + E^Tan[1 + x]/(2*(I + x)*(1 + Tan[1 + x])) + (E^Tan[1 + x]*Log[I + x]*Sec[1 + x]^2)/(2`
- `gauss-elem-2` AUTO-GAP Mathilda:RischTranscendental: `1/((1 + E^((-I)*(1 + x)))*(1 + E^(I*(2 + x)))) - (I*E^(I*(2 + x))*x)/((1 + E^((-I)*(1 + x)))*(1 + E^(I*(2 + x)))^2) + (I*E^((-I)*(1 + x))*x)/((1 + E^((-I)*(1 + x)))^2*(1 + E^(I*(2 + x)))) + ((2*I)*x*Sec[2 + x]^2)/(1 + E^`
- `gauss-elem-13` AUTO-GAP Mathilda:RischTranscendental: `(E^(-I)*Sech[1 + x]^2)/(-I + Tan[1 + x]) - (E^(-I)*Sec[1 + x]^2*Tanh[1 + x])/(-I + Tan[1 + x])^2`
- `gauss-elem-25` AUTO-GAP Mathilda:RischTranscendental: `((-2*I)*E^(I*(1 + x))*Log[I + Tanh[1 + x]])/x^2 - (2*E^(I*(1 + x))*Log[I + Tanh[1 + x]])/x - (E^I*x*Sech[1 + x]^2)/((1 + E^(I*(2 + x)))*(-I + Tanh[1 + x])^2) + E^I/((1 + E^(I*(2 + x)))*(-I + Tanh[1 + x])) - (I*E^(I + I*(`
- `lim-elem-9` AUTO-GAP Mathilda:RischTranscendental: `(-2*(1 + E^x)*Log[E^x + x]*Log[1 + Tanh[x]^2])/(E^x + x) - (2*Log[E^x + x]^2*Sech[x]^2*Tanh[x])/(1 + Tanh[x]^2)`
- `lim-elem-23` AUTO-GAP Mathilda:RischTranscendental: `(2*(1 + E^x)*Log[E^x + x])/((E^x + x)*(1 + Coth[x])) - (x*Csch[x]^2*Log[E^x + x])/((1 + Coth[x])*(1 + Coth[x]^2)) + (Csch[x]^2*Log[E^x + x]^2)/(1 + Coth[x])^2 + ((1 + E^x)*x*Log[1 + Coth[x]])/((E^x + x)*(1 + Coth[x]^2)) `
- `lim-elem-28` AUTO-GAP Mathilda:RischTranscendental: `-((Log[1 + x^2]^2*Log[1 + Tan[x]]*Sec[x]^2*Tan[x])/(1 + Tan[x]^2)^2) - (x*Log[1 + x^2]^2*Log[1 + Tan[x]]*Sec[x]^2*Tan[x])/(1 + Tan[x]^2)^2 + (2*x*Log[1 + x^2]*Log[1 + Tan[x]])/((1 + x^2)*(1 + Tan[x]^2)) + (2*x^2*Log[1 + `
- `lim-elem-40` AUTO-GAP Mathilda:RischTranscendental: `Log[1 + Tanh[2*x]^2]/x^2 + Sech[2*x]^2 - (4*Sech[2*x]^2*Tanh[2*x])/(x*(1 + Tanh[2*x]^2))`

## Only one system closes it

- Mathilda: 18
- FriCAS: 5

## Wrong answers from other systems (their bugs)

