# Integration shootout: 259-transcendental-decision (2026-11-22)

640 cases, timeout 30 s per run, wall 2.1 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | T transcendental (640) | all |
|---|---|---|
| Mathilda:Automatic | 73% (467) | 73% (467) |
| Mathilda:RischTranscendental | 65% (419) | 65% (419) |
| FriCAS | 67% (429) | 67% (429) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 354 | 113 | 172 | 0 | 1 | 0 | 0 | 0 | 0.03 |
| Mathilda:RischTranscendental | 355 | 64 | 217 | 2 | 2 | 0 | 0 | 0 | 0.03 |
| FriCAS | 333 | 96 | 0 | 192 | 0 | 19 | 0 | 0 | 0.11 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 40 | 40 | 39 |
| T | tanh-rand | 40 | 20 | 19 |
| T | tan-elem | 40 | 40 | 38 |
| T | tan-rand | 40 | 14 | 14 |
| T | tani-elem | 40 | 40 | 39 |
| T | tani-rand | 40 | 24 | 24 |
| T | gauss-elem | 40 | 40 | 28 |
| T | gauss-rand | 40 | 30 | 25 |
| T | lim-elem | 40 | 40 | 38 |
| T | lim-rand | 40 | 4 | 4 |
| T | prim-elem | 40 | 40 | 39 |
| T | prim-rand | 40 | 23 | 21 |
| T | dep-elem | 40 | 40 | 40 |
| T | dep-rand | 40 | 28 | 18 |
| T | f12-elem | 40 | 40 | 39 |
| T | f12-rand | 40 | 5 | 4 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (1)

- `tan-elem-32` AUTO-GAP Mathilda:RischTranscendental: `(x^3*Cot[x]^2*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*x^2*Cot[x]^2)/(1 + Tan[x]^2/x^2) + (2*x^3*Cot[x]*Csc[x]^2)/(1 + Tan[x]^2/x^2)`

## Only one system closes it

- Mathilda: 24
- FriCAS: 1

## Wrong answers from other systems (their bugs)

