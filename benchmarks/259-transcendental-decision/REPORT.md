# Integration shootout: 259-transcendental-decision (2026-11-21)

640 cases, timeout 30 s per run, wall 4.8 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check.

## Closed (elementary + special) by class

| run | T transcendental (640) | all |
|---|---|---|
| Mathilda:Automatic | 70% (447) | 70% (447) |
| Mathilda:ParallelMixedTower | 45% (289) | 45% (289) |
| Mathilda:RischTranscendental | 65% (413) | 65% (413) |
| FriCAS | 67% (429) | 67% (429) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 348 | 99 | 180 | 8 | 5 | 0 | 0 | 0 | 0.04 |
| Mathilda:ParallelMixedTower | 289 | 0 | 134 | 191 | 26 | 0 | 0 | 0 | 0.12 |
| Mathilda:RischTranscendental | 349 | 64 | 212 | 11 | 4 | 0 | 0 | 0 | 0.02 |
| FriCAS | 333 | 96 | 0 | 192 | 0 | 19 | 0 | 0 | 0.11 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 40 | 40 | 39 |
| T | tanh-rand | 40 | 17 | 19 |
| T | tan-elem | 40 | 39 | 38 |
| T | tan-rand | 40 | 12 | 14 |
| T | tani-elem | 40 | 40 | 39 |
| T | tani-rand | 40 | 23 | 24 |
| T | gauss-elem | 40 | 37 | 28 |
| T | gauss-rand | 40 | 30 | 25 |
| T | lim-elem | 40 | 40 | 38 |
| T | lim-rand | 40 | 2 | 4 |
| T | prim-elem | 40 | 40 | 39 |
| T | prim-rand | 40 | 18 | 21 |
| T | dep-elem | 40 | 40 | 40 |
| T | dep-rand | 40 | 28 | 18 |
| T | f12-elem | 40 | 39 | 39 |
| T | f12-rand | 40 | 3 | 4 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (1)

- `tan-elem-32` AUTO-GAP Mathilda:ParallelMixedTower: `(x^3*Cot[x]^2*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*x^2*Cot[x]^2)/(1 + Tan[x]^2/x^2) + (2*x^3*Cot[x]*Csc[x]^2)/(1 + Tan[x]^2/x^2)`

## Only one system closes it

- Mathilda: 21
- FriCAS: 2

## Wrong answers from other systems (their bugs)

