# Integration shootout: 259-transcendental-decision (2026-11-22)

Re-tabulated from the run of 2026-11-22: 79 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

561 cases, timeout 30 s per run, wall 2.1 min. Every answer was verified by its own system (derivative at three points, 10^-20). `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary; `wrong`: an answer failing its own check AND PossibleZeroQ[D[F, x] - f] in Mathilda -- the test Mathilda's own answers get (an answer failing its own check is translated into Mathilda and judged there; one PossibleZeroQ proves zero counts as verified, one it cannot decide or that does not translate as unverified).

## Closed (elementary + special) by class

| run | T transcendental (561) | all |
|---|---|---|
| Mathilda:Automatic | 70% (390) | 70% (390) |
| Mathilda:RischTranscendental | 62% (347) | 62% (347) |
| FriCAS | 64% (357) | 64% (357) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 287 | 103 | 170 | 0 | 1 | 0 | 0 | 0 | 0.04 |
| Mathilda:RischTranscendental | 288 | 59 | 210 | 2 | 2 | 0 | 0 | 0 | 0.04 |
| FriCAS | 270 | 87 | 0 | 189 | 0 | 15 | 0 | 0 | 0.11 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 37 | 37 | 37 |
| T | tanh-rand | 31 | 11 | 11 |
| T | tan-elem | 34 | 34 | 33 |
| T | tan-rand | 33 | 7 | 7 |
| T | tani-elem | 38 | 38 | 37 |
| T | tani-rand | 31 | 15 | 15 |
| T | gauss-elem | 32 | 32 | 21 |
| T | gauss-rand | 34 | 24 | 19 |
| T | lim-elem | 29 | 29 | 27 |
| T | lim-rand | 40 | 4 | 4 |
| T | prim-elem | 39 | 39 | 38 |
| T | prim-rand | 32 | 17 | 15 |
| T | dep-elem | 38 | 38 | 38 |
| T | dep-rand | 37 | 25 | 16 |
| T | f12-elem | 36 | 36 | 35 |
| T | f12-rand | 40 | 5 | 4 |

## WRONG (0)


## FALSE_CERT (0)


## DISAGREE (0)


## CORPUS? (0)


## AUTO-GAP (1)

- `tan-elem-32` AUTO-GAP Mathilda:RischTranscendental: `(x^3*Cot[x]^2*((2*Sec[x]^2*Tan[x])/x^2 - (2*Tan[x]^2)/x^3))/(1 + Tan[x]^2/x^2)^2 - (3*x^2*Cot[x]^2)/(1 + Tan[x]^2/x^2) + (2*x^3*Cot[x]*Csc[x]^2)/(1 + Tan[x]^2/x^2)`

## Only one system closes it

- Mathilda: 20
- FriCAS: 1

## Wrong answers from other systems (their bugs)

