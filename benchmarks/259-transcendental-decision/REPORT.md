# Integration shootout: 259-transcendental-decision (2026-11-22)

Re-tabulated from the run of 2026-11-22: 79 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

561 cases, timeout 30 s per run, wall 2.1 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (561) | all |
|---|---|---|
| Mathilda:Automatic | 69% (388) | 69% (388) |
| Mathilda:RischTranscendental | 62% (346) | 62% (346) |
| FriCAS | 63% (355) | 63% (355) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 285 | 103 | 170 | 0 | 1 | 0 | 2 | 0 | 0.04 |
| Mathilda:RischTranscendental | 287 | 59 | 210 | 2 | 2 | 0 | 1 | 0 | 0.04 |
| FriCAS | 270 | 85 | 0 | 193 | 0 | 1 | 6 | 6 | 0.11 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 390 | 388 | 0 | 2 | 0 | 0 | 388 | 0 |
| Mathilda:RischTranscendental | 347 | 346 | 0 | 1 | 0 | 0 | 346 | 0 |
| FriCAS | 367 | 355 | 6 | 6 | 0 | 0 | 355 | 6 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 37 | 37 | 35 |
| T | tanh-rand | 31 | 11 | 8 |
| T | tan-elem | 34 | 34 | 34 |
| T | tan-rand | 33 | 7 | 7 |
| T | tani-elem | 38 | 38 | 38 |
| T | tani-rand | 31 | 15 | 16 |
| T | gauss-elem | 32 | 31 | 21 |
| T | gauss-rand | 34 | 24 | 21 |
| T | lim-elem | 29 | 29 | 29 |
| T | lim-rand | 40 | 4 | 3 |
| T | prim-elem | 39 | 39 | 39 |
| T | prim-rand | 32 | 17 | 13 |
| T | dep-elem | 38 | 38 | 34 |
| T | dep-rand | 37 | 25 | 16 |
| T | f12-elem | 36 | 36 | 36 |
| T | f12-rand | 40 | 5 | 5 |

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

- FriCAS `prim-rand-22`: `x*Log[x]*LogIntegral[x]`
- FriCAS `prim-rand-38`: `Log[x]*LogIntegral[x]`
- FriCAS `dep-elem-9`: `-((x*Log[-1 + x])/(1 + x)) + Log[x]/(1 + x) + (2*x^3*Log[x])/(-1 + x^2) + Log[1 + x]/x - (x*Log[1 + x])/(-1 + x) - Log[-1 + x]*Log[1 + x] + x*Log[-1 + x^2] + 2*`
- FriCAS `dep-elem-12`: `(2*x^3*Log[-1 + x])/(-1 + x^2) + (2*x^2*Log[x])/(1 + x) + 2*x*Log[1 + x] + 4*x*Log[x]*Log[1 + x] + (x^2*Log[-1 + x^2])/(-1 + x) + 2*x*Log[-1 + x]*Log[-1 + x^2]`
- FriCAS `dep-elem-17`: `Log[-1 + x]/(1 + x) + Log[1 + x]/(-1 + x) + (2*x^2*Log[-1 + x]*Log[1 + x])/(1 + x) + (x^2*Log[1 + x]^2)/(-1 + x) + 2*x*Log[-1 + x]*Log[1 + x]^2 - (4*x^3*Log[2 +`
- FriCAS `dep-elem-21`: `(-2*x^2*Log[-1 + x])/(-1 + x) - (x^2*Log[-1 + x])/(1 + x) - 2*x*Log[-1 + x]^2 - (x^2*Log[1 + x])/(-1 + x) - (2*x^3*Log[1 + x])/(-1 + x^2) - 2*x*Log[-1 + x]*Log[`
