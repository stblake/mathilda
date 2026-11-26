# Integration shootout: 259-transcendental-decision (2026-11-26)

Re-tabulated from the run of 2026-11-26: 4 easy cases (easy.tsv: rational, rationalised by a substitution, or closed by DerivativeDivides) are left out.

555 cases, timeout 30 s per run, wall 2.7 min. Every answer of every system is judged by PossibleZeroQ[D[F, x] - f] in Mathilda (the other systems' answers translated into Mathilda first): zero is verified, nonzero is `wrong`, anything else -- undecided, untranslatable, or an answer text cut when the run was recorded -- is `unverified`. `elem`/`special`: verified elementary / special-function answer; `cert`: Mathilda proved non-elementary.

## Closed (elementary + special) by class

| run | T transcendental (555) | all |
|---|---|---|
| Mathilda:Automatic | 69% (382) | 69% (382) |
| Mathilda:RischTranscendental | 62% (343) | 62% (343) |
| FriCAS | 62% (345) | 62% (345) |

## Outcome counts

| run | elem | special | cert | uneval | timeout | error | unverified | wrong | median s (closed) |
|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 282 | 100 | 170 | 0 | 1 | 0 | 2 | 0 | 0.05 |
| Mathilda:RischTranscendental | 284 | 59 | 209 | 2 | 0 | 0 | 1 | 0 | 0.02 |
| FriCAS | 260 | 85 | 0 | 193 | 0 | 1 | 9 | 7 | 0.11 |

## PossibleZeroQ verdicts

Every answer, by PossibleZeroQ[D[F, x] - f] in Mathilda. `own ok`: of the verified answers, those the system's own three-point check also passed; `own disagrees`: answers where the own check and PossibleZeroQ disagree (one passes, the other does not).

| run | answers | zero (verified) | nonzero (wrong) | undecided | untranslatable | cut | own ok | own disagrees |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 384 | 382 | 0 | 2 | 0 | 0 | 382 | 0 |
| Mathilda:RischTranscendental | 344 | 343 | 0 | 1 | 0 | 0 | 343 | 0 |
| FriCAS | 361 | 345 | 7 | 9 | 0 | 0 | 345 | 7 |

## Closed by family (Mathilda = best of its methods)

| class | family | n | Mathilda | FriCAS |
|---|---|---|---|---|
| T | tanh-elem | 36 | 36 | 31 |
| T | tanh-rand | 31 | 11 | 8 |
| T | tan-elem | 33 | 33 | 33 |
| T | tan-rand | 33 | 7 | 7 |
| T | tani-elem | 38 | 38 | 37 |
| T | tani-rand | 28 | 12 | 13 |
| T | gauss-elem | 32 | 31 | 21 |
| T | gauss-rand | 34 | 24 | 21 |
| T | lim-elem | 29 | 29 | 28 |
| T | lim-rand | 39 | 3 | 3 |
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


## AUTO-GAP (2)

- `tanh-elem-39` AUTO-GAP Mathilda:RischTranscendental: `(-4*E^(x*Tanh[x^2])*x*Sech[x^2]^2)/(2 + Tanh[x^2]^2)^2 - (2*E^(x*Tanh[x^2])*x*Csch[x^2]^2)/(2 + Tanh[x^2]^2) + (E^(x*Tanh[x^2])*Coth[x^2]*(2*x^2*Sech[x^2]^2 + Tanh[x^2]))/(2 + Tanh[x^2]^2)`
- `tani-elem-26` AUTO-GAP Mathilda:RischTranscendental: `-((x*Sec[x]^2)/(I + Tan[x])^2) + (I + Tan[x])^(-1)`

## Only one system closes it

- Mathilda: 23

## Wrong answers from other systems (their bugs)

- FriCAS `tani-elem-35`: `(1 + I)*E^(2*Tan[1 + x]) + (2 + 2*I)*E^(2*Tan[1 + x])*x*Sec[1 + x]^2 - (E^Tan[1 + x]*Log[I + x]*Sec[1 + x]^2)/(2*(1 + Tan[1 + x])^2) + E^Tan[1 + x]/(2*(I + x)*(`
- FriCAS `prim-rand-22`: `x*Log[x]*LogIntegral[x]`
- FriCAS `prim-rand-38`: `Log[x]*LogIntegral[x]`
- FriCAS `dep-elem-9`: `-((x*Log[-1 + x])/(1 + x)) + Log[x]/(1 + x) + (2*x^3*Log[x])/(-1 + x^2) + Log[1 + x]/x - (x*Log[1 + x])/(-1 + x) - Log[-1 + x]*Log[1 + x] + x*Log[-1 + x^2] + 2*`
- FriCAS `dep-elem-12`: `(2*x^3*Log[-1 + x])/(-1 + x^2) + (2*x^2*Log[x])/(1 + x) + 2*x*Log[1 + x] + 4*x*Log[x]*Log[1 + x] + (x^2*Log[-1 + x^2])/(-1 + x) + 2*x*Log[-1 + x]*Log[-1 + x^2]`
- FriCAS `dep-elem-17`: `Log[-1 + x]/(1 + x) + Log[1 + x]/(-1 + x) + (2*x^2*Log[-1 + x]*Log[1 + x])/(1 + x) + (x^2*Log[1 + x]^2)/(-1 + x) + 2*x*Log[-1 + x]*Log[1 + x]^2 - (4*x^3*Log[2 +`
- FriCAS `dep-elem-21`: `(-2*x^2*Log[-1 + x])/(-1 + x) - (x^2*Log[-1 + x])/(1 + x) - 2*x*Log[-1 + x]^2 - (x^2*Log[1 + x])/(-1 + x) - (2*x^3*Log[1 + x])/(-1 + x^2) - 2*x*Log[-1 + x]*Log[`
