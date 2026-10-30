# Trager vs FriCAS stress report

Generated 2026-10-30 11:08 by `tools/trager_fricas_stress.py` over `tests/data/trager_stress_corpus.m` (378 cases, timeout 60 s per case and system, 6 jobs, 802 s wall).

Classes: see the docstring of the tool.  WRONG and DISAGREE are hard failures; FRICAS-ONLY ratchets against `BASELINE`.

## Summary

| class | Trager | TragerKauers |
|---|---:|---:|
| AGREE-ELEM | 180 | 122 |
| AGREE-NONELEM | 132 | 109 |
| MATHILDA-ONLY | 3 | 1 |
| FRICAS-ONLY | 0 | 58 |
| DISAGREE | 0 | 0 |
| WRONG | 0 | 0 |
| NONELEM-UNCONFIRMED | 63 | 13 |
| NEITHER | 0 | 70 |
| TIMEOUT | 0 | 5 |
| ERROR | 0 | 0 |

FriCAS outcomes: BAD 1, ELEM 180, ERROR 19, NONELEM 132, TIMEOUT 46 (BAD = FriCAS antiderivative failed its own derivative check).

## Per family

| family | n | Trager | TragerKauers |
|---|---:|---|---|
| binomial-xk | 45 | AGREE-ELEM 21, AGREE-NONELEM 24 | AGREE-ELEM 19, AGREE-NONELEM 24, FRICAS-ONLY 2 |
| chebyshev | 40 | AGREE-ELEM 20, AGREE-NONELEM 20 | AGREE-ELEM 14, AGREE-NONELEM 20, FRICAS-ONLY 6 |
| deriv-alg | 20 | AGREE-ELEM 20 | AGREE-ELEM 20 |
| deriv-log | 25 | AGREE-ELEM 25 | AGREE-ELEM 5, FRICAS-ONLY 20 |
| even-substitution | 6 | AGREE-ELEM 6 | AGREE-ELEM 6 |
| genus-hyper | 20 | AGREE-ELEM 4, AGREE-NONELEM 16 | AGREE-ELEM 4, AGREE-NONELEM 16 |
| genus-super | 10 | AGREE-ELEM 1, AGREE-NONELEM 7, NONELEM-UNCONFIRMED 2 | AGREE-ELEM 1, AGREE-NONELEM 7, NEITHER 2 |
| goursat | 23 | AGREE-ELEM 13, AGREE-NONELEM 5, MATHILDA-ONLY 1, NONELEM-UNCONFIRMED 4 | AGREE-ELEM 8, AGREE-NONELEM 1, FRICAS-ONLY 5, NEITHER 8, TIMEOUT 1 |
| pseudo-elliptic-recip | 11 | AGREE-ELEM 11 | AGREE-ELEM 2, FRICAS-ONLY 9 |
| random | 60 | AGREE-ELEM 16, AGREE-NONELEM 25, MATHILDA-ONLY 1, NONELEM-UNCONFIRMED 18 | AGREE-ELEM 11, AGREE-NONELEM 22, FRICAS-ONLY 5, MATHILDA-ONLY 1, NEITHER 12, NONELEM-UNCONFIRMED 6, TIMEOUT 3 |
| rational-radicand | 20 | AGREE-ELEM 12, AGREE-NONELEM 8 | AGREE-ELEM 6, AGREE-NONELEM 8, FRICAS-ONLY 6 |
| rational-times-radical | 50 | AGREE-ELEM 17, AGREE-NONELEM 11, MATHILDA-ONLY 1, NONELEM-UNCONFIRMED 21 | AGREE-ELEM 12, AGREE-NONELEM 11, FRICAS-ONLY 5, NEITHER 15, NONELEM-UNCONFIRMED 7 |
| third-kind | 30 | AGREE-ELEM 14, AGREE-NONELEM 16 | AGREE-ELEM 14, NEITHER 16 |
| torsion-heavy | 18 | NONELEM-UNCONFIRMED 18 | NEITHER 17, TIMEOUT 1 |

## Timing (seconds; Mathilda: AbsoluteTiming of the worker; FriCAS: its own `Time:` total for integrate)

| system | median | mean | max | n |
|---|---:|---:|---:|---:|
| Mathilda Trager | 0.036 | 0.282 | 21.354 | 378 |
| Mathilda TragerKauers | 0.034 | 4.125 | 60.800 | 378 |
| FriCAS | 0.090 | 7.814 | 60.068 | 359 |

### Mathilda > 10x slower than FriCAS (28)

| id | method | Mathilda s | FriCAS s |
|---|---|---:|---:|
| ratrad-045 | Trager | 3.515 | 0.10 |
| ratrad-047 | TragerKauers | 2.392 | 0.10 |
| recip-004 | TragerKauers | 1.664 | 0.16 |
| recip-005 | TragerKauers | 2.040 | 0.15 |
| recip-009 | TragerKauers | 2.561 | 0.15 |
| recip-010 | TragerKauers | 60.800 | 0.17 |
| recip-011 | TragerKauers | 2.064 | 0.15 |
| ratradicand-009 | TragerKauers | 7.909 | 0.11 |
| deriv-log-03 | TragerKauers | 60.289 | 0.12 |
| deriv-log-04 | TragerKauers | 60.279 | 0.15 |
| deriv-log-05 | TragerKauers | 60.234 | 0.49 |
| deriv-log-06 | TragerKauers | 60.322 | 0.25 |
| deriv-log-07 | TragerKauers | 60.223 | 0.26 |
| deriv-log-08 | TragerKauers | 60.226 | 0.48 |
| deriv-log-11 | TragerKauers | 60.249 | 0.11 |
| deriv-log-12 | TragerKauers | 2.059 | 0.13 |
| deriv-log-13 | TragerKauers | 60.250 | 0.54 |
| deriv-log-15 | TragerKauers | 60.258 | 0.37 |
| deriv-log-16 | TragerKauers | 60.273 | 0.17 |
| deriv-log-18 | TragerKauers | 60.246 | 0.18 |
| deriv-log-19 | TragerKauers | 60.254 | 0.36 |
| deriv-log-20 | TragerKauers | 60.243 | 0.22 |
| deriv-log-21 | TragerKauers | 60.298 | 0.23 |
| deriv-log-22 | TragerKauers | 60.278 | 0.27 |
| deriv-log-23 | TragerKauers | 60.227 | 0.33 |
| deriv-log-24 | TragerKauers | 60.260 | 0.28 |
| deriv-log-25 | TragerKauers | 60.242 | 0.49 |
| random-053 | TragerKauers | 1.693 | 0.12 |

## Hard failures and gaps

### WRONG (0)


### DISAGREE (0)


### FRICAS-ONLY (58)

- `binom-029` [TragerKauers] `x^3/(x^3 + 2)^(1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*3^(1/2)*log((((x^3+2)^(1/3))^2+x*(x^3+2)^(1/3)+x^2)/(x^2))+(2*3^(1/2)*log(((x^3+2)^(1/3)+(-1)*x)/x)+(6*atan((2*3^(1/2)*(x^3+2)^(1/3)+x*3^(1/2))/(3*x))+3*x*3^(1/2)*((x^3+2)^(1/3))^2)))/(9*3^(1/2)...`
- `binom-032` [TragerKauers] `1/(x^4 + 1)^(1/4)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(log(((x^4+1)^(1/4)+x)/x)+((-1)*log(((x^4+1)^(1/4)+(-1)*x)/x)+(-2)*atan(((x^4+1)^(1/4))/x)))/4`
- `ratrad-006` [TragerKauers] `1/x*(x^3 + 1)^(2/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*3^(1/2)*log(((x^3+1)^(1/3))^2+(x^3+1)^(1/3)+1)+(2*3^(1/2)*log((x^3+1)^(1/3)+(-1))+(6*atan((2*3^(1/2)*(x^3+1)^(1/3)+3^(1/2))/3)+3*3^(1/2)*((x^3+1)^(1/3))^2)))/(6*3^(1/2))`
- `ratrad-025` [TragerKauers] `(2*x - 1)/(x + 1)*(x + 2)^(1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(3*log(((x+2)^(1/3))^2+(x+2)^(1/3)+1)+((-6)*log((x+2)^(1/3)+(-1))+(6*3^(1/2)*atan((2*(x+2)^(1/3)+1)/(3^(1/2)))+(3*x+(-12))*(x+2)^(1/3))))/2`
- `ratrad-033` [TragerKauers] `x/(x - 2)*(x + 2)^(-1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*2^(1/3)*log(((x+2)^(1/3))^2+(2^(1/3))^2*(x+2)^(1/3)+2*2^(1/3))+(2*2^(1/3)*log((x+2)^(1/3)+(-1)*(2^(1/3))^2)+(2*2^(1/3)*3^(1/2)*atan((2*(x+2)^(1/3)+(2^(1/3))^2)/((2^(1/3))^2*3^(1/2)))+3*((x+2)^(1...`
- `ratrad-045` [TragerKauers] `1/x*(x^2 + 1)^(-1/6)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*log(((x^2+1)^(1/6))^2+(x^2+1)^(1/6)+1)+(log(((x^2+1)^(1/6))^2+(-1)*(x^2+1)^(1/6)+1)+((-2)*log((x^2+1)^(1/6)+1)+(2*log((x^2+1)^(1/6)+(-1))+(2*3^(1/2)*atan((2*(x^2+1)^(1/6)+1)/(3^(1/2)))+2*3^(1/2)...`
- `ratrad-047` [TragerKauers] `(2*x - 1)/x*(x + 2)^(-1/4)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(12*atan(2/((2*(2^(1/4))^2*((x+2)^(1/4))^2+4)^(1/2)+(2^(1/4))^3*(x+2)^(1/4)))+(3*log((2^(1/4))^3*(x+2)^(1/4)+2)+((-3)*log((2^(1/4))^3*(x+2)^(1/4)+(-2))+8*2^(1/4)*((x+2)^(1/4))^3)))/(3*2^(1/4))`
- `cheb-007` [TragerKauers] `x^(-1)*Sqrt[-x^3 + 3]` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(3*log((6*((-1)*x^3+3)^(1/2)+(x^3+(-6))*3^(1/2))/(x^3))+2*3^(1/2)*((-1)*x^3+3)^(1/2))/(3*3^(1/2))`
- `cheb-013` [TragerKauers] `x^(-1)*Sqrt[x^4 + 1]` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*log((x^4+1)^(1/2)+1)+(log((x^4+1)^(1/2)+(-1))+2*(x^4+1)^(1/2)))/4`
- `cheb-014` [TragerKauers] `x^3*(2*x^3 - 1)^(2/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*3^(1/2)*log((2^(1/3)*((2*x^3+(-1))^(1/3))^2+x*(2^(1/3))^2*(2*x^3+(-1))^(1/3)+2*x^2)/(x^2))+(2*3^(1/2)*log(((2^(1/3))^2*(2*x^3+(-1))^(1/3)+(-2)*x)/x)+(6*atan(((2^(1/3))^2*3^(1/2)*(2*x^3+(-1))^(1/...`
- `cheb-020` [TragerKauers] `x^(-1)*(2*x^4 - 1)^(2/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*log(((2*x^4+(-1))^(1/3))^2+(-1)*(2*x^4+(-1))^(1/3)+1)+(2*log((2*x^4+(-1))^(1/3)+1)+((-2)*3^(1/2)*atan((2*(2*x^4+(-1))^(1/3)+(-1))/(3^(1/2)))+3*((2*x^4+(-1))^(1/3))^2)))/8`
- `cheb-027` [TragerKauers] `x^(-1)*(-x^2 + 3)^(-1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*log(3^(1/3)*(((-1)*x^2+3)^(1/3))^2+(3^(1/3))^2*((-1)*x^2+3)^(1/3)+3)+(2*log((3^(1/3))^2*((-1)*x^2+3)^(1/3)+(-3))+2*3^(1/2)*atan((2*(3^(1/3))^2*((-1)*x^2+3)^(1/3)+3)/(3*3^(1/2)))))/(4*3^(1/3))`
- `cheb-038` [TragerKauers] `x^(-3)*(-x^4 + 3)^(3/2)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(9*x^2*atan((((-1)*x^4+3)^(1/2))/(x^2))+((-1)*x^4+(-6))*((-1)*x^4+3)^(1/2))/(4*x^2)`
- `goursat-003` [TragerKauers] `x/Sqrt[x^4 + 10*x^2 - 96*x - 71]` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `log((x^6+15*x^4+(-80)*x^3+27*x^2+(-528)*x+781)*(x^4+10*x^2+(-96)*x+(-71))^(1/2)+(x^8+20*x^6+(-128)*x^5+54*x^4+(-1408)*x^3+3124*x^2+10001))/8`
- `goursat-006` [TragerKauers] `1/(x^3 + 1)^(1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(3^(1/2)*log((((x^3+1)^(1/3))^2+x*(x^3+1)^(1/3)+x^2)/(x^2))+((-2)*3^(1/2)*log(((x^3+1)^(1/3)+(-1)*x)/x)+(-6)*atan((2*3^(1/2)*(x^3+1)^(1/3)+x*3^(1/2))/(3*x))))/(6*3^(1/2))`
- `goursat-010` [TragerKauers] `1/((x^4 + 1)^(1/4))` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(log(((x^4+1)^(1/4)+x)/x)+((-1)*log(((x^4+1)^(1/4)+(-1)*x)/x)+(-2)*atan(((x^4+1)^(1/4))/x)))/4`
- `goursat-012` [TragerKauers] `1/(x*(x^4 + 1)^(1/4))` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*log((x^4+1)^(1/4)+1)+(log((x^4+1)^(1/4)+(-1))+2*atan((x^4+1)^(1/4))))/4`
- `goursat-022` [TragerKauers] `(x - 1)/((x + 1)*Sqrt[x^3 + x^2 + x])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `atan((x^2+1)/(2*(x^3+x^2+x)^(1/2)))`
- `recip-002` [TragerKauers] `(x^2 - 1)/((x^2 + x + 1)*Sqrt[x^4 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `atan(((x^4+1)^(1/2))/(x^2+2*x+1))`
- `recip-003` [TragerKauers] `(x^2 - 1)/((x^2 + 2*x + 1)*Sqrt[x^4 + x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `log((((-12)*x^2+(-6)*x+(-12))*(x^4+x^2+1)^(1/2)+((-7)*x^4+(-4)*x^3+(-12)*x^2+(-4)*x+(-7))*3^(1/2))/(x^4+4*x^3+6*x^2+4*x+1))/(2*3^(1/2))`
- `recip-004` [TragerKauers] `(x^2 - 1)/((x^2 + 3*x + 1)*Sqrt[x^4 - x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(2^(1/2)*log((((-6)*x^2+(-6)*x+(-6))*2^(1/2)*(x^4+(-1)*x^2+1)^(1/2)+((-5)*x^4+(-6)*x^3+(-7)*x^2+(-6)*x+(-5))*3^(1/2))/(x^4+6*x^3+11*x^2+6*x+1)))/(4*3^(1/2))`
- `recip-005` [TragerKauers] `(x^2 - 1)/((x^2 - x + 1)*Sqrt[x^4 + 3*x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(2^(1/2)*log((((-2)*x^2+(-2)*x+(-2))*2^(1/2)*(x^4+3*x^2+1)^(1/2)+(3*x^4+2*x^3+9*x^2+2*x+3))/(x^4+(-2)*x^3+3*x^2+(-2)*x+1)))/4`
- `recip-007` [TragerKauers] `(x^2 - 1)/((x^2 + x + 1)*Sqrt[x^4 + 3*x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(2^(1/2)*log(((2*x^2+(-2)*x+2)*2^(1/2)*(x^4+3*x^2+1)^(1/2)+(3*x^4+(-2)*x^3+9*x^2+(-2)*x+3))/(x^4+2*x^3+3*x^2+2*x+1)))/4`
- `recip-008` [TragerKauers] `(x^2 + 1)/((x^2 + x - 1)*Sqrt[x^4 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `log(((3*x^2+(-6)*x+(-3))*(x^4+1)^(1/2)+(2*x^4+(-2)*x^3+x^2+2*x+2)*3^(1/2))/(x^4+2*x^3+(-1)*x^2+(-2)*x+1))/(2*3^(1/2))`
- `recip-009` [TragerKauers] `(x^2 + 1)/((x^2 + 2*x - 1)*Sqrt[x^4 + x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `log(((28*x^2+(-42)*x+(-28))*(x^4+x^2+1)^(1/2)+(11*x^4+(-12)*x^3+8*x^2+12*x+11)*7^(1/2))/(x^4+4*x^3+2*x^2+(-4)*x+1))/(2*7^(1/2))`
- `recip-010` [TragerKauers] `(x^2 + 1)/((x^2 + 3*x - 1)*Sqrt[x^4 + 3*x^2 + 1])` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `(2^(1/2)*log(((42*x^2+(-70)*x+(-42))*2^(1/2)*(x^4+3*x^2+1)^(1/2)+(23*x^4+(-30)*x^3+49*x^2+30*x+23)*7^(1/2))/(x^4+6*x^3+7*x^2+(-6)*x+1)))/(4*7^(1/2))`
- `recip-011` [TragerKauers] `(x^2 + 1)/((x^2 - x - 1)*Sqrt[x^4 - x^2 + 1])` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(2^(1/2)*log((((-2)*x^2+(-2)*x+2)*2^(1/2)*(x^4+(-1)*x^2+1)^(1/2)+(3*x^4+2*x^3+(-3)*x^2+(-2)*x+3))/(x^4+(-2)*x^3+(-1)*x^2+2*x+1)))/4`
- `ratradicand-002` [TragerKauers] `((x + 1)/(x - 1))^(1/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(3^(1/2)*log((((x+1)/(x+(-1)))^(1/3))^2+((x+1)/(x+(-1)))^(1/3)+1)+((-2)*3^(1/2)*log(((x+1)/(x+(-1)))^(1/3)+(-1))+(6*atan((2*3^(1/2)*((x+1)/(x+(-1)))^(1/3)+3^(1/2))/3)+(3*x+(-3))*3^(1/2)*((x+1)/(x+(-1)...`
- `ratradicand-004` [TragerKauers] `((1 - x)/(1 + x))^(1/4)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(4*2^(1/2)*atan(1/((2*((((-1)*x+1)/(x+1))^(1/4))^2+(-2)*2^(1/2)*(((-1)*x+1)/(x+1))^(1/4)+2)^(1/2)+(2^(1/2)*(((-1)*x+1)/(x+1))^(1/4)+(-1))))+(4*2^(1/2)*atan(1/((2*((((-1)*x+1)/(x+1))^(1/4))^2+2*2^(1/2)...`
- `ratradicand-007` [TragerKauers] `((x + 2)/(x - 3))^(2/3)/x` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*4^(1/3)*log(6*9^(1/3)*(((x+2)/(x+(-3)))^(1/3))^2+(-1)*(4^(1/3))^2*(9^(1/3))^2*((x+2)/(x+(-3)))^(1/3)+6*4^(1/3))+(9^(1/3)*log((((x+2)/(x+(-3)))^(1/3))^2+((x+2)/(x+(-3)))^(1/3)+1)+(2*4^(1/3)*log(2...`
- `ratradicand-008` [TragerKauers] `1/(x^2*((x + 1)/(x - 1))^(1/5))` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(((-1)*x*((-75)*rootOf((125*rootOf((625*%%H0^4+250*%%H0^3+100*%%H0^2+40*%%H0+16)/625,%%H0)^3+(125*%%H1+50)*rootOf((625*%%H0^4+250*%%H0^3+100*%%H0^2+40*%%H0+16)/625,%%H0)^2+(125*%%H1^2+50*%%H1+20)*root...`
- `ratradicand-009` [TragerKauers] `((2*x + 1)/(x + 3))^(1/6)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(20*3^(1/2)*atan((3^(1/2))/(((32^(1/6))^2*(((2*x+1)/(x+3))^(1/6))^2+(-2)*32^(1/6)*((2*x+1)/(x+3))^(1/6)+4)^(1/2)+(32^(1/6)*((2*x+1)/(x+3))^(1/6)+(-1))))+(20*3^(1/2)*atan((3^(1/2))/(((32^(1/6))^2*(((2*...`
- `ratradicand-014` [TragerKauers] `((x - 2)/(x + 2))^(3/4)/(x - 1)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(4*atan(9/((54*(27^(1/4))^2*(((x+(-2))/(x+2))^(1/4))^2+(-18)*2^(1/2)*(27^(1/4))^3*((x+(-2))/(x+2))^(1/4)+162)^(1/2)+(2^(1/2)*(27^(1/4))^3*((x+(-2))/(x+2))^(1/4)+(-9))))+(4*atan(9/((54*(27^(1/4))^2*(((...`
- `deriv-log-01` [TragerKauers] `(-27*x^5 - 45*x^4 + 126*x^3 + 9*x^2 - 21*x)/(2*(x^3 + x^2 - 3*x + 1)*(9*x^6 + 18*x^5 + 27*x^4 - 3*x^3 - 12*x^2 + 15*x - 5))*Sqrt[x^3 + x^2 - 3*x + 1]` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: logarithmic part incomplete"}`; FriCAS: ELEM `log((((-18)*x^3+(-18)*x^2+(-18)*x+12)*(x^3+x^2+(-3)*x+1)^(1/2)+((-9)*x^6+(-18)*x^5+(-27)*x^4+(-15)*x^3+(-6)*x^2+39*x+(-13)))/(9*x^6+18*x^5+27*x^4+(-3)*x^3+(-12)*x^2+15*x+(-5)))/2`
- `deriv-log-03` [TragerKauers] `(9*x^5 + 3*x^4 + 9*x^3 - 18*x^2 + 6*x - 9)/(2*(x^3 + x)*(x^6 + 2*x^5 + 5*x^4 + x^3 + 10*x^2 + 3*x + 9))*Sqrt[x^3 + x]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((6*x^3+6*x^2+12*x+18)*(x^3+x)^(1/2)+((-1)*x^6+(-2)*x^5+(-5)*x^4+(-19)*x^3+(-10)*x^2+(-21)*x+(-9)))/(x^6+2*x^5+5*x^4+x^3+10*x^2+3*x+9))/2`
- `deriv-log-04` [TragerKauers] `(-2*x^6 - 9*x^5 + 81*x^4 - 74*x^3 + 56*x^2 - 57*x + 13)/(2*(x^3 - 3*x^2 + x - 1)*(4*x^6 - 13*x^5 + 14*x^4 + 20*x^2 + x + 10))*Sqrt[x^3 - 3*x^2 + x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((4*x^4+6*x^3+(-14)*x^2+14*x+6)*(x^3+(-3)*x^2+x+(-1))^(1/2)+(4*x^6+(-11)*x^5+20*x^4+(-16)*x^3+(-24)*x^2+7*x+(-8)))/(4*x^6+(-13)*x^5+14*x^4+20*x^2+x+10))/2`
- `deriv-log-05` [TragerKauers] `(4*x^9 - 26*x^8 + 23*x^7 - 19*x^6 + 71*x^5 + 82*x^4 + 25*x^3 + 12*x^2 - 18*x + 15)/(2*(x^6 - x^5 - 3*x^4 - x^3 - 2*x^2 + x - 1)*(-4*x^8 + 8*x^7 + 8*x^6 - 13*x^5 + 16*x^4 - 17*x^3 + 28*x^2 - 5*x + 10))*Sqrt[x^6 - x^5 - 3*x^4 - x^3 - 2*x^2 + x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-4)*x^4+14*x^3+(-6)*x^2+12*x+(-6))*(x^6+(-1)*x^5+(-3)*x^4+(-1)*x^3+(-2)*x^2+x+(-1))^(1/2)+((-4)*x^8+8*x^7+6*x^6+(-1)*x^5+(-2)*x^4+(-5)*x^3+(-8)*x^2+(-5)*x+(-8)))/(4*x^8+(-8)*x^7+(-8)*x^6+13*x^5...`
- `deriv-log-06` [TragerKauers] `(12*x^5 - 36*x^4 + 6*x^3 - 33*x^2 + 54*x - 36)/(2*(x^4 + x^2 - 3*x + 3)*(-9*x^6 - 9*x^4 + 27*x^3 - 26*x^2 - 4*x + 4))*Sqrt[x^4 + x^2 - 3*x + 3]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((6*x^2+(-12)*x)*(x^4+x^2+(-3)*x+3)^(1/2)+(9*x^6+9*x^4+(-27)*x^3+28*x^2+(-4)*x+4))/(9*x^6+9*x^4+(-27)*x^3+26*x^2+4*x+(-4)))/2`
- `deriv-log-07` [TragerKauers] `(-15*x^6 + 33*x^5 - 15*x^4 + 69*x^3 + 16*x^2 - 6*x - 6)/(2*(x^4 - 3*x^3 + 2*x^2 + 2*x + 3)*(8*x^6 + 21*x^5 - 11*x^4 - 14*x^3 + 12*x^2 - 6*x + 1))*Sqrt[x^4 - 3*x^3 + 2*x^2 + 2*x + 3]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((6*x^4+6*x^3+(-6)*x^2+2*x)*(x^4+(-3)*x^3+2*x^2+2*x+3)^(1/2)+((-10)*x^6+(-15)*x^5+7*x^4+10*x^3+(-18)*x^2+6*x+(-1)))/(8*x^6+21*x^5+(-11)*x^4+(-14)*x^3+12*x^2+(-6)*x+1))/2`
- `deriv-log-08` [TragerKauers] `(4*x^8 - 3*x^7 - 31*x^6 - 15*x^5 - 27*x^4 + 4*x^3 + 18*x^2 - 36*x + 24)/(2*(x^6 + x^5 + 3*x^3 - 2*x^2 - 2*x + 2)*(-x^8 - 7*x^7 - 15*x^6 - 12*x^5 - 15*x^4 - 17*x^3 + 32*x^2 + 6*x - 18))*Sqrt[x^6 + x^5 + 3*x^3 - 2*x^2 - 2*x + 2]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-2)*x^3+(-2)*x^2+12*x)*(x^6+x^5+3*x^3+(-2)*x^2+(-2)*x+2)^(1/2)+((-1)*x^8+(-7)*x^7+(-15)*x^6+(-12)*x^5+(-17)*x^4+(-9)*x^3+24*x^2+6*x+(-18)))/(x^8+7*x^7+15*x^6+12*x^5+15*x^4+17*x^3+(-32)*x^2+(-6)...`
- `deriv-log-11` [TragerKauers] `(-4*x^6 - 4*x^5 + 42*x^4 + 32*x^3 - 14*x^2)/(2*(-2*x^3 - 3*x^2 + x - 1)*(x^6 + 12*x^5 + 22*x^4 + 8*x^3 + 13*x^2))*Sqrt[-2*x^3 - 3*x^2 + x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((4*x^2+8*x+12)*((-2)*x^3+(-3)*x^2+x+(-1))^(1/2)+((-1)*x^4+4*x^3+2*x^2+(-16)*x+(-5)))/(x^4+12*x^3+22*x^2+8*x+13))/2`
- `deriv-log-12` [TragerKauers] `(-8*x^6 - 4*x^5 - 4*x^3)/(2*(x^4 + x^3 + 2*x^2 - x)*(-4*x^6 - 4*x^5 - 4*x^4 + 4*x^3))*Sqrt[x^4 + x^3 + 2*x^2 - x]` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `log((2*(x^4+x^3+2*x^2+(-1)*x)^(1/2)+((-1)*x^3+(-1)*x^2+(-3)*x+1))/(x^3+x^2+x+(-1)))/2`
- `deriv-log-13` [TragerKauers] `(-32*x^8 - 44*x^7 + 32*x^6 + 147*x^5 + 72*x^4 - 45*x^3 + 22*x^2 + 46*x - 2)/(2*(-2*x^6 + 3*x^5 + 3*x^4 - 3*x^3 + x^2 + 2*x - 3)*(8*x^8 - 4*x^7 - 22*x^6 - 3*x^5 + 9*x^4 + 3*x^3 + 20*x^2 + 22*x + 7))*Sqrt[-2*x^6 + 3*x^5 + 3*x^4 - 3*x^3 + x^2 + 2*x - 3]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((8*x^3+16*x^2+14*x+4)*((-2)*x^6+3*x^5+3*x^4+(-3)*x^3+x^2+2*x+(-3))^(1/2)+((-8)*x^8+4*x^7+22*x^6+3*x^5+(-1)*x^4+21*x^3+14*x^2+2*x+1))/(8*x^8+(-4)*x^7+(-22)*x^6+(-3)*x^5+9*x^4+3*x^3+20*x^2+22*x+7))...`
- `deriv-log-15` [TragerKauers] `(-6*x^7 - 21*x^6 - 18*x^4 - 48*x^3 + 12*x - 3)/(2*(x^6 - x^5 + 3*x^4 + x^3 - x^2 + x)*(-9*x^6 + 9*x^5 - 26*x^4 - 5*x^3 + 15*x^2 - 5*x + 1))*Sqrt[x^6 - x^5 + 3*x^4 + x^3 - x^2 + x]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((6*x^2+12*x+6)*(x^6+(-1)*x^5+3*x^4+x^3+(-1)*x^2+x)^(1/2)+((-9)*x^6+9*x^5+(-28)*x^4+(-13)*x^3+3*x^2+(-13)*x+(-1)))/(9*x^6+(-9)*x^5+26*x^4+5*x^3+(-15)*x^2+5*x+(-1)))/2`
- `deriv-log-16` [TragerKauers] `(-9*x^5 + 27*x^4 - 43*x^3 + 27*x^2 + 2*x - 8)/(2*(x^3 - 2*x^2 + 2*x - 1)*(9*x^6 - 18*x^5 + 3*x^4 - 13*x^3 + 21*x^2 + 4*x + 10))*Sqrt[x^3 - 2*x^2 + 2*x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-6)*x^3+6*x^2+2*x+6)*(x^3+(-2)*x^2+2*x+(-1))^(1/2)+((-9)*x^6+18*x^5+(-3)*x^4+11*x^3+(-17)*x^2+(-8)*x+(-8)))/(9*x^6+(-18)*x^5+3*x^4+(-13)*x^3+21*x^2+4*x+10))/2`
- `deriv-log-18` [TragerKauers] `(-18*x^4 + 117*x^2 + 126*x - 18)/(2*(x^4 + 3*x^3 + 3*x^2 - 2*x - 3)*(-5*x^4 - 27*x^3 - 15*x^2 + 18*x + 36))*Sqrt[x^4 + 3*x^3 + 3*x^2 - 2*x - 3]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((12*x^2+18)*(x^4+3*x^3+3*x^2+(-2)*x+(-3))^(1/2)+(13*x^4+27*x^3+39*x^2+(-18)*x+(-18)))/(5*x^4+27*x^3+15*x^2+(-18)*x+(-36)))/2`
- `deriv-log-19` [TragerKauers] `(-9*x^8 - 3*x^7 - 12*x^6 + 57*x^5 + 51*x^4 - 123*x^3 + 60*x^2 + 18*x - 11)/(2*(x^5 + x^4 + x^3 + 3*x^2 + 3*x - 2)*(-9*x^7 + 6*x^6 - 4*x^5 - 16*x^4 - 16*x^3 + 34*x^2 - 17*x + 3))*Sqrt[x^5 + x^4 + x^3 + 3*x^2 + 3*x - 2]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((18*x^4+(-6)*x^3+6*x^2+(-8)*x+2)*(x^5+x^4+x^3+3*x^2+3*x+(-2))^(1/2)+((-9)*x^7+(-12)*x^6+(-4)*x^5+(-28)*x^4+(-4)*x^3+32*x^2+(-13)*x+1))/(9*x^7+(-6)*x^6+4*x^5+16*x^4+16*x^3+(-34)*x^2+17*x+(-3)))/2`
- `deriv-log-20` [TragerKauers] `(9*x^7 - 13*x^6 - 5*x^5 - 18*x^4 + 3*x^3 + 5*x^2 - 3*x - 2)/(2*(x^5 - 3*x^4 - x - 1)*(-x^7 + 3*x^6 + 9*x^4 + 7*x^3 + 8*x^2 + 2*x + 1))*Sqrt[x^5 - 3*x^4 - x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-6)*x^3+(-2)*x^2+(-2)*x)*(x^5+(-3)*x^4+(-1)*x+(-1))^(1/2)+((-1)*x^7+3*x^6+(-9)*x^4+(-5)*x^3+(-6)*x^2+(-2)*x+(-1)))/(x^7+(-3)*x^6+(-9)*x^4+(-7)*x^3+(-8)*x^2+(-2)*x+(-1)))/2`
- `deriv-log-21` [TragerKauers] `(-4*x^6 + 12*x^5 - 18*x^4 - 29*x^3 - 16*x^2 - 2*x - 5)/(2*(2*x^4 + 3*x^2 + 3*x + 2)*(2*x^6 + 8*x^4 + x^3 - x - 1))*Sqrt[2*x^4 + 3*x^2 + 3*x + 2]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-4)*x^4+(-6)*x^3+(-8)*x^2+(-8)*x+(-2))*(2*x^4+3*x^2+3*x+2)^(1/2)+((-6)*x^6+(-8)*x^5+(-18)*x^4+(-19)*x^3+(-22)*x^2+(-13)*x+(-3)))/(2*x^6+8*x^4+x^3+(-1)*x+(-1)))/2`
- `deriv-log-22` [TragerKauers] `(6*x^6 + 6*x^5 - 2*x^4 + 7*x^3 + 9*x^2 + 4*x)/(2*(x^4 + 2*x^3 - 2*x^2 - x - 1)*(-9*x^6 - 24*x^5 + 6*x^4 + 19*x^3 + 17*x^2 + 7*x + 1))*Sqrt[x^4 + 2*x^3 - 2*x^2 - x - 1]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((6*x^3+2*x^2)*(x^4+2*x^3+(-2)*x^2+(-1)*x+(-1))^(1/2)+(9*x^6+24*x^5+(-4)*x^4+(-19)*x^3+(-17)*x^2+(-7)*x+(-1)))/(9*x^6+24*x^5+(-6)*x^4+(-19)*x^3+(-17)*x^2+(-7)*x+(-1)))/2`
- `deriv-log-23` [TragerKauers] `(15*x^6 - 93*x^5 + 173*x^4 - 183*x^3 + 21*x^2 - 13*x + 6)/(2*(-x^5 + 2*x^4 - 3*x^3 - x^2 - x)*(9*x^7 - 30*x^6 + 55*x^5 - 35*x^4 + 9*x^3 - 7*x^2 - 2*x + 9))*Sqrt[-x^5 + 2*x^4 - 3*x^3 - x^2 - x]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-6)*x^2+22*x+(-12))*((-1)*x^5+2*x^4+(-3)*x^3+(-1)*x^2+(-1)*x)^(1/2)+((-9)*x^7+30*x^6+(-55)*x^5+35*x^4+(-9)*x^3+9*x^2+(-10)*x+9))/(9*x^7+(-30)*x^6+55*x^5+(-35)*x^4+9*x^3+(-7)*x^2+(-2)*x+9))/2`
- `deriv-log-24` [TragerKauers] `(12*x^8 - 40*x^7 - 104*x^6 - 30*x^5 + 45*x^4 + 21*x^3 + 3*x^2 - 3*x)/(2*(2*x^5 + x^4 - x^3 - x^2 - x)*(-18*x^7 - 17*x^6 - 7*x^5 + 6*x^4 + 28*x^3 + 16*x^2 + x))*Sqrt[2*x^5 + x^4 - x^3 - x^2 - x]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log(((12*x^3+(-8)*x^2+(-22)*x+(-6))*(2*x^5+x^4+(-1)*x^3+(-1)*x^2+(-1)*x)^(1/2)+(18*x^6+25*x^5+(-9)*x^4+(-22)*x^3+(-4)*x^2+2*x+(-1)))/(18*x^6+17*x^5+7*x^4+(-6)*x^3+(-28)*x^2+(-16)*x+(-1)))/2`
- `deriv-log-25` [TragerKauers] `(-24*x^8 + 5*x^6 - 16*x^5 + 45*x^4 - 8*x^3 - 71*x^2 + 120*x - 3)/(2*(x^6 - x^5 - x^4 + 2*x^3 + x^2 - 3*x + 3)*(-4*x^8 - 8*x^7 + 7*x^6 + 13*x^5 - 10*x^4 - 18*x^3 + 21*x^2 - 9*x - 26))*Sqrt[x^6 - x^5 - x^4 + 2*x^3 + x^2 - 3*x + 3]` -- Mathilda: TIMEOUT `$Aborted`; FriCAS: ELEM `log((((-12)*x^3+(-18)*x^2+(-4)*x+(-6))*(x^6+(-1)*x^5+(-1)*x^4+2*x^3+x^2+(-3)*x+3)^(1/2)+(4*x^8+8*x^7+(-7)*x^6+(-13)*x^5+28*x^4+18*x^3+(-9)*x^2+9*x+28))/(4*x^8+8*x^7+(-7)*x^6+(-13)*x^5+10*x^4+18*x^3+(-...`
- `random-009` [TragerKauers] `2/(-3*x)*(x - 1)^(-1/4)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(8*atan(1/((2*((x+(-1))^(1/4))^2+(-2)*2^(1/2)*(x+(-1))^(1/4)+2)^(1/2)+(2^(1/2)*(x+(-1))^(1/4)+(-1))))+(8*atan(1/((2*((x+(-1))^(1/4))^2+2*2^(1/2)*(x+(-1))^(1/4)+2)^(1/2)+(2^(1/2)*(x+(-1))^(1/4)+1)))+(2...`
- `random-027` [TragerKauers] `(-2)/(-x^2 + 2*x - 3)*(3*x + 2)^(2/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(((-4)*43^(1/3)*sin((2*atan(7/(30*2^(1/2)+(-43))))/3)+(-4)*3^(1/2)*43^(1/3)*cos((2*atan(7/(30*2^(1/2)+(-43))))/3))*atan(((12*3^(1/2)+10*2^(1/2))*(43^(1/3))^2*cos((2*atan(7/(30*2^(1/2)+(-43))))/3)*sin(...`
- `random-030` [TragerKauers] `(-3*x^2 + 3*x - 3)/(2*x^2 - x - 1)*(x + 3)^(1/4)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `(140*5^(1/4)*atan((5^(1/4))/(((2^(1/4))^2*((x+3)^(1/4))^2+(5^(1/4))^2)^(1/2)+2^(1/4)*(x+3)^(1/4)))+((-35)*5^(1/4)*log(2^(1/4)*(x+3)^(1/4)+5^(1/4))+(35*5^(1/4)*log(2^(1/4)*(x+3)^(1/4)+(-1)*5^(1/4))+(20...`
- `random-036` [TragerKauers] `3/(2*x^2 + 3*x)*(x + 1)^(2/3)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((-1)*(-1)^(1/3)*log(2*4^(1/3)*((x+1)^(1/3))^2+(-1)*((-1)^(1/3))^2*(4^(1/3))^2*(x+1)^(1/3)+(-2)*(-1)^(1/3))+((-1)*4^(1/3)*log(((x+1)^(1/3))^2+(x+1)^(1/3)+1)+(2*(-1)^(1/3)*log((4^(1/3))^2*(x+1)^(1/3)+2...`
- `random-053` [TragerKauers] `(-3*x + 2)/(2*x - 2)*(3*x)^(-4/5)` -- Mathilda: DECLINE `{"failed", "Kauers' heuristic: no principal power of the residue ideal"}`; FriCAS: ELEM `((2*81^(1/5)*(((-12)*(81^(1/5))^2*rootOf((8*(81^(1/5))^3*rootOf((16*%%J0^4*(81^(1/5))^4+8*%%J0^3*(-1)^(1/5)*(81^(1/5))^3+4*%%J0^2*((-1)^(1/5))^2*(81^(1/5))^2+2*%%J0*((-1)^(1/5))^3*81^(1/5)+((-1)^(1/5)...`

### PossibleZeroQ false negatives (7)

Verified antiderivatives (40-digit check at every point) for which `PossibleZeroQ[D[F, x] - f]` returned False -- a PossibleZeroQ bug, not an integrator bug.

- `binom-032` [Trager] `1/(x^4 + 1)^(1/4)`: `1/4 ArcTan[(-2 x^3 (1 + x^4)^(1/4) + 2 x (1 + x^4)^(3/4))/(1 + 2 x^4 - 2 x^2 Sqrt[1 + x^4])] + 1/4 Log[1 + 2 x^4 + 2 x^3 (1 + x^4)^(1/4) + 2 x^2 Sqrt[1 + x^4] + 2 x (1 + x^4)^(3/4)]`
- `goursat-010` [Trager] `1/((x^4 + 1)^(1/4))`: `1/4 ArcTan[(-2 x^3 (1 + x^4)^(1/4) + 2 x (1 + x^4)^(3/4))/(1 + 2 x^4 - 2 x^2 Sqrt[1 + x^4])] + 1/4 Log[1 + 2 x^4 + 2 x^3 (1 + x^4)^(1/4) + 2 x^2 Sqrt[1 + x^4] + 2 x (1 + x^4)^(3/4)]`
- `goursat-011` [Trager] `x^2/((x^4 - 1)^(3/4))`: `-1/4 ArcTan[(-2 x^3 (-1 + x^4)^(1/4) + 2 x (-1 + x^4)^(3/4))/(-1 + 2 x^4 - 2 x^2 Sqrt[-1 + x^4])] + 1/4 Log[-1 + 2 x^4 + 2 x^3 (-1 + x^4)^(1/4) + 2 x^2 Sqrt[-1 + x^4] + 2 x (-1 + x^4)^(3/4)]`
- `deriv-log-13` [Trager] `(-32*x^8 - 44*x^7 + 32*x^6 + 147*x^5 + 72*x^4 - 45*x^3 + 22*x^2 + 46*x - 2)/(2*(-2*x^6 + 3*x^5 + 3*x^4 - 3*x^3 + x^2 + 2*x - 3)*(8*x^8 - 4*x^7 - 22*x^6 - 3*x^5 + 9*x^4 + 3*x^3 + 20*x^2 + 22*x + 7))*Sqrt[-2*x^6 + 3*x^5 + 3*x^4 - 3*x^3 + x^2 + 2*x - 3]`: `1/2 Log[(-1 - 2 x - 14 x^2 - 21 x^3 + x^4 - 3 x^5 - 22 x^6 - 4 x^7 + 8 x^8 - 4 Sqrt[-3 + 2 x + x^2 - 3 x^3 + 3 x^4 + 3 x^5 - 2 x^6] - 14 x Sqrt[-3 + 2 x + x^2 - 3 x^3 + 3 x^4 + 3 x^5 - 2 x^6] - 16 x^2...`
- `deriv-log-16` [Trager] `(-9*x^5 + 27*x^4 - 43*x^3 + 27*x^2 + 2*x - 8)/(2*(x^3 - 2*x^2 + 2*x - 1)*(9*x^6 - 18*x^5 + 3*x^4 - 13*x^3 + 21*x^2 + 4*x + 10))*Sqrt[x^3 - 2*x^2 + 2*x - 1]`: `1/2 Log[(8 + 8 x + 17 x^2 - 11 x^3 + 3 x^4 - 18 x^5 + 9 x^6 - 6 Sqrt[-1 + 2 x - 2 x^2 + x^3] - 2 x Sqrt[-1 + 2 x - 2 x^2 + x^3] - 6 x^2 Sqrt[-1 + 2 x - 2 x^2 + x^3] + 6 x^3 Sqrt[-1 + 2 x - 2 x^2 + x^3...`
- `deriv-log-19` [Trager] `(-9*x^8 - 3*x^7 - 12*x^6 + 57*x^5 + 51*x^4 - 123*x^3 + 60*x^2 + 18*x - 11)/(2*(x^5 + x^4 + x^3 + 3*x^2 + 3*x - 2)*(-9*x^7 + 6*x^6 - 4*x^5 - 16*x^4 - 16*x^3 + 34*x^2 - 17*x + 3))*Sqrt[x^5 + x^4 + x^3 + 3*x^2 + 3*x - 2]`: `1/2 Log[(-1 + 13 x - 32 x^2 + 4 x^3 + 28 x^4 + 4 x^5 + 12 x^6 + 9 x^7 - 2 Sqrt[-2 + 3 x + 3 x^2 + x^3 + x^4 + x^5] + 8 x Sqrt[-2 + 3 x + 3 x^2 + x^3 + x^4 + x^5] - 6 x^2 Sqrt[-2 + 3 x + 3 x^2 + x^3 + ...`
- `deriv-log-20` [Trager] `(9*x^7 - 13*x^6 - 5*x^5 - 18*x^4 + 3*x^3 + 5*x^2 - 3*x - 2)/(2*(x^5 - 3*x^4 - x - 1)*(-x^7 + 3*x^6 + 9*x^4 + 7*x^3 + 8*x^2 + 2*x + 1))*Sqrt[x^5 - 3*x^4 - x - 1]`: `1/2 Log[(1 + 2 x + 6 x^2 + 5 x^3 + 9 x^4 - 3 x^6 + x^7 + 2 x Sqrt[-1 - x - 3 x^4 + x^5] + 2 x^2 Sqrt[-1 - x - 3 x^4 + x^5] + 6 x^3 Sqrt[-1 - x - 3 x^4 + x^5])/(-1 - 2 x - 8 x^2 - 7 x^3 - 9 x^4 - 3 x^6...`

### FriCAS answers failing their derivative check (1)

- `random-048` `(-3*x^2 + 3*x + 2)/(3*x^2 - x + 2)*(2*x + 1)^(-1/2)`: `((-5152)*39^(1/2)*atan((1990144*23^(1/2)*39^(1/2)+(-11481600)*23^(1/2))/((197*23^(1/2)*39^(1/2)+(-1248)*23^(1/2))*20631^(1/4)*((12608*39^(1/2)+(-79872))/(12608*39^(1/2)+(-78745)))^(1/2)*((((-172031441...` numeric [['3/7', 'bad', '0.53'], ['5/4', 'bad', '0.106'], ['7/3', 'bad', '0.0302']]

## All cases

| id | family | integrand | expect | Trager | s | TragerKauers | s | FriCAS | s | class T | class TK |
|---|---|---|---|---|---:|---|---:|---|---:|---|---|
| binom-001 | binomial-xk | `x^2/Sqrt[x^2 - 3]` | elem | ELEM | 0.026 | ELEM | 0.023 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| binom-002 | binomial-xk | `x^3/(x^2 + 2)^(1/3)` | elem | ELEM | 0.020 | ELEM | 0.019 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-003 | binomial-xk | `x^3/Sqrt[x^6 + 2]` | nonelem | NONELEM | 0.016 | NONELEM | 0.015 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| binom-004 | binomial-xk | `x^2/(x^3 + 1)^(1/3)` | elem | ELEM | 0.015 | ELEM | 0.015 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-005 | binomial-xk | `x^2/(x^4 + 1)^(1/3)` | nonelem | NONELEM | 0.015 | NONELEM | 0.016 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| binom-006 | binomial-xk | `x^3/(x^4 - 1)^(1/3)` | elem | ELEM | 0.021 | ELEM | 0.016 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-007 | binomial-xk | `x^3/(x^6 + 1)^(1/3)` | nonelem | NONELEM | 0.017 | NONELEM | 0.016 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-008 | binomial-xk | `x^3/(x^3 - 1)^(1/4)` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-009 | binomial-xk | `x/(x^4 - 1)^(1/3)` | nonelem | NONELEM | 0.018 | NONELEM | 0.015 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-010 | binomial-xk | `x^2/(x^2 + 2)^(1/3)` | nonelem | NONELEM | 0.017 | NONELEM | 0.017 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-011 | binomial-xk | `x^5/(x^2 - 3)^(1/4)` | elem | ELEM | 0.023 | ELEM | 0.022 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-012 | binomial-xk | `x^5/(x^2 + 2)^(1/4)` | elem | ELEM | 0.024 | ELEM | 0.025 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-013 | binomial-xk | `x^3/(x^4 - 1)^(1/4)` | elem | ELEM | 0.019 | ELEM | 0.017 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-014 | binomial-xk | `x^5/(x^6 - 1)^(1/3)` | elem | ELEM | 0.017 | ELEM | 0.017 | ELEM | 0.050 | AGREE-ELEM | AGREE-ELEM |
| binom-015 | binomial-xk | `x^3/(x^6 + 2)^(1/3)` | nonelem | NONELEM | 0.017 | NONELEM | 0.018 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-016 | binomial-xk | `x^2/Sqrt[x^2 + 1]` | elem | ELEM | 0.022 | ELEM | 0.022 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| binom-017 | binomial-xk | `x/Sqrt[x^2 - 1]` | elem | ELEM | 0.015 | ELEM | 0.014 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-018 | binomial-xk | `x^2/(x^3 - 1)^(1/3)` | elem | ELEM | 0.016 | ELEM | 0.016 | ELEM | 0.050 | AGREE-ELEM | AGREE-ELEM |
| binom-019 | binomial-xk | `1/(x^2 - 1)^(1/3)` | nonelem | NONELEM | 0.018 | NONELEM | 0.016 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-020 | binomial-xk | `x^3/(x^6 - 1)^(1/4)` | nonelem | NONELEM | 0.018 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-021 | binomial-xk | `1/(x^2 + 1)^(1/4)` | nonelem | NONELEM | 0.017 | NONELEM | 0.015 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-022 | binomial-xk | `x^5/Sqrt[x^3 - 1]` | elem | ELEM | 0.020 | ELEM | 0.020 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-023 | binomial-xk | `x/(x^3 + 2)^(1/3)` | nonelem | NONELEM | 0.016 | NONELEM | 0.015 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-024 | binomial-xk | `1/(x^6 - 1)^(1/3)` | nonelem | NONELEM | 0.015 | NONELEM | 0.015 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-025 | binomial-xk | `x^3/Sqrt[x^3 - 3]` | nonelem | NONELEM | 0.020 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-026 | binomial-xk | `1/Sqrt[x^6 - 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.016 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-027 | binomial-xk | `x/(x^4 - 1)^(1/4)` | nonelem | NONELEM | 0.019 | NONELEM | 0.026 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-028 | binomial-xk | `x/Sqrt[x^6 + 2]` | nonelem | NONELEM | 0.015 | NONELEM | 0.019 | NONELEM | 0.070 | AGREE-NONELEM | AGREE-NONELEM |
| binom-029 | binomial-xk | `x^3/(x^3 + 2)^(1/3)` | elem | ELEM | 0.318 | DECLINE | 0.036 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| binom-030 | binomial-xk | `x^2/Sqrt[x^6 + 1]` | elem | ELEM | 0.088 | ELEM | 0.026 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| binom-031 | binomial-xk | `x^2/(x^3 - 1)^(1/4)` | elem | ELEM | 0.018 | ELEM | 0.017 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-032 | binomial-xk | `1/(x^4 + 1)^(1/4)` | elem | ELEM | 0.396 | DECLINE | 0.036 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| binom-033 | binomial-xk | `x^5/(x^4 + 2)^(1/4)` | nonelem | NONELEM | 0.019 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-034 | binomial-xk | `x^3/(x^2 - 1)^(1/4)` | elem | ELEM | 0.020 | ELEM | 0.020 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-035 | binomial-xk | `1/(x^6 + 1)^(1/4)` | nonelem | NONELEM | 0.016 | NONELEM | 0.016 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-036 | binomial-xk | `x^5/(x^2 + 1)^(1/3)` | elem | ELEM | 0.021 | ELEM | 0.021 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-037 | binomial-xk | `x^5/(x^2 - 3)^(1/3)` | elem | ELEM | 0.020 | ELEM | 0.020 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-038 | binomial-xk | `x/(x^4 - 3)^(1/4)` | nonelem | NONELEM | 0.016 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-039 | binomial-xk | `x^3/Sqrt[x^6 - 3]` | nonelem | NONELEM | 0.016 | NONELEM | 0.018 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-040 | binomial-xk | `1/(x^6 - 3)^(1/3)` | nonelem | NONELEM | 0.015 | NONELEM | 0.015 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-041 | binomial-xk | `x^2/(x^2 + 2)^(1/4)` | nonelem | NONELEM | 0.018 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-042 | binomial-xk | `x^2/Sqrt[x^6 - 3]` | elem | ELEM | 0.090 | ELEM | 0.027 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| binom-043 | binomial-xk | `x/Sqrt[x^6 - 1]` | nonelem | NONELEM | 0.014 | NONELEM | 0.017 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-044 | binomial-xk | `1/Sqrt[x^6 - 3]` | nonelem | NONELEM | 0.015 | NONELEM | 0.016 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| binom-045 | binomial-xk | `x^5/(x^6 - 3)^(1/4)` | elem | ELEM | 0.018 | ELEM | 0.020 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-001 | third-kind | `1/((x)*Sqrt[x^2 + 1])` | elem | ELEM | 0.022 | ELEM | 0.023 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-002 | third-kind | `1/((x - 1)*Sqrt[x^2 + 1])` | elem | ELEM | 0.031 | ELEM | 0.029 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-003 | third-kind | `1/((x + 1)*Sqrt[x^2 + 1])` | elem | ELEM | 0.031 | ELEM | 0.029 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-004 | third-kind | `1/((x - 2)*Sqrt[x^2 + 1])` | elem | ELEM | 0.031 | ELEM | 0.027 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-005 | third-kind | `1/((x - 3)*Sqrt[x^2 + 1])` | elem | ELEM | 0.030 | ELEM | 0.029 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-006 | third-kind | `1/((x)*Sqrt[x^2 - 2*x + 3])` | elem | ELEM | 0.030 | ELEM | 0.030 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-007 | third-kind | `1/((x - 1)*Sqrt[x^2 - 2*x + 3])` | elem | ELEM | 0.032 | ELEM | 0.027 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-008 | third-kind | `1/((x + 1)*Sqrt[x^2 - 2*x + 3])` | elem | ELEM | 0.031 | ELEM | 0.029 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-009 | third-kind | `1/((x - 2)*Sqrt[x^2 - 2*x + 3])` | elem | ELEM | 0.031 | ELEM | 0.031 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-010 | third-kind | `1/((x - 3)*Sqrt[x^2 - 2*x + 3])` | elem | ELEM | 0.031 | ELEM | 0.028 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-011 | third-kind | `1/((x)*Sqrt[x^3 + 1])` |  | ELEM | 0.065 | ELEM | 0.029 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-012 | third-kind | `1/((x - 1)*Sqrt[x^3 + 1])` |  | NONELEM | 0.048 | DECLINE | 0.058 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-013 | third-kind | `1/((x - 2)*Sqrt[x^3 + 1])` |  | NONELEM | 0.070 | DECLINE | 0.050 | NONELEM | 0.110 | AGREE-NONELEM | NEITHER |
| third-014 | third-kind | `1/((x - 3)*Sqrt[x^3 + 1])` |  | NONELEM | 0.054 | DECLINE | 0.092 | NONELEM | 0.120 | AGREE-NONELEM | NEITHER |
| third-015 | third-kind | `1/((x - 2)*Sqrt[x^3 - x])` |  | NONELEM | 0.051 | DECLINE | 0.092 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-016 | third-kind | `1/((x - 3)*Sqrt[x^3 - x])` |  | NONELEM | 0.049 | DECLINE | 0.113 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-017 | third-kind | `1/((x)*Sqrt[x^4 + 1])` |  | ELEM | 0.080 | ELEM | 0.023 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-018 | third-kind | `1/((x - 1)*Sqrt[x^4 + 1])` |  | NONELEM | 0.093 | DECLINE | 0.062 | NONELEM | 0.140 | AGREE-NONELEM | NEITHER |
| third-019 | third-kind | `1/((x + 1)*Sqrt[x^4 + 1])` |  | NONELEM | 0.088 | DECLINE | 0.059 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-020 | third-kind | `1/((x - 2)*Sqrt[x^4 + 1])` |  | NONELEM | 0.084 | DECLINE | 0.185 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-021 | third-kind | `1/((x - 3)*Sqrt[x^4 + 1])` |  | NONELEM | 0.072 | DECLINE | 0.212 | NONELEM | 0.140 | AGREE-NONELEM | NEITHER |
| third-022 | third-kind | `1/((x)*Sqrt[x^4 - 1])` |  | ELEM | 0.066 | ELEM | 0.031 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| third-023 | third-kind | `1/((x - 2)*Sqrt[x^4 - 1])` |  | NONELEM | 0.049 | DECLINE | 0.159 | NONELEM | 0.150 | AGREE-NONELEM | NEITHER |
| third-024 | third-kind | `1/((x - 3)*Sqrt[x^4 - 1])` |  | NONELEM | 0.049 | DECLINE | 0.203 | NONELEM | 0.150 | AGREE-NONELEM | NEITHER |
| third-025 | third-kind | `1/((x)*Sqrt[x^4 + x^2 + 1])` |  | ELEM | 0.067 | ELEM | 0.026 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| third-026 | third-kind | `1/((x - 1)*Sqrt[x^4 + x^2 + 1])` |  | NONELEM | 0.100 | DECLINE | 0.065 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-027 | third-kind | `1/((x + 1)*Sqrt[x^4 + x^2 + 1])` |  | NONELEM | 0.096 | DECLINE | 0.063 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-028 | third-kind | `1/((x - 2)*Sqrt[x^4 + x^2 + 1])` |  | NONELEM | 0.064 | DECLINE | 0.222 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-029 | third-kind | `1/((x - 3)*Sqrt[x^4 + x^2 + 1])` |  | NONELEM | 0.076 | DECLINE | 0.287 | NONELEM | 0.130 | AGREE-NONELEM | NEITHER |
| third-030 | third-kind | `1/((x)*Sqrt[x^3 + x + 1])` |  | NONELEM | 0.038 | DECLINE | 0.040 | NONELEM | 0.100 | AGREE-NONELEM | NEITHER |
| ratrad-001 | rational-times-radical | `x^3/x*(x^3 + 1)^(-2/3)` |  | ELEM | 0.015 | ELEM | 0.016 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| ratrad-002 | rational-times-radical | `x^3/(x - 2)*(x^2 + 1)^(1/6)` |  | NONELEM | 0.265 | DECLINE | 5.422 | TIMEOUT | 60.068 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-003 | rational-times-radical | `(x^2 + 1)/(x + 1)*(x^2 - x + 1)^(-1/2)` | elem | ELEM | 0.040 | ELEM | 0.041 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-004 | rational-times-radical | `x^3/(x - 2)*(x^4 + 3)^(-3/5)` |  | NONELEM | 0.031 | NONELEM | 0.031 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-005 | rational-times-radical | `x^3/(x - 2)*(x^3 + 1)^(3/2)` |  | NONELEM | 0.024 | NONELEM | 0.024 | NONELEM | 0.120 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-006 | rational-times-radical | `1/x*(x^3 + 1)^(2/3)` |  | ELEM | 0.343 | DECLINE | 0.088 | ELEM | 0.070 | AGREE-ELEM | FRICAS-ONLY |
| ratrad-007 | rational-times-radical | `1/(x + 1)*(x^3 - 2*x)^(1/4)` |  | NONELEM | 0.236 | DECLINE | 0.072 | TIMEOUT | 60.059 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-008 | rational-times-radical | `(2*x - 1)/(x + 1)*(x^2 + 1)^(-2/3)` |  | NONELEM | 0.056 | DECLINE | 0.057 | TIMEOUT | 60.040 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-009 | rational-times-radical | `x^3*(x^2 + 1)^(1/5)` |  | ELEM | 0.022 | ELEM | 0.021 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| ratrad-010 | rational-times-radical | `(2*x - 1)/(x + 1)*(x^2 - x + 1)^(1/5)` |  | NONELEM | 0.218 | DECLINE | 1.506 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-011 | rational-times-radical | `x^3/x*(x^2 - x + 1)^(1/5)` |  | NONELEM | 0.022 | NONELEM | 0.021 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-012 | rational-times-radical | `x/(x - 2)*(x^3 + 1)^(-2/3)` |  | NONELEM | 0.062 | DECLINE | 0.066 | TIMEOUT | 60.035 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-013 | rational-times-radical | `1/x*(x^2 - x + 1)^(3/4)` |  | NONELEM | 0.019 | NONELEM | 0.021 | TIMEOUT | 60.050 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-014 | rational-times-radical | `x^3/(x + 1)*Sqrt[x^3 - 2*x]` |  | NONELEM | 0.025 | NONELEM | 0.024 | NONELEM | 0.140 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-015 | rational-times-radical | `(x^2 + 1)/(x - 2)*(x^3 - 2*x)^(-3/5)` |  | NONELEM | 0.034 | NONELEM | 0.034 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-016 | rational-times-radical | `(x^2 + 1)/x*(x^4 + 3)^(1/6)` |  | NONELEM | 7.256 | DECLINE | 0.121 | TIMEOUT | 60.063 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-017 | rational-times-radical | `x*(x^4 + 3)^(3/2)` |  | ELEM | 0.094 | ELEM | 0.035 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| ratrad-018 | rational-times-radical | `1/x*(x^3 - 2*x)^(3/2)` |  | NONELEM | 0.023 | NONELEM | 0.024 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-019 | rational-times-radical | `1/(x^2 + 2)*(x^3 + 1)^(-1/4)` |  | NONELEM | 0.903 | DECLINE | 5.956 | TIMEOUT | 60.039 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-020 | rational-times-radical | `x^3/(x - 2)*(x^3 - 2*x)^(-1/2)` |  | NONELEM | 0.021 | NONELEM | 0.021 | NONELEM | 0.160 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-021 | rational-times-radical | `x^3/(x - 2)*(x^2 - x + 1)^(-1/6)` |  | NONELEM | 0.040 | NONELEM | 0.040 | TIMEOUT | 60.060 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-022 | rational-times-radical | `x/x*(x + 2)^(-1/2)` | elem | ELEM | 0.019 | ELEM | 0.020 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-023 | rational-times-radical | `(x^2 + 1)*(x^2 + 1)^(1/5)` |  | NONELEM | 0.025 | NONELEM | 0.024 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-024 | rational-times-radical | `(2*x - 1)*(x^2 - x + 1)^(1/5)` |  | ELEM | 0.031 | ELEM | 0.033 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-025 | rational-times-radical | `(2*x - 1)/(x + 1)*(x + 2)^(1/3)` | elem | ELEM | 0.198 | DECLINE | 0.060 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| ratrad-026 | rational-times-radical | `x^3*(x + 2)^(1/4)` | elem | ELEM | 0.037 | ELEM | 0.035 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| ratrad-027 | rational-times-radical | `x^3*(x + 2)^(-1/2)` | elem | ELEM | 0.027 | ELEM | 0.028 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-028 | rational-times-radical | `x^3*(x^2 + 1)^(1/6)` |  | ELEM | 0.033 | ELEM | 0.033 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-029 | rational-times-radical | `(x^2 + 1)/(x^2 + 2)*(x^2 + 1)^(3/4)` |  | NONELEM | 0.027 | NONELEM | 0.025 | TIMEOUT | 60.047 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-030 | rational-times-radical | `x*(x + 2)^(-1/4)` | elem | ELEM | 0.028 | ELEM | 0.028 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| ratrad-031 | rational-times-radical | `x/(x + 1)*(x^2 + 1)^(1/5)` |  | NONELEM | 0.206 | DECLINE | 1.688 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-032 | rational-times-radical | `(2*x - 1)/(x + 1)*(x^2 + 1)^(2/5)` |  | NONELEM | 0.265 | DECLINE | 1.500 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-033 | rational-times-radical | `x/(x - 2)*(x + 2)^(-1/3)` | elem | ELEM | 0.047 | DECLINE | 0.075 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| ratrad-034 | rational-times-radical | `x^3/x*(x^2 + 1)^(5/6)` |  | NONELEM | 0.025 | NONELEM | 0.030 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-035 | rational-times-radical | `x/x*(x^4 + 3)^(5/6)` |  | NONELEM | 0.027 | NONELEM | 0.033 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-036 | rational-times-radical | `(2*x - 1)/(x^2 + 2)*(x^4 + 3)^(5/6)` |  | NONELEM | 0.040 | NONELEM | 0.039 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-037 | rational-times-radical | `1/(x + 1)*(x^3 - 2*x)^(-1/4)` |  | NONELEM | 0.295 | DECLINE | 1.080 | TIMEOUT | 60.034 | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-038 | rational-times-radical | `1/(x^2 + 2)*(x^3 + 1)^(2/5)` |  | NONELEM | 0.025 | NONELEM | 0.026 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-039 | rational-times-radical | `(2*x - 1)/x*(x^3 - 2*x)^(1/6)` |  | NONELEM | 0.028 | NONELEM | 0.029 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-040 | rational-times-radical | `(2*x - 1)/(x + 1)*(x^4 + 3)^(2/5)` |  | NONELEM | 0.030 | NONELEM | 0.030 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| ratrad-041 | rational-times-radical | `x^3*Sqrt[x + 2]` | elem | ELEM | 0.028 | ELEM | 0.030 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-042 | rational-times-radical | `(x^2 + 1)*(x^4 + 3)^(-1/3)` |  | NONELEM | 0.021 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratrad-043 | rational-times-radical | `x^3*(x^3 - 2*x)^(-2/3)` |  | ELEM | 0.354 | DECLINE | 0.060 | TIMEOUT | 60.032 | MATHILDA-ONLY | NEITHER |
| ratrad-044 | rational-times-radical | `x/(x - 2)*(x^3 - 2*x)^(-1/3)` |  | NONELEM | 0.493 | DECLINE | 0.160 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-045 | rational-times-radical | `1/x*(x^2 + 1)^(-1/6)` |  | ELEM | 3.515 | DECLINE | 0.062 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| ratrad-046 | rational-times-radical | `x^3/x*(x^3 + 1)^(1/3)` |  | ELEM | 0.027 | ELEM | 0.027 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratrad-047 | rational-times-radical | `(2*x - 1)/x*(x + 2)^(-1/4)` | elem | ELEM | 0.061 | DECLINE | 2.392 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| ratrad-048 | rational-times-radical | `(2*x - 1)/(x^2 + 2)*(x^2 - x + 1)^(1/3)` |  | NONELEM | 0.111 | DECLINE | 0.922 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-049 | rational-times-radical | `(2*x - 1)/(x - 2)*(x^2 + 1)^(-3/5)` |  | NONELEM | 0.330 | DECLINE | 1.072 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| ratrad-050 | rational-times-radical | `x/(x + 1)*(x^2 - x + 1)^(1/3)` |  | NONELEM | 0.607 | DECLINE | 0.034 | TIMEOUT | 60.036 | NONELEM-UNCONFIRMED | NEITHER |
| cheb-001 | chebyshev | `x*(x^2 + 1)^(1/3)` | elem | ELEM | 0.023 | ELEM | 0.023 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-002 | chebyshev | `x^2*Sqrt[2*x^4 - 1]` | nonelem | NONELEM | 0.019 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-003 | chebyshev | `x^4*(-x^3 + 3)^(-1/2)` | nonelem | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-004 | chebyshev | `x*(2*x^3 - 1)^(-3/4)` | nonelem | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-005 | chebyshev | `(2*x^4 - 1)^(1/4)` | nonelem | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-006 | chebyshev | `(-x^4 + 3)^(1/6)` | nonelem | NONELEM | 0.022 | NONELEM | 0.023 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-007 | chebyshev | `x^(-1)*Sqrt[-x^3 + 3]` | elem | ELEM | 0.124 | DECLINE | 0.081 | ELEM | 0.070 | AGREE-ELEM | FRICAS-ONLY |
| cheb-008 | chebyshev | `x^4*(-x^4 + 3)^(2/3)` | nonelem | NONELEM | 0.021 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-009 | chebyshev | `x^3*(-x^4 + 3)^(-1/2)` | elem | ELEM | 0.020 | ELEM | 0.019 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-010 | chebyshev | `x^5*(x^3 + 1)^(1/4)` | elem | ELEM | 0.028 | ELEM | 0.028 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-011 | chebyshev | `x^(-2)*(x^4 + 1)^(-3/4)` | elem | ELEM | 0.023 | ELEM | 0.023 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| cheb-012 | chebyshev | `(-x^2 + 3)^(1/4)` | nonelem | NONELEM | 0.019 | NONELEM | 0.019 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-013 | chebyshev | `x^(-1)*Sqrt[x^4 + 1]` | elem | ELEM | 0.118 | DECLINE | 0.061 | ELEM | 0.070 | AGREE-ELEM | FRICAS-ONLY |
| cheb-014 | chebyshev | `x^3*(2*x^3 - 1)^(2/3)` | elem | ELEM | 0.609 | DECLINE | 0.051 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| cheb-015 | chebyshev | `x^3*(-x^2 + 3)^(1/3)` | elem | ELEM | 0.030 | ELEM | 0.028 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| cheb-016 | chebyshev | `x^5*(x^2 + 1)^(1/6)` | elem | ELEM | 0.031 | ELEM | 0.030 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-017 | chebyshev | `x^5*(2*x^3 - 1)^(-3/4)` | elem | ELEM | 0.028 | ELEM | 0.027 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-018 | chebyshev | `x^(-2)*(x^3 + 1)^(2/5)` | nonelem | NONELEM | 0.025 | NONELEM | 0.025 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-019 | chebyshev | `(3*x^2 + 2)^(-1/4)` | nonelem | NONELEM | 0.021 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-020 | chebyshev | `x^(-1)*(2*x^4 - 1)^(2/3)` | elem | ELEM | 0.516 | DECLINE | 0.555 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| cheb-021 | chebyshev | `x*(3*x^4 + 2)^(1/3)` | nonelem | NONELEM | 0.019 | NONELEM | 0.019 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-022 | chebyshev | `x*(2*x^3 - 1)^(3/2)` | nonelem | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-023 | chebyshev | `(3*x^4 + 2)^(1/6)` | nonelem | NONELEM | 0.022 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-024 | chebyshev | `x^4*(x^3 + 1)^(-3/4)` | nonelem | NONELEM | 0.022 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-025 | chebyshev | `x^3*(3*x^4 + 2)^(-5/6)` | elem | ELEM | 0.031 | ELEM | 0.028 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| cheb-026 | chebyshev | `x^(-3)*(-x^3 + 3)^(3/2)` | nonelem | NONELEM | 0.023 | NONELEM | 0.023 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-027 | chebyshev | `x^(-1)*(-x^2 + 3)^(-1/3)` | elem | ELEM | 0.594 | DECLINE | 0.051 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| cheb-028 | chebyshev | `x*(-x^3 + 3)^(-5/6)` | nonelem | NONELEM | 0.025 | NONELEM | 0.025 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-029 | chebyshev | `x^3*(2*x^2 - 1)^(3/2)` | elem | ELEM | 0.027 | ELEM | 0.026 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| cheb-030 | chebyshev | `x^4*(2*x^2 - 1)^(1/4)` | nonelem | NONELEM | 0.021 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-031 | chebyshev | `x^(-1)*Sqrt[-x^2 + 3]` | elem | ELEM | 0.038 | ELEM | 0.048 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| cheb-032 | chebyshev | `x^3*(-x^2 + 3)^(-5/6)` | elem | ELEM | 0.035 | ELEM | 0.035 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-033 | chebyshev | `x^4*(3*x^4 + 2)^(3/2)` | nonelem | NONELEM | 0.023 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-034 | chebyshev | `x^5*(3*x^2 + 2)^(3/2)` | elem | ELEM | 0.027 | ELEM | 0.028 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| cheb-035 | chebyshev | `x^(-3)*(3*x^3 + 2)^(1/4)` | nonelem | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-036 | chebyshev | `x^5*(3*x^3 + 2)^(-1/2)` | elem | ELEM | 0.024 | ELEM | 0.023 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| cheb-037 | chebyshev | `x^5*(2*x^3 - 1)^(-1/3)` | elem | ELEM | 0.023 | ELEM | 0.023 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| cheb-038 | chebyshev | `x^(-3)*(-x^4 + 3)^(3/2)` | elem | ELEM | 0.206 | DECLINE | 0.085 | ELEM | 0.080 | AGREE-ELEM | FRICAS-ONLY |
| cheb-039 | chebyshev | `x^4*(3*x^4 + 2)^(2/5)` | nonelem | NONELEM | 0.023 | NONELEM | 0.023 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| cheb-040 | chebyshev | `x^4*(2*x^4 - 1)^(-1/3)` | nonelem | NONELEM | 0.021 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| goursat-001 | goursat | `(x^2 - 1)/((x^2 + 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.121 | ELEM | 0.037 | ELEM | 0.130 | AGREE-ELEM | AGREE-ELEM |
| goursat-002 | goursat | `(x^2 + 1)/((x^2 - 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.118 | ELEM | 0.034 | ELEM | 0.140 | AGREE-ELEM | AGREE-ELEM |
| goursat-003 | goursat | `x/Sqrt[x^4 + 10*x^2 - 96*x - 71]` | elem | ELEM | 0.155 | DECLINE | 0.222 | ELEM | 0.190 | AGREE-ELEM | FRICAS-ONLY |
| goursat-004 | goursat | `x/Sqrt[x^4 + 10*x^2 - 96*x - 72]` |  | NONELEM | 0.054 | DECLINE | 0.213 | NONELEM | 0.170 | AGREE-NONELEM | NEITHER |
| goursat-005 | goursat | `(1 - x^2)/((1 + x^2)*Sqrt[1 + x^2 + x^4])` | elem | ELEM | 0.124 | ELEM | 0.037 | ELEM | 0.140 | AGREE-ELEM | AGREE-ELEM |
| goursat-006 | goursat | `1/(x^3 + 1)^(1/3)` | elem | ELEM | 0.297 | DECLINE | 0.043 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| goursat-007 | goursat | `x/((x^3 + 1)^(1/3)*(x^3 - 2))` |  | NONELEM | 0.067 | DECLINE | 2.201 | NONELEM | 2.830 | AGREE-NONELEM | NEITHER |
| goursat-008 | goursat | `x/((x^3 + 8)*(x^3 + 1)^(1/3))` |  | NONELEM | 0.080 | DECLINE | 3.188 | NONELEM | 3.510 | AGREE-NONELEM | NEITHER |
| goursat-009 | goursat | `(x^3 - 2)/((x^3 + 1)^(2/3)*(x^3 + 4))` |  | NONELEM | 0.077 | DECLINE | 0.295 | NONELEM | 3.120 | AGREE-NONELEM | NEITHER |
| goursat-010 | goursat | `1/((x^4 + 1)^(1/4))` | elem | ELEM | 0.496 | DECLINE | 0.043 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| goursat-011 | goursat | `x^2/((x^4 - 1)^(3/4))` |  | ELEM | 0.390 | ELEM | 0.051 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| goursat-012 | goursat | `1/(x*(x^4 + 1)^(1/4))` | elem | ELEM | 0.440 | DECLINE | 0.043 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| goursat-013 | goursat | `(x^4 - 1)/((x^4 + 1)^(1/4)*x^2)` |  | NONELEM | 0.029 | NONELEM | 0.029 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| goursat-014 | goursat | `1/((x^2 + 1)*(x^4 + 1)^(1/4))` |  | NONELEM | 0.320 | DECLINE | 1.718 | TIMEOUT | 60.047 | NONELEM-UNCONFIRMED | NEITHER |
| goursat-015 | goursat | `(x^2 - 1)/((x^2 + 1)*(x^4 + 1)^(1/4))` |  | NONELEM | 1.524 | DECLINE | 0.047 | TIMEOUT | 60.046 | NONELEM-UNCONFIRMED | NEITHER |
| goursat-016 | goursat | `x/((x^2 + 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.115 | ELEM | 0.034 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| goursat-017 | goursat | `(3*x^2 + 1)/((x^2 - 1)*Sqrt[x^4 - x^2])` |  | ELEM | 0.041 | ELEM | 0.046 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| goursat-018 | goursat | `(x^6 + 2)/((x^6 - 1)*Sqrt[x^6 + 1])` |  | NONELEM | 0.129 | DECLINE | 2.719 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| goursat-019 | goursat | `(5*x^4 + 1)/Sqrt[x^5 + x + 1]` | elem | ELEM | 0.020 | ELEM | 0.021 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| goursat-020 | goursat | `(x^2 + 2)/((x^2 + 1)*(x^3 + x)^(1/3))` |  | ELEM | 0.337 | DECLINE | 0.218 | TIMEOUT | 60.035 | MATHILDA-ONLY | NEITHER |
| goursat-021 | goursat | `1/((x^2 - 2)*(x^4 - 4*x^2 + 1)^(1/4))` |  | NONELEM | 0.428 | TIMEOUT | 60.033 | TIMEOUT | 60.045 | NONELEM-UNCONFIRMED | TIMEOUT |
| goursat-022 | goursat | `(x - 1)/((x + 1)*Sqrt[x^3 + x^2 + x])` |  | ELEM | 0.104 | DECLINE | 0.066 | ELEM | 0.120 | AGREE-ELEM | FRICAS-ONLY |
| goursat-023 | goursat | `(x^2 - 2)/((x^2 + 2)*Sqrt[x^4 + 4])` | elem | ELEM | 0.099 | ELEM | 0.034 | ELEM | 0.140 | AGREE-ELEM | AGREE-ELEM |
| recip-001 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.121 | ELEM | 0.036 | ELEM | 0.140 | AGREE-ELEM | AGREE-ELEM |
| recip-002 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + x + 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.137 | DECLINE | 0.153 | ELEM | 0.140 | AGREE-ELEM | FRICAS-ONLY |
| recip-003 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + 2*x + 1)*Sqrt[x^4 + x^2 + 1])` | elem | ELEM | 0.127 | DECLINE | 0.082 | ELEM | 0.140 | AGREE-ELEM | FRICAS-ONLY |
| recip-004 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + 3*x + 1)*Sqrt[x^4 - x^2 + 1])` | elem | ELEM | 0.201 | DECLINE | 1.664 | ELEM | 0.160 | AGREE-ELEM | FRICAS-ONLY |
| recip-005 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 - x + 1)*Sqrt[x^4 + 3*x^2 + 1])` | elem | ELEM | 0.231 | DECLINE | 2.040 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| recip-006 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + 1)*Sqrt[x^4 + 5*x^2 + 1])` | elem | ELEM | 0.186 | ELEM | 0.037 | ELEM | 0.140 | AGREE-ELEM | AGREE-ELEM |
| recip-007 | pseudo-elliptic-recip | `(x^2 - 1)/((x^2 + x + 1)*Sqrt[x^4 + 3*x^2 + 1])` | elem | ELEM | 0.227 | DECLINE | 1.471 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| recip-008 | pseudo-elliptic-recip | `(x^2 + 1)/((x^2 + x - 1)*Sqrt[x^4 + 1])` | elem | ELEM | 0.199 | DECLINE | 1.248 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| recip-009 | pseudo-elliptic-recip | `(x^2 + 1)/((x^2 + 2*x - 1)*Sqrt[x^4 + x^2 + 1])` | elem | ELEM | 0.136 | DECLINE | 2.561 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| recip-010 | pseudo-elliptic-recip | `(x^2 + 1)/((x^2 + 3*x - 1)*Sqrt[x^4 + 3*x^2 + 1])` | elem | ELEM | 0.197 | TIMEOUT | 60.800 | ELEM | 0.170 | AGREE-ELEM | FRICAS-ONLY |
| recip-011 | pseudo-elliptic-recip | `(x^2 + 1)/((x^2 - x - 1)*Sqrt[x^4 - x^2 + 1])` | elem | ELEM | 0.214 | DECLINE | 2.064 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| evensub-001 | even-substitution | `x/((x^2 + 3)*Sqrt[x^4 + x^2 + 1])` | elem | ELEM | 0.121 | ELEM | 0.036 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| evensub-002 | even-substitution | `x^3/((x^2 - 2)*Sqrt[x^4 + 1])` | elem | ELEM | 0.395 | ELEM | 0.049 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| evensub-003 | even-substitution | `x/Sqrt[x^4 - 3*x^2 + 1]` | elem | ELEM | 0.090 | ELEM | 0.036 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| evensub-004 | even-substitution | `x^5/((x^4 + 1)*Sqrt[x^4 + 2])` | elem | ELEM | 0.326 | ELEM | 0.045 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| evensub-005 | even-substitution | `x/((x^4 + x^2 + 1)*Sqrt[x^4 + 1])` | elem | ELEM | 1.010 | ELEM | 0.050 | ELEM | 0.110 | AGREE-ELEM | AGREE-ELEM |
| evensub-006 | even-substitution | `(x^3 + x)/((x^2 + 5)*Sqrt[x^4 - 2*x^2 + 7])` | elem | ELEM | 0.264 | ELEM | 0.062 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| hyper-001 | genus-hyper | `1/Sqrt[x + 3]` | elem | ELEM | 0.018 | ELEM | 0.018 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| hyper-002 | genus-hyper | `x^1/Sqrt[x + 3]` | elem | ELEM | 0.022 | ELEM | 0.022 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| hyper-003 | genus-hyper | `1/Sqrt[x^2 + x + 1]` | elem | ELEM | 0.027 | ELEM | 0.028 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| hyper-004 | genus-hyper | `x^1/Sqrt[x^2 + x + 1]` | elem | ELEM | 0.028 | ELEM | 0.031 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| hyper-005 | genus-hyper | `1/Sqrt[x^3 - x + 1]` | nonelem | NONELEM | 0.017 | NONELEM | 0.017 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-006 | genus-hyper | `x^2/Sqrt[x^3 - x + 1]` |  | NONELEM | 0.019 | NONELEM | 0.019 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-007 | genus-hyper | `1/Sqrt[x^4 + x + 1]` | nonelem | NONELEM | 0.017 | NONELEM | 0.018 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-008 | genus-hyper | `x^2/Sqrt[x^4 + x + 1]` |  | NONELEM | 0.019 | NONELEM | 0.019 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-009 | genus-hyper | `1/Sqrt[x^5 - x + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-010 | genus-hyper | `x^1/Sqrt[x^5 - x + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-011 | genus-hyper | `x^3/Sqrt[x^5 - x + 1]` |  | NONELEM | 0.019 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-012 | genus-hyper | `1/Sqrt[x^6 + x^2 + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-013 | genus-hyper | `x^1/Sqrt[x^6 + x^2 + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-014 | genus-hyper | `x^3/Sqrt[x^6 + x^2 + 1]` |  | NONELEM | 0.019 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-015 | genus-hyper | `1/Sqrt[x^7 + x + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-016 | genus-hyper | `x^2/Sqrt[x^7 + x + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-017 | genus-hyper | `x^4/Sqrt[x^7 + x + 1]` |  | NONELEM | 0.020 | NONELEM | 0.019 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-018 | genus-hyper | `1/Sqrt[x^8 + x^3 + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-019 | genus-hyper | `x^2/Sqrt[x^8 + x^3 + 1]` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| hyper-020 | genus-hyper | `x^4/Sqrt[x^8 + x^3 + 1]` |  | NONELEM | 0.020 | NONELEM | 0.019 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-001 | genus-super | `1/(x^4 + 1)^(1/3)` | nonelem | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-002 | genus-super | `x/(x^4 + 1)^(2/3)` |  | NONELEM | 0.018 | NONELEM | 0.018 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-003 | genus-super | `1/(x^5 + 2)^(1/3)` |  | NONELEM | 0.019 | NONELEM | 0.026 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| super-004 | genus-super | `x/(x^5 + 2)^(2/3)` |  | NONELEM | 0.022 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-005 | genus-super | `1/(x^3 + x + 1)^(1/3)` |  | NONELEM | 0.119 | DECLINE | 0.487 | TIMEOUT | 60.046 | NONELEM-UNCONFIRMED | NEITHER |
| super-006 | genus-super | `x/(x^3 + x + 1)^(2/3)` |  | NONELEM | 0.117 | DECLINE | 0.061 | TIMEOUT | 60.031 | NONELEM-UNCONFIRMED | NEITHER |
| super-007 | genus-super | `1/(x^4 + x + 1)^(1/5)` |  | NONELEM | 0.025 | NONELEM | 0.025 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-008 | genus-super | `x/(x^4 + x + 1)^(2/5)` |  | NONELEM | 0.025 | NONELEM | 0.026 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| super-009 | genus-super | `1/(x^2 + 1)^(1/3)` |  | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| super-010 | genus-super | `x/(x^2 + 1)^(2/3)` |  | ELEM | 0.020 | ELEM | 0.020 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-001 | rational-radicand | `Sqrt[(x + 1)/(x - 1)]` | elem | ELEM | 0.029 | ELEM | 0.032 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-002 | rational-radicand | `((x + 1)/(x - 1))^(1/3)` | elem | ELEM | 0.217 | DECLINE | 0.059 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-003 | rational-radicand | `Sqrt[x/(x + 1)]/x` | elem | ELEM | 0.030 | ELEM | 0.032 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-004 | rational-radicand | `((1 - x)/(1 + x))^(1/4)` | elem | ELEM | 0.093 | DECLINE | 0.085 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-005 | rational-radicand | `x*Sqrt[(x - 1)/(x + 2)]` | elem | ELEM | 0.037 | ELEM | 0.036 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-006 | rational-radicand | `Sqrt[(x^2 + 1)/(x^2 - 1)]` |  | NONELEM | 0.019 | NONELEM | 0.020 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-007 | rational-radicand | `((x + 2)/(x - 3))^(2/3)/x` | elem | ELEM | 0.148 | DECLINE | 0.077 | ELEM | 0.090 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-008 | rational-radicand | `1/(x^2*((x + 1)/(x - 1))^(1/5))` | elem | ELEM | 0.220 | DECLINE | 0.178 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-009 | rational-radicand | `((2*x + 1)/(x + 3))^(1/6)` | elem | ELEM | 0.185 | DECLINE | 7.909 | ELEM | 0.110 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-010 | rational-radicand | `Sqrt[(x - 1)/(x^2 + 1)]` |  | NONELEM | 0.020 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-011 | rational-radicand | `((x^2 + 1)/x)^(1/3)` |  | NONELEM | 0.021 | NONELEM | 0.021 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-012 | rational-radicand | `Sqrt[(x^3 + 1)/x]` |  | NONELEM | 0.020 | NONELEM | 0.020 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-013 | rational-radicand | `(x + 1)/(x^2*Sqrt[(x - 1)/(x + 1)])` | elem | ELEM | 0.065 | ELEM | 0.076 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-014 | rational-radicand | `((x - 2)/(x + 2))^(3/4)/(x - 1)` | elem | ELEM | 0.223 | DECLINE | 0.075 | ELEM | 0.120 | AGREE-ELEM | FRICAS-ONLY |
| ratradicand-015 | rational-radicand | `Sqrt[x/(x^2 + x + 1)]` |  | NONELEM | 0.021 | NONELEM | 0.022 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-016 | rational-radicand | `1/(x - 2)*Sqrt[(x + 1)/(x - 1)]` | elem | ELEM | 0.051 | ELEM | 0.065 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-017 | rational-radicand | `(x/(x^3 + 1))^(1/3)` |  | NONELEM | 0.023 | NONELEM | 0.022 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-018 | rational-radicand | `Sqrt[(2*x - 1)/(3*x + 5)]/(x + 1)` | elem | ELEM | 0.108 | ELEM | 0.089 | ELEM | 0.090 | AGREE-ELEM | AGREE-ELEM |
| ratradicand-019 | rational-radicand | `((x^2 - 1)/(x^2 + 1))^(1/4)` |  | NONELEM | 0.023 | NONELEM | 0.024 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| ratradicand-020 | rational-radicand | `Sqrt[(x^2 + 1)/x]/x` |  | NONELEM | 0.022 | NONELEM | 0.021 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| torsion-001 | torsion-heavy | `1/((x - 3)*(x^3 + 1)^(1/3))` |  | NONELEM | 0.084 | DECLINE | 2.273 | TIMEOUT | 60.049 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-002 | torsion-heavy | `x/((x - 3)*(x^3 + 1)^(2/3))` |  | NONELEM | 0.083 | DECLINE | 0.112 | TIMEOUT | 60.035 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-003 | torsion-heavy | `1/((x + 2)*(x^3 + 2)^(1/3))` |  | NONELEM | 0.074 | DECLINE | 1.088 | TIMEOUT | 60.049 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-004 | torsion-heavy | `1/((2*x - 1)*(x^3 - 1)^(1/3))` |  | NONELEM | 0.076 | DECLINE | 1.616 | TIMEOUT | 60.051 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-005 | torsion-heavy | `x/((x - 2)*(x^3 + 2)^(2/3))` |  | NONELEM | 0.084 | DECLINE | 0.076 | TIMEOUT | 60.029 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-006 | torsion-heavy | `1/((x^2 + 1)*(x^3 + 1)^(1/4))` |  | NONELEM | 0.266 | DECLINE | 6.743 | TIMEOUT | 60.051 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-007 | torsion-heavy | `1/((x^2 + 3)*(x^3 + 1)^(1/4))` |  | NONELEM | 0.211 | DECLINE | 9.794 | TIMEOUT | 60.048 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-008 | torsion-heavy | `x/((x^2 + 2)*(x^3 + 1)^(3/4))` |  | NONELEM | 0.855 | DECLINE | 0.541 | TIMEOUT | 60.037 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-009 | torsion-heavy | `(x^2 + 1)/x*(x^4 + 2)^(1/6)` |  | NONELEM | 9.780 | DECLINE | 0.121 | TIMEOUT | 60.052 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-010 | torsion-heavy | `(x^2 - 1)/x*(x^4 + 5)^(1/6)` |  | NONELEM | 10.348 | DECLINE | 0.918 | TIMEOUT | 60.055 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-011 | torsion-heavy | `1/((x^2 - 3)*(x^4 - 6*x^2 + 1)^(1/4))` |  | NONELEM | 0.362 | TIMEOUT | 60.032 | TIMEOUT | 60.068 | NONELEM-UNCONFIRMED | TIMEOUT |
| torsion-012 | torsion-heavy | `1/((x^2 + 2)*(x^4 + 4)^(1/4))` |  | NONELEM | 0.343 | DECLINE | 2.485 | TIMEOUT | 60.043 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-013 | torsion-heavy | `(2*x + 1)/(-3*x^2 - 2*x + 1)*(x^2 - x + 3)^(1/4)` |  | NONELEM | 0.148 | DECLINE | 0.445 | TIMEOUT | 60.054 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-014 | torsion-heavy | `1/(2*x^2 + 3*x + 1)*(x^3 - x^2 + x + 1)^(-2/3)` |  | NONELEM | 0.096 | DECLINE | 0.095 | TIMEOUT | 60.053 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-015 | torsion-heavy | `x/(x + 2)*(x^2 + 3)^(1/5)` |  | NONELEM | 0.472 | DECLINE | 2.336 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| torsion-016 | torsion-heavy | `(x - 1)/(x + 1)*(x^2 + 2)^(2/5)` |  | NONELEM | 0.246 | DECLINE | 2.328 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| torsion-017 | torsion-heavy | `(-x^2 - x - 1)/(2*x + 3)*(x^4 + x^3 + 2*x^2 - x - 1)^(-5/6)` |  | NONELEM | 5.403 | DECLINE | 3.876 | TIMEOUT | 60.056 | NONELEM-UNCONFIRMED | NEITHER |
| torsion-018 | torsion-heavy | `1/(x - 1)*(x^3 - 3*x)^(1/4)` |  | NONELEM | 0.187 | DECLINE | 0.090 | TIMEOUT | 60.056 | NONELEM-UNCONFIRMED | NEITHER |
| deriv-alg-01 | deriv-alg | `(-68*x^4 - 28*x^3 - 36*x^2 - 28*x)/(4*(-2*x^3 - x^2 - 1))*(-2*x^3 - x^2 - 1)^(3/4)` | elem | ELEM | 0.042 | ELEM | 0.032 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-02 | deriv-alg | `(-4*x^5 - 10*x^4 - 32*x^3 + 27*x^2 + 36*x + 6)/(3*(-2*x^4 + 2*x^3 + 3*x^2 - 2))*(-2*x^4 + ...` | elem | ELEM | 0.064 | ELEM | 0.062 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-03 | deriv-alg | `(10*x^3 + 10*x^2 - 22*x - 22)/(4*(x^3 - 2*x^2 - x + 2))*(x^3 - 2*x^2 - x + 2)^(-3/4)` | elem | ELEM | 0.047 | ELEM | 0.047 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-04 | deriv-alg | `(4*x^2 + 16*x + 12)/(4*(x^2 - 2*x - 3))*(x^2 - 2*x - 3)^(-3/4)` | elem | ELEM | 0.034 | ELEM | 0.034 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-05 | deriv-alg | `(46*x^6 + 97*x^5 + 67*x^4 - 57*x^3 - 33*x^2 + 97*x + 3)/(4*(x^5 + x^4 + x^3 - 3*x^2 + 3*x ...` | elem | ELEM | 0.043 | ELEM | 0.044 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-06 | deriv-alg | `(14*x^5 + 46*x^4 + 18*x^3 - 10*x^2 + 24*x - 4)/(3*(x^4 + 2*x^3 - x^2 + 2*x))*(x^4 + 2*x^3 ...` | elem | ELEM | 0.041 | ELEM | 0.040 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-07 | deriv-alg | `(-48*x^5 + 42*x^4 - 24*x + 6)/(6*(3*x^4 - 1))*(3*x^4 - 1)^(-5/6)` | elem | ELEM | 0.050 | ELEM | 0.050 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-08 | deriv-alg | `(-21*x^4 - 67*x^3 - 10*x^2 + 16*x - 10)/(5*(x^4 + 3*x^3 - 2*x^2 + 2))*(x^4 + 3*x^3 - 2*x^2...` | elem | ELEM | 0.035 | ELEM | 0.034 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-09 | deriv-alg | `(12*x^5 + 15*x^4 + 8*x^3 - 19*x^2 - 4*x)/(3*(3*x^4 + 3*x^3 + x^2 - 2*x))*(3*x^4 + 3*x^3 + ...` | elem | ELEM | 0.049 | ELEM | 0.047 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-10 | deriv-alg | `(6*x^3 - 21*x^2 + 24*x + 39)/(6*(x^3 - 3*x^2 + 3*x + 3))*(x^3 - 3*x^2 + 3*x + 3)^(-1/6)` | elem | ELEM | 0.092 | ELEM | 0.092 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-11 | deriv-alg | `(-3*x^4 + 36*x^3 + 27*x^2 - 66*x - 39)/(5*(x^3 - 2*x^2 - 2*x + 3))*(x^3 - 2*x^2 - 2*x + 3)...` | elem | ELEM | 0.085 | ELEM | 0.085 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-12 | deriv-alg | `(48*x^6 + 13*x^5 + 42*x^4 + 78*x^3 + 2*x^2 - 55*x - 12)/(3*(-x^5 - 2*x^3 - 2*x^2 - x + 2))...` | elem | ELEM | 0.035 | ELEM | 0.036 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-13 | deriv-alg | `(20*x^3 - 24*x^2 - 24*x + 4)/(4*(x^2 - 2*x))*(x^2 - 2*x)^(1/4)` | elem | ELEM | 0.031 | ELEM | 0.031 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-14 | deriv-alg | `(-18*x^3 - 33*x^2 + 36*x + 17)/(4*(3*x^2 + 3*x - 1))*(3*x^2 + 3*x - 1)^(-3/4)` | elem | ELEM | 0.044 | ELEM | 0.046 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-15 | deriv-alg | `(102*x^6 - 31*x^5 + 10*x^4 + 12*x^3 - 70*x^2 - 65*x - 18)/(6*(3*x^5 - 2*x^4 + x^2 - 3*x - ...` | elem | ELEM | 0.093 | ELEM | 0.094 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-16 | deriv-alg | `(-18*x^3 + 7*x^2 - 28*x - 13)/(2*(x^2 - x + 3))*(x^2 - x + 3)^(1/2)` | elem | ELEM | 0.027 | ELEM | 0.029 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-17 | deriv-alg | `(8*x^3 + 25*x^2 + 12*x)/(5*(-x^2 - 2*x))*(-x^2 - 2*x)^(-3/5)` | elem | ELEM | 0.051 | ELEM | 0.050 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-18 | deriv-alg | `(-12*x^5 - 10*x^4 - 30*x^3 - 32*x^2 - 6*x + 22)/(4*(-3*x^4 - 2*x^3 - 2*x^2 - 2*x + 2))*(-3...` | elem | ELEM | 0.062 | ELEM | 0.061 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-19 | deriv-alg | `(24*x^3 - 36*x^2 - 36*x + 6)/(3*(x^4 - 2*x^3 - 3*x^2 + x + 2))*(x^4 - 2*x^3 - 3*x^2 + x + ...` | elem | ELEM | 0.041 | ELEM | 0.041 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| deriv-alg-20 | deriv-alg | `(16*x + 12)/(5*(2*x^2 + 3*x + 2))*(2*x^2 + 3*x + 2)^(-4/5)` | elem | ELEM | 0.051 | ELEM | 0.051 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-01 | deriv-log | `(-27*x^5 - 45*x^4 + 126*x^3 + 9*x^2 - 21*x)/(2*(x^3 + x^2 - 3*x + 1)*(9*x^6 + 18*x^5 + 27*...` | elem | ELEM | 0.167 | DECLINE | 0.029 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-02 | deriv-log | `(-24*x^6 + 6*x^5 - 4*x^4 + 20*x^3 + 24*x^2 - 26*x - 12)/(2*(x^6 + x^5 + 2*x^4 + 2*x^3 - x^...` | elem | ELEM | 0.182 | ELEM | 0.043 | ELEM | 0.280 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-03 | deriv-log | `(9*x^5 + 3*x^4 + 9*x^3 - 18*x^2 + 6*x - 9)/(2*(x^3 + x)*(x^6 + 2*x^5 + 5*x^4 + x^3 + 10*x^...` | elem | ELEM | 0.151 | TIMEOUT | 60.289 | ELEM | 0.120 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-04 | deriv-log | `(-2*x^6 - 9*x^5 + 81*x^4 - 74*x^3 + 56*x^2 - 57*x + 13)/(2*(x^3 - 3*x^2 + x - 1)*(4*x^6 - ...` | elem | ELEM | 0.184 | TIMEOUT | 60.279 | ELEM | 0.150 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-05 | deriv-log | `(4*x^9 - 26*x^8 + 23*x^7 - 19*x^6 + 71*x^5 + 82*x^4 + 25*x^3 + 12*x^2 - 18*x + 15)/(2*(x^6...` | elem | ELEM | 0.273 | TIMEOUT | 60.234 | ELEM | 0.490 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-06 | deriv-log | `(12*x^5 - 36*x^4 + 6*x^3 - 33*x^2 + 54*x - 36)/(2*(x^4 + x^2 - 3*x + 3)*(-9*x^6 - 9*x^4 + ...` | elem | ELEM | 0.171 | TIMEOUT | 60.322 | ELEM | 0.250 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-07 | deriv-log | `(-15*x^6 + 33*x^5 - 15*x^4 + 69*x^3 + 16*x^2 - 6*x - 6)/(2*(x^4 - 3*x^3 + 2*x^2 + 2*x + 3)...` | elem | ELEM | 0.175 | TIMEOUT | 60.223 | ELEM | 0.260 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-08 | deriv-log | `(4*x^8 - 3*x^7 - 31*x^6 - 15*x^5 - 27*x^4 + 4*x^3 + 18*x^2 - 36*x + 24)/(2*(x^6 + x^5 + 3*...` | elem | ELEM | 0.291 | TIMEOUT | 60.226 | ELEM | 0.480 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-09 | deriv-log | `(8*x^4 - 48*x^3 - 28*x - 4)/(2*(2*x^4 + 2*x^2 + 2*x - 2)*(-8*x^4 - 7*x^2 - 14*x + 17))*Sqr...` | elem | ELEM | 0.128 | ELEM | 0.036 | ELEM | 0.160 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-10 | deriv-log | `(-4*x^4 + 16*x^3 + 6*x - 4)/(2*(2*x^4 + 2*x^2 - 2*x)*(-2*x^4 - x^2 - 2*x + 4))*Sqrt[2*x^4 ...` | elem | ELEM | 0.115 | ELEM | 0.035 | ELEM | 0.150 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-11 | deriv-log | `(-4*x^6 - 4*x^5 + 42*x^4 + 32*x^3 - 14*x^2)/(2*(-2*x^3 - 3*x^2 + x - 1)*(x^6 + 12*x^5 + 22...` | elem | ELEM | 0.127 | TIMEOUT | 60.249 | ELEM | 0.110 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-12 | deriv-log | `(-8*x^6 - 4*x^5 - 4*x^3)/(2*(x^4 + x^3 + 2*x^2 - x)*(-4*x^6 - 4*x^5 - 4*x^4 + 4*x^3))*Sqrt...` | elem | ELEM | 0.115 | DECLINE | 2.059 | ELEM | 0.130 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-13 | deriv-log | `(-32*x^8 - 44*x^7 + 32*x^6 + 147*x^5 + 72*x^4 - 45*x^3 + 22*x^2 + 46*x - 2)/(2*(-2*x^6 + 3...` | elem | ELEM | 1.073 | TIMEOUT | 60.250 | ELEM | 0.540 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-14 | deriv-log | `(6*x^5 - 12*x^4 - 2*x^3 - 4*x + 12)/(2*(-x^5 + 3*x^4 + x^3 - 2*x^2 - 2*x + 3)*(4*x^5 - 12*...` | elem | ELEM | 0.179 | ELEM | 0.048 | ELEM | 0.160 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-15 | deriv-log | `(-6*x^7 - 21*x^6 - 18*x^4 - 48*x^3 + 12*x - 3)/(2*(x^6 - x^5 + 3*x^4 + x^3 - x^2 + x)*(-9*...` | elem | ELEM | 0.195 | TIMEOUT | 60.258 | ELEM | 0.370 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-16 | deriv-log | `(-9*x^5 + 27*x^4 - 43*x^3 + 27*x^2 + 2*x - 8)/(2*(x^3 - 2*x^2 + 2*x - 1)*(9*x^6 - 18*x^5 +...` | elem | ELEM | 0.159 | TIMEOUT | 60.273 | ELEM | 0.170 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-17 | deriv-log | `(6*x^3 + 4*x + 4)/(2*(-3*x^3 - x^2 + 2*x + 1)*(3*x^3 + 5*x^2 - 2*x - 1))*Sqrt[-3*x^3 - x^2...` | elem | ELEM | 0.116 | ELEM | 0.041 | ELEM | 0.110 | AGREE-ELEM | AGREE-ELEM |
| deriv-log-18 | deriv-log | `(-18*x^4 + 117*x^2 + 126*x - 18)/(2*(x^4 + 3*x^3 + 3*x^2 - 2*x - 3)*(-5*x^4 - 27*x^3 - 15*...` | elem | ELEM | 0.129 | TIMEOUT | 60.246 | ELEM | 0.180 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-19 | deriv-log | `(-9*x^8 - 3*x^7 - 12*x^6 + 57*x^5 + 51*x^4 - 123*x^3 + 60*x^2 + 18*x - 11)/(2*(x^5 + x^4 +...` | elem | ELEM | 0.362 | TIMEOUT | 60.254 | ELEM | 0.360 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-20 | deriv-log | `(9*x^7 - 13*x^6 - 5*x^5 - 18*x^4 + 3*x^3 + 5*x^2 - 3*x - 2)/(2*(x^5 - 3*x^4 - x - 1)*(-x^7...` | elem | ELEM | 0.309 | TIMEOUT | 60.243 | ELEM | 0.220 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-21 | deriv-log | `(-4*x^6 + 12*x^5 - 18*x^4 - 29*x^3 - 16*x^2 - 2*x - 5)/(2*(2*x^4 + 3*x^2 + 3*x + 2)*(2*x^6...` | elem | ELEM | 0.179 | TIMEOUT | 60.298 | ELEM | 0.230 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-22 | deriv-log | `(6*x^6 + 6*x^5 - 2*x^4 + 7*x^3 + 9*x^2 + 4*x)/(2*(x^4 + 2*x^3 - 2*x^2 - x - 1)*(-9*x^6 - 2...` | elem | ELEM | 0.163 | TIMEOUT | 60.278 | ELEM | 0.270 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-23 | deriv-log | `(15*x^6 - 93*x^5 + 173*x^4 - 183*x^3 + 21*x^2 - 13*x + 6)/(2*(-x^5 + 2*x^4 - 3*x^3 - x^2 -...` | elem | ELEM | 0.207 | TIMEOUT | 60.227 | ELEM | 0.330 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-24 | deriv-log | `(12*x^8 - 40*x^7 - 104*x^6 - 30*x^5 + 45*x^4 + 21*x^3 + 3*x^2 - 3*x)/(2*(2*x^5 + x^4 - x^3...` | elem | ELEM | 0.222 | TIMEOUT | 60.260 | ELEM | 0.280 | AGREE-ELEM | FRICAS-ONLY |
| deriv-log-25 | deriv-log | `(-24*x^8 + 5*x^6 - 16*x^5 + 45*x^4 - 8*x^3 - 71*x^2 + 120*x - 3)/(2*(x^6 - x^5 - x^4 + 2*x...` | elem | ELEM | 0.297 | TIMEOUT | 60.242 | ELEM | 0.490 | AGREE-ELEM | FRICAS-ONLY |
| random-001 | random | `(2*x^2)/(-2*x^2 + 3*x - 1)*(x + 2)^(-1/2)` |  | ELEM | 0.109 | ELEM | 0.078 | ELEM | 0.100 | AGREE-ELEM | AGREE-ELEM |
| random-002 | random | `(3*x + 3)/(-x - 2)*(3*x^3 + x^2 + 3*x + 2)^(2/3)` |  | NONELEM | 0.028 | NONELEM | 0.027 | TIMEOUT | 60.044 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-003 | random | `(3*x + 2)/(-2*x + 3)*(x^4 + 2*x^2 - x + 3)^(-3/4)` |  | NONELEM | 0.504 | DECLINE | 2.872 | TIMEOUT | 60.046 | NONELEM-UNCONFIRMED | NEITHER |
| random-004 | random | `(3*x + 1)/(x^2 - 2*x - 2)*(x^3 + 2*x^2 - 1)^(1/3)` |  | NONELEM | 0.103 | DECLINE | 0.065 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| random-005 | random | `(x^2 + 1)/(2*x + 2)*Sqrt[x^4 - 3*x^3 + 3*x^2 + x + 2]` |  | NONELEM | 0.027 | NONELEM | 0.027 | NONELEM | 0.230 | AGREE-NONELEM | AGREE-NONELEM |
| random-006 | random | `(-2)/(-3)*(x^2 - 2*x)^(2/3)` |  | NONELEM | 0.022 | NONELEM | 0.022 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-007 | random | `(-1)/(-2)*(x^2 + 3*x - 2)^(2/3)` |  | NONELEM | 0.024 | NONELEM | 0.024 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-008 | random | `1/(-3)*(-x^2 - x + 1)^(-1/2)` |  | ELEM | 0.047 | ELEM | 0.052 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| random-009 | random | `2/(-3*x)*(x - 1)^(-1/4)` |  | ELEM | 0.079 | DECLINE | 0.076 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| random-010 | random | `(-3*x^2 - 3*x)/(-2*x^2 + 2*x - 1)*(-2*x^4 + 3*x^3 - 3*x^2 + 3*x + 2)^(-5/6)` |  | NONELEM | 21.354 | TIMEOUT | 60.373 | TIMEOUT | 60.031 | NONELEM-UNCONFIRMED | TIMEOUT |
| random-011 | random | `(-1)/(x + 2)*(-x^4 - 3*x^3 - x^2 - 1)^(-1/2)` |  | NONELEM | 0.073 | DECLINE | 0.209 | NONELEM | 0.160 | AGREE-NONELEM | NEITHER |
| random-012 | random | `(3*x - 3)/(2*x + 1)*(x^2 + x)^(-5/6)` |  | NONELEM | 3.141 | DECLINE | 0.053 | TIMEOUT | 60.045 | NONELEM-UNCONFIRMED | NEITHER |
| random-013 | random | `(x^2 - 3*x + 2)/3*(x^2 + x + 3)^(1/4)` |  | NONELEM | 0.026 | NONELEM | 0.025 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-014 | random | `(3*x^2 + 3*x + 2)/(-2*x + 1)*(x^2 + x + 3)^(4/5)` |  | NONELEM | 0.035 | NONELEM | 0.034 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-015 | random | `(2*x^2 - x)/(-2)*Sqrt[x^4 - x^3 - x]` |  | NONELEM | 0.025 | NONELEM | 0.026 | NONELEM | 0.150 | AGREE-NONELEM | AGREE-NONELEM |
| random-016 | random | `3/(-3)*(x + 1)^(-1/2)` |  | ELEM | 0.019 | ELEM | 0.019 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| random-017 | random | `2*(x + 2)^(1/3)` |  | ELEM | 0.024 | ELEM | 0.023 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| random-018 | random | `2/(-2)*(-2*x)^(3/4)` |  | ELEM | 0.024 | ELEM | 0.024 | ELEM | 0.070 | AGREE-ELEM | AGREE-ELEM |
| random-019 | random | `(-2*x^2 - 3*x + 2)/(-2)*(x^3 - x^2 + x + 1)^(-1/6)` |  | NONELEM | 0.034 | NONELEM | 0.034 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| random-020 | random | `1/(-2*x^2 - 2*x)*(-3*x^2 + 3*x - 2)^(2/3)` |  | NONELEM | 0.023 | NONELEM | 0.023 | TIMEOUT | 60.044 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-021 | random | `(-3)/2*x^(-1/3)` |  | ELEM | 0.020 | ELEM | 0.019 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| random-022 | random | `(x^2 - 2*x - 2)/(2*x^2 - x + 3)*(-x^4 - 3*x^3 + 3*x^2 - 3*x - 2)^(-1/2)` |  | NONELEM | 0.164 | TIMEOUT | 60.280 | TIMEOUT | 60.029 | NONELEM-UNCONFIRMED | TIMEOUT |
| random-023 | random | `(-x^2 + 2*x + 1)/(-3)*(x^3 + 2*x^2 + 2*x)^(1/4)` |  | NONELEM | 0.028 | NONELEM | 0.030 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-024 | random | `(2*x^2 - 3*x)/(-1)*(x + 2)^(1/4)` |  | ELEM | 0.032 | ELEM | 0.032 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| random-025 | random | `(3*x^2 - x - 2)/(-x + 2)*(x^3 - 3*x^2 + 2*x - 3)^(-1/2)` |  | NONELEM | 0.022 | NONELEM | 0.022 | NONELEM | 0.180 | AGREE-NONELEM | AGREE-NONELEM |
| random-026 | random | `(-3*x^2 + x + 2)/2*(x^3 - x^2 - 3*x + 2)^(-1/5)` |  | NONELEM | 0.032 | NONELEM | 0.032 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| random-027 | random | `(-2)/(-x^2 + 2*x - 3)*(3*x + 2)^(2/3)` |  | ELEM | 0.120 | DECLINE | 0.152 | ELEM | 0.170 | AGREE-ELEM | FRICAS-ONLY |
| random-028 | random | `(x^2 - 2*x - 2)/(x^2 - 2*x - 2)*(-3*x^3 - 2*x^2 - x - 3)^(-3/4)` |  | NONELEM | 0.025 | NONELEM | 0.026 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-029 | random | `(3*x^2 - 2*x + 3)/x*(2*x^3 + 2*x^2 + 3*x + 3)^(2/5)` |  | NONELEM | 0.034 | NONELEM | 0.034 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-030 | random | `(-3*x^2 + 3*x - 3)/(2*x^2 - x - 1)*(x + 3)^(1/4)` |  | ELEM | 0.159 | DECLINE | 0.271 | ELEM | 0.120 | AGREE-ELEM | FRICAS-ONLY |
| random-031 | random | `(-1)/(-3*x - 1)*x^(-2/3)` |  | ELEM | 0.046 | ELEM | 0.032 | ELEM | 0.080 | AGREE-ELEM | AGREE-ELEM |
| random-032 | random | `(2*x + 3)/(-3*x^2 - 2*x + 3)*(x^2 - x + 2)^(1/4)` |  | NONELEM | 0.240 | TIMEOUT | 60.314 | TIMEOUT | 60.043 | NONELEM-UNCONFIRMED | TIMEOUT |
| random-033 | random | `3/(-2*x^2 + 3*x - 3)*(-3*x^2 + 2*x)^(1/4)` |  | NONELEM | 0.160 | DECLINE | 0.359 | TIMEOUT | 60.046 | NONELEM-UNCONFIRMED | NEITHER |
| random-034 | random | `(-3*x - 1)/(-2)*Sqrt[x + 3]` |  | ELEM | 0.029 | ELEM | 0.026 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| random-035 | random | `(x^2 - 3*x + 3)/(-1)*(x^4 + x^3 - x^2 + 3)^(-1/2)` |  | NONELEM | 0.021 | NONELEM | 0.022 | NONELEM | 0.200 | AGREE-NONELEM | AGREE-NONELEM |
| random-036 | random | `3/(2*x^2 + 3*x)*(x + 1)^(2/3)` |  | ELEM | 0.203 | DECLINE | 0.055 | ELEM | 0.100 | AGREE-ELEM | FRICAS-ONLY |
| random-037 | random | `(-x^2 - 3*x - 3)/(2*x + 1)*(x^4 + 3*x^3 + 3*x^2 - 3*x - 2)^(-5/6)` |  | NONELEM | 7.280 | DECLINE | 2.094 | TIMEOUT | 60.037 | NONELEM-UNCONFIRMED | NEITHER |
| random-038 | random | `(-1)/(-3*x + 1)*(-2*x^4 - 2*x^3 + x^2 - x - 3)^(-3/4)` |  | NONELEM | 0.397 | DECLINE | 1.230 | TIMEOUT | 60.051 | NONELEM-UNCONFIRMED | NEITHER |
| random-039 | random | `(2*x + 1)/(-1)*(x^2 + x + 1)^(-1/2)` |  | ELEM | 0.023 | ELEM | 0.023 | ELEM | 0.060 | AGREE-ELEM | AGREE-ELEM |
| random-040 | random | `2/(2*x^2 + 3*x + 2)*(x^3 - 2*x^2 + x + 1)^(-2/3)` |  | NONELEM | 0.128 | DECLINE | 2.356 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
| random-041 | random | `(-3)/(3*x^2 - 2*x)*(x^3 - 3*x^2 - 2*x - 3)^(3/4)` |  | NONELEM | 0.030 | NONELEM | 0.030 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-042 | random | `(x^2 + 2*x + 3)/(-2*x + 2)*(x^3 + x - 3)^(-1/2)` |  | NONELEM | 0.022 | NONELEM | 0.022 | NONELEM | 0.160 | AGREE-NONELEM | AGREE-NONELEM |
| random-043 | random | `(-3*x^2 - 2)/2*(x^3 + 2*x^2 + 3*x + 2)^(-2/3)` |  | NONELEM | 0.089 | DECLINE | 0.068 | TIMEOUT | 60.025 | NONELEM-UNCONFIRMED | NEITHER |
| random-044 | random | `(-3)/(-3)*Sqrt[x^4 + 3*x^3 - 2*x^2 + 2*x + 3]` |  | NONELEM | 0.024 | NONELEM | 0.024 | NONELEM | 0.160 | AGREE-NONELEM | AGREE-NONELEM |
| random-045 | random | `(-2)/2*(-3*x^4 - 2*x^3 - x^2 - x + 2)^(3/5)` |  | NONELEM | 0.034 | NONELEM | 0.034 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-046 | random | `(-3)/(-3)*(-3*x^4 + x^2 + 2*x + 2)^(-1/2)` |  | NONELEM | 0.019 | NONELEM | 0.020 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-047 | random | `(-2*x + 1)/(-1)*x^(-2/3)` |  | ELEM | 0.024 | ELEM | 0.025 | ELEM | 0.050 | AGREE-ELEM | AGREE-ELEM |
| random-048 | random | `(-3*x^2 + 3*x + 2)/(3*x^2 - x + 2)*(2*x + 1)^(-1/2)` |  | ELEM | 0.076 | ELEM | 0.071 | BAD | 0.130 | MATHILDA-ONLY | MATHILDA-ONLY |
| random-049 | random | `(x + 1)/2*Sqrt[-2*x^3 + 2*x^2 - 2]` |  | NONELEM | 0.023 | NONELEM | 0.023 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-050 | random | `(2*x - 3)/(-2)*(2*x^2 - 3*x - 2)^(1/4)` |  | NONELEM | 0.024 | NONELEM | 0.025 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-051 | random | `(2*x - 1)/(-3*x^2 + 2*x - 1)*(x^3 - 2*x^2 + 3*x - 2)^(2/3)` |  | NONELEM | 0.027 | NONELEM | 0.027 | ERROR | - | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-052 | random | `(-2)/3*(x^3 + 3*x^2 - 3*x - 3)^(1/4)` |  | NONELEM | 0.024 | NONELEM | 0.024 | NONELEM | 0.100 | AGREE-NONELEM | AGREE-NONELEM |
| random-053 | random | `(-3*x + 2)/(2*x - 2)*(3*x)^(-4/5)` |  | ELEM | 0.074 | DECLINE | 1.693 | ELEM | 0.120 | AGREE-ELEM | FRICAS-ONLY |
| random-054 | random | `(-3*x^2 + 2*x + 1)/(-3)*(x^4 - x^3 - 3*x^2 - 2*x - 2)^(1/6)` |  | NONELEM | 0.036 | NONELEM | 0.035 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-055 | random | `(x^2 + 1)/(-3*x^2 - 2*x - 1)*(x^3 + 3*x - 2)^(5/6)` |  | NONELEM | 0.067 | NONELEM | 0.059 | NONELEM | 0.090 | AGREE-NONELEM | AGREE-NONELEM |
| random-056 | random | `(x + 3)/(-3)*(2*x^4 + 2*x^3 + x^2 + 2*x + 3)^(-2/3)` |  | NONELEM | 0.023 | NONELEM | 0.023 | NONELEM | 0.080 | AGREE-NONELEM | AGREE-NONELEM |
| random-057 | random | `(2*x + 2)/(-3*x - 1)*(x^3 + x^2 + x)^(2/3)` |  | NONELEM | 0.027 | NONELEM | 0.027 | TIMEOUT | 60.041 | NONELEM-UNCONFIRMED | NONELEM-UNCONFIRMED |
| random-058 | random | `(-3*x - 3)/2*(x^4 + x^3 + 3*x^2 - 3*x + 2)^(-1/2)` |  | NONELEM | 0.079 | DECLINE | 0.098 | NONELEM | 0.140 | AGREE-NONELEM | NEITHER |
| random-059 | random | `(-3)/(3*x + 2)*(x^3 + x^2 + 2*x - 3)^(-1/2)` |  | NONELEM | 0.060 | DECLINE | 0.154 | NONELEM | 1.320 | AGREE-NONELEM | NEITHER |
| random-060 | random | `(2*x^2 - 1)/(2*x^2 - 3*x + 2)*(-2*x^2 - 2*x + 1)^(-2/3)` |  | NONELEM | 0.113 | DECLINE | 0.379 | ERROR | - | NONELEM-UNCONFIRMED | NEITHER |
