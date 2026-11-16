# 258 — parametric A/B of the Trager and Bronstein cascade stages

1661 cases; A = `Integrate`$CascadeParameters = False`, B = `{"Trager"}`; timeout 30 s, check at a = 1/7.

- hard failures in B: **1**
- B worse than A: **10**
- B better than A: **239**

| A | B | cases |
|---|---|---|
| ELEM | ELEM | 772 |
| UNEVAL | UNEVAL | 385 |
| TIMEOUT | TIMEOUT | 272 |
| UNEVAL | CERT | 167 |
| TIMEOUT | CERT | 20 |
| UNEVAL | ELEM | 15 |
| TIMEOUT | ELEM | 13 |
| SPECIAL | SPECIAL | 5 |
| UNEVAL | TIMEOUT | 5 |
| ELEM | TIMEOUT | 3 |
| WRONG | ELEM | 2 |
| UNVERIFIED | UNVERIFIED | 1 |
| WRONG | WRONG | 1 |

## Hard failures

- `trager-cheb-002` (unknown): WRONG — `x^2*Sqrt[-1 + a + 2*x^4]`

## Worse

- `genus1-1-2` (elem): lost ELEM -> TIMEOUT — `(3*x*Log[1 + x^2]^2)/(2 + a + Log[1 + x^2]^3 + 2*Sqrt[2 + a + Log[1 + x^2]^3] + 2*x^2*Sqrt[2 + a + Log[1 + x^2]^3] + x^2*(2 + a + Log[1 + x^2]^3))`
- `genus1-1-8` (elem): lost ELEM -> TIMEOUT — `(3*x*Log[1 + x^2]^2)/(3 + a + Log[1 + x^2]^3 + Sqrt[3 + a + Log[1 + x^2]^3] + x^2*Sqrt[3 + a + Log[1 + x^2]^3] + x^2*(3 + a + Log[1 + x^2]^3))`
- `trager-evensub-005` (unknown): larger 107 -> 280 — `x/(Sqrt[1 + a + x^4]*(1 + x^2 + x^4))`
- `genus1-3-12` (elem): slower 30.0s -> 148.8s — `(3*E^x*Log[1 + E^x]^2)/(6 + 6*E^x + 2*Log[1 + E^x]^3 + 2*E^x*Log[1 + E^x]^3 + 4*Sqrt[3 + Log[1 + E^x]^3] + 4*E^x*Sqrt[3 + Log[1 + E^x]^3]) + (E^x*Sqrt[3 + Log[1`
- `algres-3-2` (elem): slower 30.0s -> 246.2s — `(-2*x*Sqrt[2 + Log[1 + x^2]^2]*(2 + Log[1 + x^2]^2 + Sqrt[2 + Log[1 + x^2]^2] + Log[1 + x^2]*Sqrt[2 + Log[1 + x^2]^2]) + 2*x*(Log[1 + x^2]^2 + Log[1 + x^2]^3 + `
- `algres-3-10` (elem): slower 30.0s -> 244.0s — `((-1 + Log[x] + Sqrt[2 + a + Log[x]^2])*((Log[x] + Sqrt[2 + Log[x]^2])/(-(x*Sqrt[2 + Log[x]^2]) + x*Log[x]*Sqrt[2 + Log[x]^2] + x*Sqrt[2 + Log[x]^2]*Sqrt[2 + a `
- `ratb-3-4` (elem): slower 30.1s -> 246.7s — `(6 + 2*a + 18*x + 6*a*x + 3*Log[1 + 3*x] + a*Log[1 + 3*x] + 9*x*Log[1 + 3*x] + 3*a*x*Log[1 + 3*x] - 12*(1 + Log[1 + (3 + a)*x]) - 36*x*(1 + Log[1 + (3 + a)*x]) `
- `cA-1-263` (elem): slower 0.1s -> 10.1s — `(2*a*x^2 - 2*a*x^2*(1 + x^3)^(1/6) - 2*a*x^2*(1 + x^3)^(1/3) - a^2*x^2*(1 + x^3)^(1/3) + 4*a*x^2*Sqrt[1 + x^3] + a^2*x^2*Sqrt[1 + x^3] - 2*a*x^2*(1 + x^3)^(2/3)`
- `cM-1-292` (elem): lost ELEM -> TIMEOUT — `(2*(1 + Log[1 + 2*x])^(5/6) + 4*x*(1 + Log[1 + 2*x])^(5/6) - 2*(1 + Log[1 + 2*x]) - 4*x*(1 + Log[1 + 2*x]) + 2*(1 + Log[1 + 2*x])^(7/6) + 4*x*(1 + Log[1 + 2*x])`
- `cM-2-67` (unknown): larger 78 -> 179 — `-((E^x*(-2 + E^x))/(Sqrt[2 + a + E^x]*(-1 + E^x - 2*E^(2*x))))`

## Better

- `agen-1-5` (elem): UNEVAL -> ELEM — `(-4 + 14*x - 2*a*x^2 - 15*x^3 + 4*a*x^3 + 27*x^4 + 9*a*x^4 + 6*Sqrt[-3 - 2*x + 2*x^2 + 3*x^3] + 2*a*Sqrt[-3 - 2*x + 2*x^2 + 3*x^3] - 4*x*Sqrt[-3 - 2*x + 2*x^2 +`
- `agen-1-11` (elem): UNEVAL -> ELEM — `(2 + 8*x + 9*x^2 + 3*x^3 - 2*(-2 + a + 2*x + 3*x^2 + x^3))/(2*Sqrt[-2 + a + 2*x + 3*x^2 + x^3] + 4*x*Sqrt[-2 + a + 2*x + 3*x^2 + x^3] + 2*x^2*Sqrt[-2 + a + 2*x `
- `bbottomp-1-12` (unknown): TIMEOUT -> CERT — `E^x*Sqrt[x] - E^x/(2*(2 + x)^(3/2)) + 2*E^x*Sqrt[2 + x] + (2*E^x)/Sqrt[2 + a + x]`
- `trager-binom-003` (unknown): UNEVAL -> CERT — `x^3/Sqrt[2 + a + x^6]`
- `trager-binom-005` (unknown): UNEVAL -> CERT — `x^2/(1 + a + x^4)^(1/3)`
- `trager-binom-007` (unknown): UNEVAL -> CERT — `x^3/(1 + a + x^6)^(1/3)`
- `trager-binom-008` (unknown): UNEVAL -> CERT — `x^3/(-1 + a + x^3)^(1/4)`
- `trager-binom-009` (unknown): UNEVAL -> CERT — `x/(-1 + a + x^4)^(1/3)`
- `trager-binom-010` (unknown): UNEVAL -> CERT — `x^2/(2 + a + x^2)^(1/3)`
- `trager-binom-015` (unknown): UNEVAL -> CERT — `x^3/(2 + a + x^6)^(1/3)`
- `trager-binom-019` (unknown): UNEVAL -> CERT — `1/(-1 + a + x^2)^(1/3)`
- `trager-binom-020` (unknown): UNEVAL -> CERT — `x^3/(-1 + a + x^6)^(1/4)`
- `trager-binom-021` (unknown): UNEVAL -> CERT — `1/(1 + a + x^2)^(1/4)`
- `trager-binom-023` (unknown): UNEVAL -> CERT — `x/(2 + a + x^3)^(1/3)`
- `trager-binom-024` (unknown): UNEVAL -> CERT — `1/(-1 + a + x^6)^(1/3)`
- `trager-binom-025` (unknown): UNEVAL -> CERT — `x^3/Sqrt[-3 + a + x^3]`
- `trager-binom-026` (unknown): UNEVAL -> CERT — `1/Sqrt[-1 + a + x^6]`
- `trager-binom-027` (unknown): UNEVAL -> CERT — `x/(-1 + a + x^4)^(1/4)`
- `trager-binom-033` (unknown): UNEVAL -> CERT — `x^5/(2 + a + x^4)^(1/4)`
- `trager-binom-035` (unknown): UNEVAL -> CERT — `1/(1 + a + x^6)^(1/4)`
- `trager-binom-038` (unknown): UNEVAL -> CERT — `x/(-3 + a + x^4)^(1/4)`
- `trager-binom-039` (unknown): UNEVAL -> CERT — `x^3/Sqrt[-3 + a + x^6]`
- `trager-binom-040` (unknown): UNEVAL -> CERT — `1/(-3 + a + x^6)^(1/3)`
- `trager-binom-041` (unknown): UNEVAL -> CERT — `x^2/(2 + a + x^2)^(1/4)`
- `trager-binom-042` (unknown): WRONG -> ELEM — `x^2/Sqrt[-3 + a + x^6]`
- `trager-binom-044` (unknown): UNEVAL -> CERT — `1/Sqrt[-3 + a + x^6]`
- `trager-binom-045` (unknown): WRONG -> ELEM — `x^5/(-3 + a + x^6)^(1/4)`
- `trager-third-012` (unknown): UNEVAL -> CERT — `1/((-1 + a + x)*Sqrt[1 + x^3])`
- `trager-third-013` (unknown): TIMEOUT -> CERT — `1/((-2 + a + x)*Sqrt[1 + x^3])`
- `trager-third-014` (unknown): TIMEOUT -> CERT — `1/((-3 + a + x)*Sqrt[1 + x^3])`
- `trager-third-015` (unknown): TIMEOUT -> CERT — `1/((-2 + a + x)*Sqrt[-x + x^3])`
- `trager-third-016` (unknown): UNEVAL -> CERT — `1/((-3 + a + x)*Sqrt[-x + x^3])`
- `trager-third-018` (unknown): UNEVAL -> CERT — `1/((-1 + a + x)*Sqrt[1 + x^4])`
- `trager-third-019` (unknown): UNEVAL -> CERT — `1/((1 + a + x)*Sqrt[1 + x^4])`
- `trager-third-020` (unknown): UNEVAL -> CERT — `1/((-2 + a + x)*Sqrt[1 + x^4])`
- `trager-third-021` (unknown): UNEVAL -> CERT — `1/((-3 + a + x)*Sqrt[1 + x^4])`
- `trager-third-023` (unknown): UNEVAL -> CERT — `1/((-2 + a + x)*Sqrt[-1 + x^4])`
- `trager-third-024` (unknown): UNEVAL -> CERT — `1/((-3 + a + x)*Sqrt[-1 + x^4])`
- `trager-third-026` (unknown): UNEVAL -> CERT — `1/((-1 + a + x)*Sqrt[1 + x^2 + x^4])`
- `trager-third-027` (unknown): UNEVAL -> CERT — `1/((1 + a + x)*Sqrt[1 + x^2 + x^4])`
- `trager-third-028` (unknown): UNEVAL -> CERT — `1/((-2 + a + x)*Sqrt[1 + x^2 + x^4])`
- `trager-third-029` (unknown): UNEVAL -> CERT — `1/((-3 + a + x)*Sqrt[1 + x^2 + x^4])`
- `trager-third-030` (unknown): UNEVAL -> CERT — `1/(x*Sqrt[1 + a + x + x^3])`
- `trager-ratrad-004` (unknown): UNEVAL -> CERT — `x^3/((-2 + a + x)*(3 + x^4)^(3/5))`
- `trager-ratrad-005` (unknown): UNEVAL -> CERT — `(x^3*(1 + x^3)^(3/2))/(-2 + a + x)`
- `trager-ratrad-006` (unknown): smaller 294 -> 60 — `(1 + a + x^3)^(2/3)/x`
- `trager-ratrad-007` (unknown): UNEVAL -> CERT — `(-2*x + x^3)^(1/4)/(1 + a + x)`
- `trager-ratrad-008` (unknown): TIMEOUT -> CERT — `(-1 + a + 2*x)/((1 + x)*(1 + x^2)^(2/3))`
- `trager-ratrad-011` (unknown): UNEVAL -> CERT — `x^2*(1 + a - x + x^2)^(1/5)`
- `trager-ratrad-012` (unknown): UNEVAL -> CERT — `x/((-2 + a + x)*(1 + x^3)^(2/3))`
- `trager-ratrad-013` (unknown): TIMEOUT -> CERT — `(1 + a - x + x^2)^(3/4)/x`
- `trager-ratrad-014` (unknown): UNEVAL -> CERT — `(x^3*Sqrt[-2*x + x^3])/(1 + a + x)`
- `trager-ratrad-015` (unknown): UNEVAL -> CERT — `(1 + a + x^2)/((-2 + x)*(-2*x + x^3)^(3/5))`
- `trager-ratrad-017` (unknown): smaller 155 -> 66 — `x*(3 + a + x^4)^(3/2)`
- `trager-ratrad-018` (unknown): UNEVAL -> CERT — `((-2 + a)*x + x^3)^(3/2)/x`
- `trager-ratrad-020` (unknown): UNEVAL -> CERT — `x^3/((-2 + a + x)*Sqrt[-2*x + x^3])`
- `trager-ratrad-021` (unknown): TIMEOUT -> CERT — `x^3/((-2 + x)*(1 + a - x + x^2)^(1/6))`
- `trager-ratrad-023` (unknown): UNEVAL -> CERT — `(1 + a + x^2)^(6/5)`
- `trager-ratrad-024` (unknown): UNEVAL -> CERT — `(-1 + (2 + a)*x)*(1 - x + x^2)^(1/5)`
- `trager-ratrad-029` (unknown): TIMEOUT -> CERT — `(1 + a + x^2)^(7/4)/(2 + x^2)`
- `trager-ratrad-034` (unknown): UNEVAL -> CERT — `x^2*(1 + a + x^2)^(5/6)`
- `trager-ratrad-035` (unknown): UNEVAL -> CERT — `(3 + a + x^4)^(5/6)`
- `trager-ratrad-036` (unknown): UNEVAL -> CERT — `((-1 + (2 + a)*x)*(3 + x^4)^(5/6))/(2 + x^2)`
- `trager-ratrad-037` (unknown): UNEVAL -> CERT — `1/((1 + a + x)*(-2*x + x^3)^(1/4))`
- `trager-ratrad-038` (unknown): UNEVAL -> CERT — `(1 + x^3)^(2/5)/(2 + a + x^2)`
- `trager-ratrad-039` (unknown): UNEVAL -> CERT — `((-1 + (2 + a)*x)*(-2*x + x^3)^(1/6))/x`
- `trager-ratrad-040` (unknown): UNEVAL -> CERT — `((-1 + a + 2*x)*(3 + x^4)^(2/5))/(1 + x)`
- `trager-ratrad-042` (unknown): UNEVAL -> CERT — `(1 + a + x^2)/(3 + x^4)^(1/3)`
- `trager-ratrad-045` (unknown): smaller 435 -> 48 — `1/(x*(1 + a + x^2)^(1/6))`
- `trager-ratrad-050` (unknown): TIMEOUT -> CERT — `(x*(1 + a - x + x^2)^(1/3))/(1 + x)`
- `trager-cheb-003` (unknown): UNEVAL -> CERT — `x^4/Sqrt[3 + a - x^3]`
- `trager-cheb-004` (unknown): UNEVAL -> CERT — `x/(-1 + a + 2*x^3)^(3/4)`
- `trager-cheb-005` (unknown): UNEVAL -> CERT — `(-1 + a + 2*x^4)^(1/4)`
- `trager-cheb-006` (unknown): UNEVAL -> CERT — `(3 + a - x^4)^(1/6)`
- `trager-cheb-007` (unknown): UNEVAL -> ELEM — `Sqrt[3 + a - x^3]/x`
- `trager-cheb-008` (unknown): UNEVAL -> CERT — `x^4*(3 + a - x^4)^(2/3)`
- `trager-cheb-012` (unknown): UNEVAL -> CERT — `(3 + a - x^2)^(1/4)`
- `trager-cheb-014` (unknown): smaller 599 -> 110 — `x^3*(-1 + a + 2*x^3)^(2/3)`
- `trager-cheb-018` (unknown): UNEVAL -> CERT — `(1 + a + x^3)^(2/5)/x^2`
- `trager-cheb-019` (unknown): UNEVAL -> CERT — `1/(2 + a + 3*x^2)^(1/4)`
- `trager-cheb-020` (unknown): UNEVAL -> ELEM — `(-1 + a + 2*x^4)^(2/3)/x`
- `trager-cheb-021` (unknown): UNEVAL -> CERT — `x*(2 + a + 3*x^4)^(1/3)`
- `trager-cheb-022` (unknown): UNEVAL -> CERT — `x*(-1 + a + 2*x^3)^(3/2)`
- `trager-cheb-023` (unknown): UNEVAL -> CERT — `(2 + a + 3*x^4)^(1/6)`
- `trager-cheb-024` (unknown): UNEVAL -> CERT — `x^4/(1 + a + x^3)^(3/4)`
- `trager-cheb-026` (unknown): UNEVAL -> CERT — `(3 + a - x^3)^(3/2)/x^3`
- `trager-cheb-027` (unknown): smaller 257 -> 50 — `1/(x*(3 + a - x^2)^(1/3))`
- `trager-cheb-028` (unknown): UNEVAL -> CERT — `x/(3 + a - x^3)^(5/6)`
- `trager-cheb-030` (unknown): TIMEOUT -> CERT — `x^4*(-1 + a + 2*x^2)^(1/4)`
- `trager-cheb-035` (unknown): UNEVAL -> CERT — `(2 + a + 3*x^3)^(1/4)/x^3`
- `trager-cheb-038` (unknown): smaller 405 -> 69 — `(3 + a - x^4)^(3/2)/x^3`
- `trager-cheb-039` (unknown): UNEVAL -> CERT — `x^4*(2 + a + 3*x^4)^(2/5)`
- `trager-cheb-040` (unknown): UNEVAL -> CERT — `x^4/(-1 + a + 2*x^4)^(1/3)`
- `trager-goursat-001` (unknown): UNEVAL -> CERT — `(-1 + x^2)/((1 + a + x^2)*Sqrt[1 + x^4])`
- `trager-goursat-002` (unknown): UNEVAL -> CERT — `(1 + a + x^2)/((-1 + x^2)*Sqrt[1 + x^4])`
- `trager-goursat-005` (unknown): UNEVAL -> CERT — `(1 + (-1 + a)*x^2)/((1 + x^2)*Sqrt[1 + x^2 + x^4])`
- `trager-goursat-007` (unknown): UNEVAL -> CERT — `x/((1 + x^3)^(1/3)*(-2 + a + x^3))`
- `trager-goursat-008` (unknown): UNEVAL -> CERT — `x/((8 + x^3)*(1 + a + x^3)^(1/3))`
- `trager-goursat-009` (unknown): UNEVAL -> CERT — `(-2 + x^3)/((4 + x^3)*(1 + a + x^3)^(2/3))`
- `trager-goursat-012` (unknown): smaller 122 -> 48 — `1/(x*(1 + a + x^4)^(1/4))`
- `trager-goursat-013` (unknown): UNEVAL -> CERT — `(-1 + a + x^4)/(x^2*(1 + x^4)^(1/4))`
- `trager-goursat-014` (unknown): UNEVAL -> CERT — `1/((1 + a + x^2)*(1 + x^4)^(1/4))`
- `trager-goursat-017` (unknown): UNEVAL -> ELEM — `(1 + a + 3*x^2)/((-1 + x^2)*Sqrt[-x^2 + x^4])`
- `trager-goursat-019` (unknown): UNEVAL -> CERT — `(1 + (5 + a)*x^4)/Sqrt[1 + x + x^5]`
- `trager-goursat-020` (unknown): UNEVAL -> ELEM — `(2 + x^2)/((1 + a + x^2)*(x + x^3)^(1/3))`
- `trager-goursat-021` (unknown): UNEVAL -> CERT — `1/((-2 + x^2)*(1 + a - 4*x^2 + x^4)^(1/4))`
- `trager-goursat-022` (unknown): UNEVAL -> CERT — `(-1 + a + x)/((1 + x)*Sqrt[x + x^2 + x^3])`
- `trager-goursat-023` (unknown): UNEVAL -> CERT — `(-2 + x^2)/((2 + a + x^2)*Sqrt[4 + x^4])`
- `trager-recip-005` (unknown): smaller 1120 -> 122 — `(-1 + x^2)/((1 + (-1 + a)*x + x^2)*Sqrt[1 + 3*x^2 + x^4])`
- `trager-recip-006` (unknown): UNEVAL -> CERT — `(-1 + x^2)/((1 + a + x^2)*Sqrt[1 + 5*x^2 + x^4])`
- `trager-recip-010` (unknown): UNEVAL -> ELEM — `(1 + x^2)/((-1 + (3 + a)*x + x^2)*Sqrt[1 + 3*x^2 + x^4])`
- `trager-recip-011` (unknown): UNEVAL -> ELEM — `(1 + x^2)/((-1 + (-1 + a)*x + x^2)*Sqrt[1 - x^2 + x^4])`
- `trager-hyper-005` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a - x + x^3]`
- `trager-hyper-006` (unknown): UNEVAL -> CERT — `x^2/Sqrt[1 + a - x + x^3]`
- `trager-hyper-007` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a + x + x^4]`
- `trager-hyper-008` (unknown): UNEVAL -> CERT — `x^2/Sqrt[1 + a + x + x^4]`
- `trager-hyper-009` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a - x + x^5]`
- `trager-hyper-010` (unknown): UNEVAL -> CERT — `x/Sqrt[1 + a - x + x^5]`
- `trager-hyper-011` (unknown): UNEVAL -> CERT — `x^3/Sqrt[1 + a - x + x^5]`
- `trager-hyper-012` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a + x^2 + x^6]`
- `trager-hyper-013` (unknown): UNEVAL -> CERT — `x/Sqrt[1 + a + x^2 + x^6]`
- `trager-hyper-014` (unknown): UNEVAL -> CERT — `x^3/Sqrt[1 + a + x^2 + x^6]`
- `trager-hyper-015` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a + x + x^7]`
- `trager-hyper-016` (unknown): UNEVAL -> CERT — `x^2/Sqrt[1 + a + x + x^7]`
- `trager-hyper-017` (unknown): UNEVAL -> CERT — `x^4/Sqrt[1 + a + x + x^7]`
- `trager-hyper-018` (unknown): UNEVAL -> CERT — `1/Sqrt[1 + a + x^3 + x^8]`
- `trager-hyper-019` (unknown): UNEVAL -> CERT — `x^2/Sqrt[1 + a + x^3 + x^8]`
- `trager-hyper-020` (unknown): UNEVAL -> CERT — `x^4/Sqrt[1 + a + x^3 + x^8]`
- `trager-super-001` (unknown): UNEVAL -> CERT — `1/(1 + a + x^4)^(1/3)`
- `trager-super-002` (unknown): UNEVAL -> CERT — `x/(1 + a + x^4)^(2/3)`
- `trager-super-003` (unknown): UNEVAL -> CERT — `1/(2 + a + x^5)^(1/3)`
- `trager-super-004` (unknown): UNEVAL -> CERT — `x/(2 + a + x^5)^(2/3)`
- `trager-super-005` (unknown): UNEVAL -> CERT — `1/(1 + a + x + x^3)^(1/3)`
- `trager-super-006` (unknown): UNEVAL -> CERT — `x/(1 + a + x + x^3)^(2/3)`
- `trager-super-007` (unknown): UNEVAL -> CERT — `1/(1 + a + x + x^4)^(1/5)`
- `trager-super-008` (unknown): UNEVAL -> CERT — `x/(1 + a + x + x^4)^(2/5)`
- `trager-super-009` (unknown): UNEVAL -> CERT — `1/(1 + a + x^2)^(1/3)`
- `trager-ratradicand-006` (unknown): UNEVAL -> CERT — `Sqrt[(1 + x^2)/(-1 + a + x^2)]`
- `trager-ratradicand-010` (unknown): UNEVAL -> CERT — `Sqrt[(-1 + a + x)/(1 + x^2)]`
- `trager-ratradicand-011` (unknown): UNEVAL -> CERT — `((1 + a + x^2)/x)^(1/3)`
- `trager-ratradicand-012` (unknown): UNEVAL -> CERT — `Sqrt[(1 + a + x^3)/x]`
- `trager-ratradicand-015` (unknown): UNEVAL -> CERT — `Sqrt[x/(1 + a + x + x^2)]`
- `trager-ratradicand-017` (unknown): UNEVAL -> CERT — `(x/(1 + a + x^3))^(1/3)`
- `trager-ratradicand-019` (unknown): TIMEOUT -> CERT — `((-1 + a + x^2)/(1 + x^2))^(1/4)`
- `trager-ratradicand-020` (unknown): UNEVAL -> CERT — `Sqrt[(1 + a + x^2)/x]/x`
- `trager-torsion-001` (unknown): UNEVAL -> CERT — `1/((-3 + a + x)*(1 + x^3)^(1/3))`
- `trager-torsion-003` (unknown): UNEVAL -> CERT — `1/((2 + a + x)*(2 + x^3)^(1/3))`
- `trager-torsion-004` (unknown): UNEVAL -> CERT — `1/((-1 + (2 + a)*x)*(-1 + x^3)^(1/3))`
- `trager-torsion-011` (unknown): UNEVAL -> CERT — `1/((-3 + x^2)*(1 + a - 6*x^2 + x^4)^(1/4))`
- `trager-torsion-012` (unknown): UNEVAL -> CERT — `1/((2 + a + x^2)*(4 + x^4)^(1/4))`
- `trager-torsion-018` (unknown): UNEVAL -> CERT — `(-3*x + x^3)^(1/4)/(-1 + a + x)`
- `trager-deriv-alg-02` (unknown): UNEVAL -> CERT — `(6 + a + 36*x + 27*x^2 - 32*x^3 - 10*x^4 - 4*x^5)/(3*(-2 + 3*x^2 + 2*x^3 - 2*x^4)^(5/3))`
- `trager-deriv-alg-03` (unknown): UNEVAL -> CERT — `(-22 + a - 22*x + 10*x^2 + 10*x^3)/(4*(2 - x - 2*x^2 + x^3)^(7/4))`
- `trager-deriv-alg-04` (unknown): UNEVAL -> CERT — `(12 + a + 16*x + 4*x^2)/(4*(-3 - 2*x + x^2)^(7/4))`
- `trager-deriv-alg-05` (unknown): UNEVAL -> CERT — `(3 + (97 + a)*x - 33*x^2 - 57*x^3 + 67*x^4 + 97*x^5 + 46*x^6)/(4*(1 + 3*x - 3*x^2 + x^3 + x^4 + x^5)^(1/4))`
- `trager-deriv-alg-06` (unknown): UNEVAL -> CERT — `(-4 + (24 + a)*x - 10*x^2 + 18*x^3 + 46*x^4 + 14*x^5)/(3*(2*x - x^2 + 2*x^3 + x^4)^(1/3))`
- `trager-deriv-alg-07` (unknown): UNEVAL -> CERT — `(6 + a - 24*x + 42*x^4 - 48*x^5)/(6*(-1 + 3*x^4)^(11/6))`
- `trager-deriv-alg-08` (unknown): UNEVAL -> CERT — `(-10 + 16*x - 10*x^2 + (-67 + a)*x^3 - 21*x^4)/(5*(2 - 2*x^2 + 3*x^3 + x^4)^(1/5))`
- `trager-deriv-alg-09` (unknown): UNEVAL -> CERT — `((-4 + a)*x - 19*x^2 + 8*x^3 + 15*x^4 + 12*x^5)/(3*(-2*x + x^2 + 3*x^3 + 3*x^4)^(4/3))`
- `trager-deriv-alg-10` (unknown): UNEVAL -> CERT — `(39 + a + 24*x - 21*x^2 + 6*x^3)/(6*(3 + 3*x - 3*x^2 + x^3)^(7/6))`
- `trager-deriv-alg-11` (unknown): UNEVAL -> CERT — `(-39 + a - 66*x + 27*x^2 + 36*x^3 - 3*x^4)/(5*(3 - 2*x - 2*x^2 + x^3)^(8/5))`
- `trager-deriv-alg-12` (unknown): UNEVAL -> CERT — `(-12 + a - 55*x + 2*x^2 + 78*x^3 + 42*x^4 + 13*x^5 + 48*x^6)/(3*(2 - x - 2*x^2 - 2*x^3 - x^5)^(1/3))`
- `trager-deriv-alg-13` (unknown): UNEVAL -> CERT — `(4 + (-24 + a)*x - 24*x^2 + 20*x^3)/(4*(-2*x + x^2)^(3/4))`
- `trager-deriv-alg-14` (unknown): UNEVAL -> CERT — `(17 + a + 36*x - 33*x^2 - 18*x^3)/(4*(-1 + 3*x + 3*x^2)^(7/4))`
- `trager-deriv-alg-15` (unknown): UNEVAL -> CERT — `(-18 + (-65 + a)*x - 70*x^2 + 12*x^3 + 10*x^4 - 31*x^5 + 102*x^6)/(6*(-2 - 3*x + x^2 - 2*x^4 + 3*x^5)^(5/6))`
- `trager-deriv-alg-18` (unknown): UNEVAL -> CERT — `(22 + a - 6*x - 32*x^2 - 30*x^3 - 10*x^4 - 12*x^5)/(4*(2 - 2*x - 2*x^2 - 2*x^3 - 3*x^4)^(5/4))`
- `trager-deriv-alg-19` (unknown): UNEVAL -> CERT — `(6 - 36*x - 36*x^2 + (24 + a)*x^3)/(3*(2 + x - 3*x^2 - 2*x^3 + x^4)^(5/3))`
- `trager-deriv-log-14` (unknown): UNEVAL -> ELEM — `(12 - 4*x - 2*x^3 - 12*x^4 + 6*x^5)/(2*Sqrt[3 - 2*x - 2*x^2 + x^3 + 3*x^4 - x^5]*(-12 + 8*x + (9 + a)*x^2 - 4*x^3 - 12*x^4 + 4*x^5))`
- `trager-deriv-log-22` (unknown): UNEVAL -> ELEM — `(4*x + 9*x^2 + 7*x^3 - 2*x^4 + 6*x^5 + 6*x^6)/(2*Sqrt[-1 - x - 2*x^2 + 2*x^3 + x^4]*(1 + 7*x + 17*x^2 + 19*x^3 + (6 + a)*x^4 - 24*x^5 - 9*x^6))`
- `trager-random-002` (unknown): UNEVAL -> CERT — `((3 + (3 + a)*x)*(2 + 3*x + x^2 + 3*x^3)^(2/3))/(-2 - x)`
- `trager-random-003` (unknown): UNEVAL -> CERT — `(2 + (3 + a)*x)/((3 - 2*x)*(3 - x + 2*x^2 + x^4)^(3/4))`
- `trager-random-005` (unknown): UNEVAL -> CERT — `((1 + a + x^2)*Sqrt[2 + x + 3*x^2 - 3*x^3 + x^4])/(2 + 2*x)`
- `trager-random-006` (unknown): UNEVAL -> CERT — `(2*((-2 + a)*x + x^2)^(2/3))/3`
- `trager-random-007` (unknown): UNEVAL -> CERT — `(-2 + a + 3*x + x^2)^(2/3)/2`
- `trager-random-013` (unknown): UNEVAL -> CERT — `((3 + x + x^2)^(1/4)*(2 + (-3 + a)*x + x^2))/3`
- `trager-random-014` (unknown): TIMEOUT -> CERT — `((3 + a + x + x^2)^(4/5)*(2 + 3*x + 3*x^2))/(1 - 2*x)`
- `trager-random-015` (unknown): UNEVAL -> CERT — `((x - 2*x^2)*Sqrt[(-1 + a)*x - x^3 + x^4])/2`
- `trager-random-019` (unknown): UNEVAL -> CERT — `(-2 + 3*x + (2 + a)*x^2)/(2*(1 + x - x^2 + x^3)^(1/6))`
- `trager-random-023` (unknown): UNEVAL -> CERT — `((-1 + (-2 + a)*x + x^2)*(2*x + 2*x^2 + x^3)^(1/4))/3`
- `trager-random-025` (unknown): UNEVAL -> CERT — `(-2 + (-1 + a)*x + 3*x^2)/((2 - x)*Sqrt[-3 + 2*x - 3*x^2 + x^3])`
- `trager-random-026` (unknown): UNEVAL -> CERT — `(2 + x - 3*x^2)/(2*(2 + a - 3*x - x^2 + x^3)^(1/5))`
- `trager-random-029` (unknown): UNEVAL -> CERT — `((3 - 2*x + 3*x^2)*(3 + a + 3*x + 2*x^2 + 2*x^3)^(2/5))/x`
- `trager-random-035` (unknown): UNEVAL -> CERT — `(-3 + 3*x + (-1 + a)*x^2)/Sqrt[3 - x^2 + x^3 + x^4]`
- `trager-random-041` (unknown): UNEVAL -> CERT — `(-3*(-3 - 2*x - 3*x^2 + x^3)^(3/4))/(-2*x + (3 + a)*x^2)`
- `trager-random-042` (unknown): UNEVAL -> CERT — `(3 + a + 2*x + x^2)/((2 - 2*x)*Sqrt[-3 + x + x^3])`
- `trager-random-043` (unknown): UNEVAL -> CERT — `(-2 - 3*x^2)/(2*(2 + a + 3*x + 2*x^2 + x^3)^(2/3))`
- `trager-random-044` (unknown): UNEVAL -> CERT — `Sqrt[3 + (2 + a)*x - 2*x^2 + 3*x^3 + x^4]`
- `trager-random-045` (unknown): UNEVAL -> CERT — `-(2 + (-1 + a)*x - x^2 - 2*x^3 - 3*x^4)^(3/5)`
- `trager-random-046` (unknown): UNEVAL -> CERT — `1/Sqrt[2 + (2 + a)*x + x^2 - 3*x^4]`
- `trager-random-050` (unknown): UNEVAL -> CERT — `((3 - 2*x)*(-2 + a - 3*x + 2*x^2)^(1/4))/2`
- `trager-random-051` (unknown): UNEVAL -> CERT — `((-1 + 2*x)*(-2 + 3*x - 2*x^2 + x^3)^(2/3))/(-1 + (2 + a)*x - 3*x^2)`
- `trager-random-052` (unknown): UNEVAL -> CERT — `(-2*(-3 + (-3 + a)*x + 3*x^2 + x^3)^(1/4))/3`
- `trager-random-054` (unknown): UNEVAL -> CERT — `((-1 - 2*x + 3*x^2)*(-2 + a - 2*x - 3*x^2 - x^3 + x^4)^(1/6))/3`
- `trager-random-055` (unknown): UNEVAL -> CERT — `((1 + x^2)*(-2 + 3*x + x^3)^(5/6))/(-1 - 2*x + (-3 + a)*x^2)`
- `trager-random-056` (unknown): UNEVAL -> CERT — `(-3 - x)/(3*(3 + a + 2*x + x^2 + 2*x^3 + 2*x^4)^(2/3))`
- `trager-random-057` (unknown): UNEVAL -> CERT — `((2 + 2*x)*(x + x^2 + x^3)^(2/3))/(-1 + (-3 + a)*x)`
- `trager-random-058` (unknown): UNEVAL -> CERT — `(-3 - 3*x)/(2*Sqrt[2 + a - 3*x + 3*x^2 + x^3 + x^4])`
- `trager-random-059` (unknown): UNEVAL -> CERT — `-3/((2 + (3 + a)*x)*Sqrt[-3 + 2*x + x^2 + x^3])`
- `agen-2-5` (elem): UNEVAL -> ELEM — `(4 - 2*a + 10*x + a*x - 18*x^2 + 6*x^3 + 2*(-1 - 2*x + a*x - 3*x^2 + 2*x^3))/(2*Sqrt[-1 - 2*x + a*x - 3*x^2 + 2*x^3])`
- `agen-2-10` (elem): UNEVAL -> ELEM — `(3*(-3 + x - 2*x^2 + 3*x^3)^(1/4) - 11*x*(-3 + x - 2*x^2 + 3*x^3)^(1/4) + 23*x^2*(-3 + x - 2*x^2 + 3*x^3)^(1/4) + 9*a*x^2*(-3 + x - 2*x^2 + 3*x^3)^(1/4) + 9*x^3`
- `bbottomp-2-11` (unknown): TIMEOUT -> CERT — `E^x*Sqrt[x] - E^x/(2*(3 + x)^(3/2)) + 2*E^x*Sqrt[3 + x] + (2*E^x)/Sqrt[3 + a + x]`
- `bbottomp-2-12` (unknown): TIMEOUT -> CERT — `E^x*Sqrt[x] - E^x/(2*(1 + x)^(3/2)) + 2*E^x*Sqrt[1 + x] + (2*E^x)/Sqrt[1 + a + x]`
- `ratb-2-8` (elem): TIMEOUT -> ELEM — `(-2*E^x^2*(-1 + a + 2*E^x^2)^(1/3)*x)/(1 + 2*E^x^2 + E^(2*x^2)) + (4*E^x^2*x)/(3*(-1 + a + 2*E^x^2)^(2/3) + 3*E^x^2*(-1 + a + 2*E^x^2)^(2/3))`
- `agen-3-3` (elem): UNEVAL -> ELEM — `(-2 - a + 2*x - 6*x^2 - 2*Sqrt[(2 + a)*x - x^2 + 2*x^3])/((2 + a)*x - x^2 + 2*x^3 + Sqrt[(2 + a)*x - x^2 + 2*x^3] + x*Sqrt[(2 + a)*x - x^2 + 2*x^3])`
- `agen-3-6` (elem): UNEVAL -> ELEM — `(3 + 9*x^2 + 6*x^4 + 4*x*(-2 + a + 3*x + 2*x^3))/Sqrt[-2 + a + 3*x + 2*x^3]`
- `agen-3-8` (elem): TIMEOUT -> ELEM — `(6*(-1 + 3*x + 3*x^2)^(1/4) + 15*x*(-1 + 3*x + 3*x^2)^(1/4) + 4*a*x*(-1 + 3*x + 3*x^2)^(1/4) + 6*x^2*(-1 + 3*x + 3*x^2)^(1/4) + 2*a*x^2*(-1 + 3*x + 3*x^2)^(1/4)`
- `agen-3-9` (elem): TIMEOUT -> ELEM — `(4*x + 4*x^2 - 3*(1 + a + 2*x^2))/(3*(1 + a + 2*x^2)^(2/3) + 6*x*(1 + a + 2*x^2)^(2/3) + 3*x^2*(1 + a + 2*x^2)^(2/3))`
- `ratb-3-8` (elem): TIMEOUT -> ELEM — `(2 + 6*E^x)/(3*(3 + 3*E^x + x)^(1/3)) - (E^x*(3 + 3*E^x + a*E^x + x)^(1/3))/(4 + 4*E^x + E^(2*x)) + (1 + 3*E^x + a*E^x)/(6*(3 + 3*E^x + a*E^x + x)^(2/3) + 3*E^x`
- `constk-3-3` (elem): TIMEOUT -> ELEM — `(-E^(1 + x) + a*E^(1 + x) - Sqrt[1 + 2*E^(1 + x)] + a*Sqrt[1 + 2*E^(1 + x)] + 2*(1 + 2*E^(1 + x)) + (1 + 2*E^(1 + x))^(3/2) + 2*E^(1 + x)*x + E^(1 + x)*Sqrt[1 +`
- `cA-1-290` (elem): smaller 514 -> 96 — `(-1 - a + a^2 - x + 2*a*x + x^2 + (1 + x)^(1/3)*(1 + a + x)^(2/3))/(-(1 + x)^(1/3) - a*(1 + x)^(1/3) + a^2*(1 + x)^(1/3) - x*(1 + x)^(1/3) + 2*a*x*(1 + x)^(1/3)`
- `cA-1-298` (elem): smaller 396 -> 67 — `(-2*(3 + x)^(1/4) + a^2*(3 + x)^(1/4) + 2*a*x*(3 + x)^(1/4) + x^2*(3 + x)^(1/4) + 4*(3 + a + x)^(1/4))/(-2 + a^2 + 2*a*x + x^2)`
- `cB-1-426` (unknown): TIMEOUT -> CERT — `(2*(1 + x + x^2)^(1/4))/(1 + a + 2*x)`
- `cB-1-1188` (unknown): TIMEOUT -> CERT — `(2*(1 + x + x^2)^(2/3))/(1 + a + 2*x)`
- `cA-2-6` (elem): smaller 326 -> 96 — `(32*Sqrt[1 + x + x^2] + 64*x*Sqrt[1 + x + x^2] - 21*Sqrt[1 + a + x + x^2] - 45*x*Sqrt[1 + a + x + x^2] + 27*x^2*Sqrt[1 + a + x + x^2] + 16*x^3*Sqrt[1 + a + x + `
- `cA-2-30` (elem): smaller 456 -> 126 — `(-2560*Sqrt[1 + x + x^2] - 5120*x*Sqrt[1 + x + x^2] + 2501*Sqrt[1 + a + x + x^2] + 4829*x*Sqrt[1 + a + x + x^2] + 117*x^2*Sqrt[1 + a + x + x^2] + 688*x^3*Sqrt[1`
- `cA-2-100` (elem): TIMEOUT -> ELEM — `(-2*(3 + x)^(1/6) + a^2*(3 + x)^(1/6) + 2*a*x*(3 + x)^(1/6) + x^2*(3 + x)^(1/6) + 4*(3 + a + x)^(1/6))/(-2 + a^2 + 2*a*x + x^2)`
- `cA-2-112` (elem): smaller 258 -> 71 — `(-24*Sqrt[1 + x] - 18*x*Sqrt[1 + x] + 15*x^2*Sqrt[1 + x] + 20*Sqrt[1 + x]*(1 + a + x) - 10*x*Sqrt[1 + x]*(1 + a + x) + 4*(1 + a + x)^(5/2))/(-18*Sqrt[1 + x]*(1 `
- `cA-2-114` (elem): smaller 396 -> 69 — `(-3*(3 + x)^(1/4) + a^2*(3 + x)^(1/4) + 2*a*x*(3 + x)^(1/4) + x^2*(3 + x)^(1/4) + (3 + a + x)^(1/4))/(-3 + a^2 + 2*a*x + x^2)`
- `cA-2-115` (elem): smaller 360 -> 96 — `(-42*Sqrt[2 + x + x^2] - 84*x*Sqrt[2 + x + x^2] + 88*Sqrt[2 + a + x + x^2] + 123*x*Sqrt[2 + a + x + x^2] + 47*x^2*Sqrt[2 + a + x + x^2] + 16*x^3*Sqrt[2 + a + x `
- `cA-2-118` (elem): smaller 408 -> 111 — `(340*Sqrt[2 + x + x^2] + 680*x*Sqrt[2 + x + x^2] - 442*Sqrt[2 + a + x + x^2] - 491*x*Sqrt[2 + a + x + x^2] + 181*x^2*Sqrt[2 + a + x + x^2] + 304*x^3*Sqrt[2 + a `
- `cB-2-371` (unknown): TIMEOUT -> CERT — `2/((1 + a + 2*x)*(2 + x + x^2)^(1/4))`
- `cB-2-461` (unknown): TIMEOUT -> CERT — `(2*(2 + x + x^2)^(2/3))/(1 + a + 2*x)`
- `cB-3-193` (unknown): TIMEOUT -> CERT — `(2*(2 + x + x^2)^(1/4))/(1 + a + 2*x)`
- `cA-3-209` (elem): UNEVAL -> ELEM — `(1 + x^4 + x^4*(1 + x^4) + x^4*(1 + x^4)^(3/2)*Sqrt[1 + a + x^4])/(x*(1 + x^4)^(3/2))`
- `cA-3-258` (elem): smaller 548 -> 103 — `(2*Sqrt[2 + x] + a*Sqrt[2 + x] - x*Sqrt[2 + x] - a*x*Sqrt[2 + x] - x^2*Sqrt[2 + x] + 7*(2 + a + x)^(3/2) + 7*a*(2 + a + x)^(3/2) + 9*x*(2 + a + x)^(3/2) + 2*a*x`
- `cA-3-284` (elem): smaller 410 -> 108 — `(-3*(1 + x)^(2/3) + a^2*(1 + x)^(2/3) + 2*a*x*(1 + x)^(2/3) + x^2*(1 + x)^(2/3) + (1 + a + x)^(2/3))/(-3 + a^2 + 2*a*x + x^2)`
- `cB-3-1140` (elem): smaller 114 -> 51 — `((4*x)/(3 + x^2)^(2/3) + (8*x)/(-2*(3 + a + x^2)^(2/3) + (3 + a + x^2)^(5/3)))/2`
- `cA-4-94` (elem): TIMEOUT -> ELEM — `(-3*(3 + x)^(1/6) + a^2*(3 + x)^(1/6) + 2*a*x*(3 + x)^(1/6) + x^2*(3 + x)^(1/6) + 5*(3 + a + x)^(1/6))/(-3 + a^2 + 2*a*x + x^2)`
- `cA-4-170` (unknown): UNEVAL -> CERT — `x/((-2 + (2 + a)*x^2)*(1 + x^6)^(1/3))`
- `cA-4-244` (elem): smaller 558 -> 113 — `(-2*Sqrt[2 + x] - a*Sqrt[2 + x] + x*Sqrt[2 + x] + a*x*Sqrt[2 + x] + x^2*Sqrt[2 + x] - 7*(2 + a + x)^(3/2) - 7*a*(2 + a + x)^(3/2) - 9*x*(2 + a + x)^(3/2) - 2*a*`
- `cM-1-36` (unknown): TIMEOUT -> CERT — `(2*E^(2*x))/((1 + a + E^(2*x))*(-1 + E^(6*x))^(1/3))`
- `cM-1-142` (elem): TIMEOUT -> ELEM — `-4/(x*Log[x]) + (-4 - 8*Log[x] + 4*Sqrt[1 + Log[x] + Log[x]^2])/(2*x*Sqrt[1 + Log[x] + Log[x]^2] + x*Log[x]*Sqrt[1 + Log[x] + Log[x]^2] - 2*x*(1 + Log[x] + Log[`
- `cM-2-79` (unknown): smaller 336 -> 115 — `1/(x*Log[x]*(1 + a + Log[x]^2)^(2/3))`
- `cM-2-136` (elem): TIMEOUT -> ELEM — `(45*Log[1 + x]*Sqrt[3 + Log[1 + x]^4] + 15*a*Log[1 + x]*Sqrt[3 + Log[1 + x]^4] + 30*Log[1 + x]^5*Sqrt[3 + Log[1 + x]^4] + 27*Log[1 + x]*Sqrt[3 + a + Log[1 + x]^`
- `cM-2-222` (elem): TIMEOUT -> ELEM — `(5*Log[1 + x]^2*Sqrt[1 + Log[1 + x]^6] + 5*a*Log[1 + x]^2*Sqrt[1 + Log[1 + x]^6] + 10*Log[1 + x]^8*Sqrt[1 + Log[1 + x]^6] + 3*Log[1 + x]^2*Sqrt[1 + a + Log[1 + `
- `cM-3-104` (unknown): TIMEOUT -> ELEM — `(2*E^x^2*x)/(Sqrt[3 + a + E^x^2]*(3 + E^(2*x^2)))`
- `cM-3-240` (unknown): TIMEOUT -> ELEM — `(E^x^2*(2 + E^x^2)*x)/(Sqrt[3 + a + E^x^2]*(1 + E^(2*x^2)))`
- `cM-3-322` (elem): TIMEOUT -> ELEM — `(1 + Log[x]^6 - Log[x]^12 - Log[x]^18 - Sqrt[1 + Log[x]^6] + a*Sqrt[1 + Log[x]^6] - 3*Log[x]^6*Sqrt[1 + Log[x]^6] + 2*a*Log[x]^6*Sqrt[1 + Log[x]^6] - 2*Log[x]^1`
- `cM-4-295` (unknown): smaller 337 -> 115 — `1/(x*Log[x]*(2 + a + Log[x]^3)^(2/3))`
