# 258 — parametric A/B of the Trager and Bronstein cascade stages

153 cases; A = shipped (`Integrate`$CascadeParameters = False`), B = `True`; timeout 30 s, check at a = 1/7.

- hard failures in B: **1**
- B worse than A: **110**
- B better than A: **27**

| A | B | cases |
|---|---|---|
| ELEM | ELEM | 117 |
| TIMEOUT | ELEM | 24 |
| UNEVAL | UNEVAL | 5 |
| UNEVAL | ELEM | 3 |
| ELEM | TIMEOUT | 2 |
| WRONG | WRONG | 1 |
| TIMEOUT | TIMEOUT | 1 |

## Hard failures

- `trager-cheb-002` (unknown): WRONG — `x^2*Sqrt[-1 + a + 2*x^4]`

## Worse

- `logxy-1-3` (elem): larger 78 -> 2397 — `(-2 - a + 6*E^x + 3*a*E^x + 6*(-1 + 3*E^x - x)^(2/3) + 3*a*(-1 + 3*E^x - x)^(2/3) + 12*(-1 + 3*E^x - x) + 6*(-1 + 3*E^x - x)^(4/3) - 4*x + 12*E^x*x - 2*(-1 + 3*`
- `logxy-1-8` (elem): larger 97 -> 177 — `Sqrt[-2 - 2*x + 3*Log[1 + E^x]] + (-2*x + E^x*x)/(2*Sqrt[-2 - 2*x + 3*Log[1 + E^x]] + 2*E^x*Sqrt[-2 - 2*x + 3*Log[1 + E^x]]) + (-2 + E^x + 2*Sqrt[-2 - 2*x + 3*L`
- `logxy-1-10` (elem): larger 73 -> 2750 — `(-2 + E^x + 3*(3 + E^x - 2*x)^(2/3) + 3*(3 + E^x - 2*x) + 3*a*(3 + E^x - 2*x) + 3*(3 + E^x - 2*x)^(4/3) - 2*x - 2*a*x + E^x*x + a*E^x*x - 2*(3 + E^x - 2*x)^(1/3`
- `ratb-1-3` (elem): larger 220 -> 554 — `(3*E^x^2*x)/(2*(-2 + E^x^2)^(1/4)) - (2*E^x^2*(-2 + a + E^x^2)^(1/4)*x)/(4 + 4*E^x^2 + E^(2*x^2)) + (E^x^2*x)/(4*(-2 + a + E^x^2)^(3/4) + 2*E^x^2*(-2 + a + E^x^`
- `constk-1-3` (elem): larger 82 -> 2627 — `(-4 - 2*a + 2*E^(1 + x) + a*E^(1 + x) + 6*(1 + E^(1 + x) - 2*x)^(2/3) + 3*a*(1 + E^(1 + x) - 2*x)^(2/3) + 6*(1 + E^(1 + x) - 2*x) + 3*(1 + E^(1 + x) - 2*x)^(4/3`
- `tower2-1-5` (elem): larger 15 -> 210 — `(1 + E^x*x + 3*x*(1 + a + E^x + Log[x])^(2/3))/(3*x^2*(1 + a + E^x + Log[x])^(2/3) + 3*x*(1 + a + E^x + Log[x]))`
- `tower2-1-6` (elem): slower 1.0s -> 21.3s — `(2 + E^x*x)/(4*x*(E^x + 2*Log[x])^(3/4)) + (2 + E^x*x)/(4*x*(E^x + 2*Log[x])^(3/4) + 4*a*x*(E^x + 2*Log[x])^(3/4) + 4*x*(E^x + 2*Log[x]))`
- `tower2-1-8` (elem): larger 15 -> 1990 — `(1 + E^x*x + 4*x*(-1 + a + E^x + Log[x])^(3/4))/(4*x^2*(-1 + a + E^x + Log[x])^(3/4) + 4*x*(-1 + a + E^x + Log[x]))`
- `tower2-1-12` (elem): larger 17 -> 236 — `(2 + E^x*x + 3*x*(-1 + a + E^x + 2*Log[x])^(2/3))/(3*x^2*(-1 + a + E^x + 2*Log[x])^(2/3) + 3*x*(-1 + a + E^x + 2*Log[x]))`
- `expcube-1-2` (elem): larger 17 -> 999 — `(-4*E^x^2*x)/(6*(3 + E^x^2)^(2/3) + 3*a*(3 + E^x^2)^(2/3) + 3*(3 + E^x^2))`
- `expcube-1-5` (elem): larger 15 -> 990 — `(2*E^x^2*x)/(6*(1 + E^x^2)^(2/3) + 3*a*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2))`
- `expcube-1-7` (elem): larger 15 -> 871 — `(2*E^x)/(6*(3 + E^x)^(2/3) + 3*a*(3 + E^x)^(2/3) + 3*(3 + E^x))`
- `expcube-1-8` (elem): larger 15 -> 466 — `(2*E^(2*x))/(-3*(3 + a + E^(2*x))^(2/3) + 3*(3 + a + E^(2*x)))`
- `expcube-1-11` (elem): larger 27 -> 1003 — `(2*E^x^2*x)/(3*(2 + E^x^2)^(2/3)) + (2*E^x^2*x)/(-3*(2 + E^x^2)^(2/3) + 3*a*(2 + E^x^2)^(2/3) + 3*(2 + E^x^2))`
- `trager-ratrad-019` (unknown): slower 2.3s -> 12.5s — `1/((2 + a + x^2)*(1 + x^3)^(1/4))`
- `trager-evensub-005` (unknown): larger 107 -> 342 — `x/(Sqrt[1 + a + x^4]*(1 + x^2 + x^4))`
- `trager-torsion-006` (unknown): slower 2.2s -> 12.4s — `1/((1 + a + x^2)*(1 + x^3)^(1/4))`
- `trager-torsion-007` (unknown): slower 2.2s -> 12.4s — `1/((3 + a + x^2)*(1 + x^3)^(1/4))`
- `trager-torsion-008` (unknown): slower 3.2s -> 13.4s — `x/((2 + a + x^2)*(1 + x^3)^(3/4))`
- `trager-deriv-log-25` (unknown): slower 4.7s -> 14.8s — `(-3 + 120*x - 71*x^2 - 8*x^3 + 45*x^4 + (-16 + a)*x^5 + 5*x^6 - 24*x^8)/(2*Sqrt[3 - 3*x + x^2 + 2*x^3 - x^4 - x^5 + x^6]*(-26 - 9*x + 21*x^2 - 18*x^3 - 10*x^4 +`
- `logxy-2-2` (elem): larger 65 -> 2842 — `(E^x + 2*(-3 + a + 2*E^x)^(3/4))/(2*(-3 + a + 2*E^x)^(3/4) + 2*(-3 + a + 2*E^x) + 2*(-3 + a + 2*E^x)^(3/4)*x)`
- `logxy-2-6` (elem): larger 32 -> 5782, slower 1.1s -> 9.1s — `(6 + 3*a + 12*x + 6*x^2 + 6*x*(-3 + 3*Log[x])^(1/4) + 8*x*(-3 + 3*Log[x])^(3/4) + 4*a*x*(-3 + 3*Log[x])^(3/4) + 16*x*(-3 + 3*Log[x]) + 8*x^2*(-3 + 3*Log[x]) + 8`
- `logxy-2-10` (elem): larger 77 -> 1782 — `(2 - 2*E^x - 3*(2 + a + 2*E^x - 2*x)^(2/3))/(3*(2 + a + 2*E^x - 2*x)^(2/3) + 3*(2 + a + 2*E^x - 2*x) + 3*(2 + a + 2*E^x - 2*x)^(2/3)*x)`
- `ratb-2-6` (elem): larger 32 -> 726 — `4/(3*x*(3 + 2*Log[x])^(1/3)) - (3 + 2*Log[x] + a*Log[x])^(1/3)/(x + 2*x*Log[x] + x*Log[x]^2) + (2 + a)/(3*x*(3 + 2*Log[x] + a*Log[x])^(2/3) + 3*x*Log[x]*(3 + 2*`
- `constk-2-3` (elem): slower 5.2s -> 15.6s — `(-2 + 2*E^x + x)^(1/4) + (x + 2*E^x*x)/(4*(-2 + 2*E^x + x)^(3/4)) + (1 + 2*E^x + 4*(-2 + 2*E^x + x)^(3/4))/(8*(-2 + 2*E^x + x)^(3/4) + 4*a*(-2 + 2*E^x + x)^(3/4`
- `constk-2-4` (elem): larger 78 -> 2477 — `(-2 + a - 3*x + 2*a*x + 3*x^2 + 2*x^3 + x*(-3 + 2*x + Log[3*x])^(1/3) + 2*x^2*(-3 + 2*x + Log[3*x])^(1/3) - 6*x*(-3 + 2*x + Log[3*x])^(2/3) + 3*a*x*(-3 + 2*x + `
- `constk-2-5` (elem): larger 83 -> 2382 — `(2 - 2*a - 6*E^x + 6*a*E^x - 3*(2 + 6*E^x - 2*x)^(2/3) + 3*a*(2 + 6*E^x - 2*x)^(2/3) + 6*(2 + 6*E^x - 2*x) + 3*(2 + 6*E^x - 2*x)^(4/3) - 4*x + 12*E^x*x - 2*(2 +`
- `constk-2-6` (elem): larger 93 -> 2814 — `(-7 + 7*a + 17*x - 3*a*x + x^2 - 3*x^3 + 7*x*(2 - x + 3*Log[2 + 3*x])^(1/3) - 3*x^2*(2 - x + 3*Log[2 + 3*x])^(1/3) - 6*(2 - x + 3*Log[2 + 3*x])^(2/3) + 6*a*(2 -`
- `constk-2-8` (elem): slower 1.7s -> 11.9s — `(10 + 5*a + 24*x + 2*a*x + 18*x^2 + 4*x^3 + 10*x*(-1 + 2*x + 3*Log[3 + 3*x])^(1/4) + 4*x^2*(-1 + 2*x + 3*Log[3 + 3*x])^(1/4) + 8*(-1 + 2*x + 3*Log[3 + 3*x])^(3/`
- `tower2-2-1` (elem): larger 67 -> 1319 — `(1 + E^x*x)/(3*x*(1 + E^x + Log[x])^(2/3)) + (1 + E^x*x)/(3*x*(1 + E^x + Log[x])^(2/3) + 3*a*x*(1 + E^x + Log[x])^(2/3) + 3*x*(1 + E^x + Log[x]))`
- `tower2-2-6` (elem): larger 17 -> 236 — `(2 + E^x*x + 3*x*(-1 + a + E^x + 2*Log[x])^(2/3))/(3*x^2*(-1 + a + E^x + 2*Log[x])^(2/3) + 3*x*(-1 + a + E^x + 2*Log[x]))`
- `tower2-2-7` (elem): larger 17 -> 2052 — `(1 + 2*E^x*x + 4*x*(2 + a + 2*E^x + Log[x])^(3/4))/(4*x^2*(2 + a + 2*E^x + Log[x])^(3/4) + 4*x*(2 + a + 2*E^x + Log[x]))`
- `lccontent-2-5` (elem): larger 14 -> 7434 — `(1 + Log[x])/(3*(2 + x*Log[x])^(2/3) + 3*a*(2 + x*Log[x])^(2/3) + 3*(2 + x*Log[x]))`
- `lccontent-2-9` (elem): larger 22 -> 8763 — `(4*x*Log[1 + x^2] + 4*x^2*Log[1 + x^2] + Log[1 + x^2]^2 + x^2*Log[1 + x^2]^2)/(3*(1 + a + (1 + x)*Log[1 + x^2]^2)^(2/3) + 3*x^2*(1 + a + (1 + x)*Log[1 + x^2]^2)`
- `expcube-2-2` (elem): larger 15 -> 992 — `(2*E^x^2*x)/(3*(2 + E^x^2)^(2/3) + 3*a*(2 + E^x^2)^(2/3) + 3*(2 + E^x^2))`
- `expcube-2-3` (elem): larger 15 -> 511 — `(2*E^(2*x))/(6*(3 + a + E^(2*x))^(2/3) + 3*(3 + a + E^(2*x)))`
- `expcube-2-5` (elem): larger 17 -> 999 — `(-4*E^(2*x))/(3*(3 + E^(2*x))^(2/3) + 3*a*(3 + E^(2*x))^(2/3) + 3*(3 + E^(2*x)))`
- `expcube-2-7` (elem): larger 15 -> 864 — `-(E^x/(3*(2 + E^x)^(2/3) + 3*a*(2 + E^x)^(2/3) + 3*(2 + E^x)))`
- `expcube-2-8` (elem): larger 15 -> 871 — `(2*E^x)/(3*(3 + E^x)^(2/3) + 3*a*(3 + E^x)^(2/3) + 3*(3 + E^x))`
- `expcube-2-9` (elem): larger 17 -> 992 — `(-2*E^x^2*x)/(3*(2 + E^x^2)^(2/3) + 3*a*(2 + E^x^2)^(2/3) + 3*(2 + E^x^2))`
- `expcube-2-10` (elem): larger 27 -> 1003 — `(2*E^(2*x))/(3*(3 + E^(2*x))^(2/3)) + (2*E^(2*x))/(6*(3 + E^(2*x))^(2/3) + 3*a*(3 + E^(2*x))^(2/3) + 3*(3 + E^(2*x)))`
- `expcube-2-11` (elem): larger 17 -> 990 — `(-2*E^(2*x))/(6*(1 + E^(2*x))^(2/3) + 3*a*(1 + E^(2*x))^(2/3) + 3*(1 + E^(2*x)))`
- `logxy-3-4` (elem): larger 79 -> 1916 — `(-4 - 4*E^x - 6*(-3 + a + 2*E^x + 2*x)^(2/3))/(6*(-3 + a + 2*E^x + 2*x)^(2/3) + 3*x*(-3 + a + 2*E^x + 2*x)^(2/3) + 3*(-3 + a + 2*E^x + 2*x))`
- `logxy-3-5` (elem): larger 38 -> 2248 — `(4 + 2*a + 12*x + 2*a*x + 12*x^2 + 4*x^3 + 4*x*(1 + 2*x + 2*Log[x])^(1/3) + 4*x^2*(1 + 2*x + 2*Log[x])^(1/3) + 6*x*(1 + 2*x + 2*Log[x])^(2/3) + 3*a*x*(1 + 2*x +`
- `logxy-3-6` (elem): larger 87 -> 1560 — `(-2 - 4*E^x^2*x - 3*(-2 + a + 2*E^x^2 + 2*x)^(2/3))/(6*(-2 + a + 2*E^x^2 + 2*x)^(2/3) + 3*x*(-2 + a + 2*E^x^2 + 2*x)^(2/3) + 3*(-2 + a + 2*E^x^2 + 2*x))`
- `logxy-3-7` (elem): larger 76 -> 170 — `(9 + 11*x + 9*a*x + 11*x^2 + 2*a*x^2 + 2*x^3 + 6*Sqrt[-2 + x + 3*Log[3 + 2*x]] + 13*x*Sqrt[-2 + x + 3*Log[3 + 2*x]] + 2*x^2*Sqrt[-2 + x + 3*Log[3 + 2*x]] + 6*(-`
- `logxy-3-9` (elem): larger 19 -> 1506 — `(1 + 6*E^x^2*x + 3*(3 + a + 3*E^x^2 + x)^(2/3))/(6*(3 + a + 3*E^x^2 + x)^(2/3) + 3*x*(3 + a + 3*E^x^2 + x)^(2/3) + 3*(3 + a + 3*E^x^2 + x))`
- `constk-3-5` (elem): larger 98 -> 3105 — `(4 + 2*a + 8*x + 2*a*x + 10*x^2 + 2*a*x^2 + 6*x^3 + 2*x^4 + 2*x*(2 + 2*x + Log[2*(1 + x^2)])^(1/3) + 2*x^2*(2 + 2*x + Log[2*(1 + x^2)])^(1/3) + 2*x^3*(2 + 2*x +`
- `constk-3-9` (elem): larger 86 -> 2544 — `(-4 + 2*a - 2*x + 2*a*x + 4*x^2 + 2*x^3 + 2*x*(-3 + 2*x + 2*Log[3*x])^(1/3) + 2*x^2*(-3 + 2*x + 2*Log[3*x])^(1/3) - 6*x*(-3 + 2*x + 2*Log[3*x])^(2/3) + 3*a*x*(-`
- `tower2-3-2` (elem): larger 87 -> 2980, slower 2.5s -> 19.6s — `(1 + E^x*x)/(2*x*(1 + 2*E^x + 2*Log[x])^(3/4)) + (1 + E^x*x)/(2*x + 4*E^x*x + 4*x*Log[x] + 2*x*(1 + 2*E^x + 2*Log[x])^(3/4) + 2*a*x*(1 + 2*E^x + 2*Log[x])^(3/4)`
- `tower2-3-3` (elem): larger 15 -> 1990 — `(1 + E^x*x + 4*x*(-1 + a + E^x + Log[x])^(3/4))/(4*x^2*(-1 + a + E^x + Log[x])^(3/4) + 4*x*(-1 + a + E^x + Log[x]))`
- `tower2-3-6` (elem): larger 87 -> 1441 — `(2 + 2*E^x*x)/(3*x*(1 + 2*E^x + 2*Log[x])^(2/3)) + (2 + 2*E^x*x)/(3*x*(1 + 2*E^x + 2*Log[x])^(2/3) + 3*a*x*(1 + 2*E^x + 2*Log[x])^(2/3) + 3*x*(1 + 2*E^x + 2*Log`
- `tower2-3-7` (elem): larger 77 -> 1380 — `(1 + 2*E^x*x)/(3*x*(1 + 2*E^x + Log[x])^(2/3)) + (1 + 2*E^x*x)/(3*x*(1 + 2*E^x + Log[x])^(2/3) + 3*a*x*(1 + 2*E^x + Log[x])^(2/3) + 3*x*(1 + 2*E^x + Log[x]))`
- `tower2-3-8` (elem): larger 67 -> 1210 — `(1 + E^x*x)/(3*x*(-1 + E^x + Log[x])^(2/3)) + (1 + E^x*x)/(3*x*(-1 + E^x + Log[x])^(2/3) + 3*a*x*(-1 + E^x + Log[x])^(2/3) + 3*x*(-1 + E^x + Log[x]))`
- `lccontent-3-4` (elem): larger 18 -> 6047 — `(2*Log[x] + 2*x*Log[x] + x*Log[x]^2)/(3*x*(1 + a + (1 + x)*Log[x]^2)^(2/3) + 3*x*(1 + a + (1 + x)*Log[x]^2))`
- `expcube-3-1` (elem): larger 13 -> 864 — `E^x/(6*(3 + E^x)^(2/3) + 3*a*(3 + E^x)^(2/3) + 3*(3 + E^x))`
- `expcube-3-2` (elem): larger 15 -> 992 — `(2*E^x^2*x)/(-3*(2 + E^x^2)^(2/3) + 3*a*(2 + E^x^2)^(2/3) + 3*(2 + E^x^2))`
- `expcube-3-5` (elem): larger 17 -> 992 — `(-2*E^(2*x))/(-3*(2 + E^(2*x))^(2/3) + 3*a*(2 + E^(2*x))^(2/3) + 3*(2 + E^(2*x)))`
- `expcube-3-7` (elem): larger 15 -> 871 — `(2*E^x)/(6*(2 + E^x)^(2/3) + 3*a*(2 + E^x)^(2/3) + 3*(2 + E^x))`
- `expcube-3-8` (elem): larger 17 -> 920 — `(-4*E^x^2*x)/(-3*(1 + E^x^2)^(2/3) + 3*a*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2))`
- `expcube-3-10` (elem): larger 15 -> 992 — `(2*E^x^2*x)/(-3*(3 + E^x^2)^(2/3) + 3*a*(3 + E^x^2)^(2/3) + 3*(3 + E^x^2))`
- `expcube-3-11` (elem): larger 17 -> 999 — `(4*E^x^2*x)/(3*(3 + E^x^2)^(2/3) + 3*a*(3 + E^x^2)^(2/3) + 3*(3 + E^x^2))`
- `cA-1-113` (elem): larger 90 -> 218 — `5/(2*Sqrt[1 + x^2]) + ((2*x)/(1 + x^2)^(3/2) + (2*x^2)/(1 + x^2)^(3/2) - 2/Sqrt[1 + x^2])/(1 - 1/(2*(1 + x^2)) - x/(1 + x^2) - x^2/(2*(1 + x^2))) + (4*x*Sqrt[1 `
- `cA-1-132` (elem): larger 129 -> 246 — `4/Sqrt[1 + x^2] + ((2*x)/(1 + x^2)^(3/2) + (2*x^2)/(1 + x^2)^(3/2) - 2/Sqrt[1 + x^2])/(1 - 1/(2*(1 + x^2)) - x/(1 + x^2) - x^2/(2*(1 + x^2))) + (-9 + 9*x^2 + 9*`
- `cA-1-263` (elem): slower 0.1s -> 10.1s — `(2*a*x^2 - 2*a*x^2*(1 + x^3)^(1/6) - 2*a*x^2*(1 + x^3)^(1/3) - a^2*x^2*(1 + x^3)^(1/3) + 4*a*x^2*Sqrt[1 + x^3] + a^2*x^2*Sqrt[1 + x^3] - 2*a*x^2*(1 + x^3)^(2/3)`
- `cA-3-147` (elem): larger 90 -> 217 — `3/(2*Sqrt[1 + x^2]) + (x/(1 + x^2)^(3/2) + x^2/(1 + x^2)^(3/2) - 1/Sqrt[1 + x^2])/(1 - 1/(2*(1 + x^2)) - x/(1 + x^2) - x^2/(2*(1 + x^2))) + (4*x*Sqrt[1 + x^2] +`
- `cM-1-98` (unknown): larger 220 -> 760 — `(E^x*(-1 + Log[1 + E^x])*(2 + a + Log[1 + E^x])^(3/2))/((1 + E^x)*(1 + Log[1 + E^x] + Log[1 + E^x]^2))`
- `cM-1-108` (unknown): slower 0.2s -> 10.2s — `1/(1 + a + E^(6*x))^(1/4)`
- `cM-1-111` (elem): larger 52 -> 273 — `(1 + a + Log[x])^(1/3)/(3*x - 2*x*(1 + a + Log[x]) + x*(1 + a + Log[x])^2)`
- `cM-1-176` (unknown): slower 0.1s -> 10.2s — `2/(1 + a + E^(12*x))^(1/4)`
- `cM-1-230` (unknown): larger 83 -> 367 — `(E^x*(1 + E^x)^(1/4))/(-2 + a + E^x)`
- `cM-1-259` (unknown): larger 540 -> 2773 — `-((E^x*(2 + E^x))/((1 + E^x + E^(2*x))^(3/2)*(-1 + E^x + (-1 + a)*E^(2*x))))`
- `cM-1-327` (unknown): larger 72 -> 474 — `(2*x*Sqrt[3 + a + Log[1 + x^2]])/((1 + x^2)*(2 + Log[1 + x^2]^2))`
- `cM-1-341` (unknown): larger 773 -> 2293 — `(-2*E^(2*x)*(1 + E^(2*x)))/((1 + (2 + a)*E^(2*x) - E^(4*x))*(3 + E^(4*x))^(3/2))`
- `cM-1-362` (unknown): larger 38 -> 394 — `(2*x)/((1 + x^2)*(2 + Log[1 + x^2])^(1/4)*(-1 + a + Log[1 + x^2]))`
- `cM-2-67` (unknown): larger 78 -> 689 — `-((E^x*(-2 + E^x))/(Sqrt[2 + a + E^x]*(-1 + E^x - 2*E^(2*x))))`
- `cM-2-97` (unknown): slower 1.6s -> 11.6s — `((1 + E^(-x))^(1/4)*(1 + (-2 + a)*E^x))/(-2 - 2*E^x + E^(2*x))`
- `cM-2-105` (unknown): larger 182 -> 1545 — `(2*E^(2*x)*(1 + E^(2*x))*(2 + a + E^(2*x))^(1/3))/(2 + E^(2*x) + E^(4*x))`
- `cM-2-130` (unknown): larger 119 -> 327 — `1/(x*Log[x]*(3 + a + Log[x])^(1/4))`
- `cM-2-155` (elem): slower 6.6s -> 27.0s — `(Log[x]^4 + Log[x]^2*(1 + Log[x]^6)^(1/3) - 2*(1 + Log[x]^6)^(2/3))/(6*x*Log[x]*(1 + Log[x]^6)^(2/3)) + (Log[x]^4 + Log[x]^2*(1 + Log[x]^6)^(1/3) + (1 + Log[x]^`
- `cM-2-157` (unknown): larger 105 -> 619 — `(2*x*(3 + a + Log[1 + x^2])^(3/2))/((1 + x^2)*(2 + Log[1 + x^2]^2))`
- `cM-2-194` (unknown): slower 0.2s -> 10.2s — `2/(1 + a + E^(6*x))^(1/4)`
- `cM-2-197` (unknown): larger 114 -> 515 — `2/(2 + a + E^(2*x))^(1/3)`
- `cM-2-286` (unknown): larger 120 -> 804 — `(E^x*(2 + E^(-x))^(3/2))/(-1 + a + E^x - E^(2*x))`
- `cM-2-307` (unknown): slower 0.3s -> 21.2s — `((2 + Log[x])*(1 + Log[x] + Log[x]^2)^(3/2))/(x*(-1 + a + Log[x] + Log[x]^2))`
- `cM-2-329` (unknown): larger 279 -> 442 — `(2*E^x^2*Sqrt[1 + a + E^x^2]*x)/(2 + E^(2*x^2))`
- `cM-2-352` (unknown): larger 56 -> 374 — `(E^x*Sqrt[2 + Log[1 + E^x]])/((1 + E^x)*(1 + a - Log[1 + E^x] + 2*Log[1 + E^x]^2))`
- `cM-3-10` (elem): larger 51 -> 428 — `E^x/(2*(2 + a + E^x)^(1/3) + a^2*(2 + a + E^x)^(1/3) + 2*a*E^x*(2 + a + E^x)^(1/3) + E^(2*x)*(2 + a + E^x)^(1/3))`
- `cM-3-11` (unknown): slower 0.0s -> 10.0s — `-(1/(1 + a + E^(-2*x))^(1/4))`
- `cM-3-21` (elem): larger 1082 -> 4517 — `E^x/(1 + E^x)^(5/6) - (-1/6*E^x/(1 + E^x)^(5/6) + E^x/(3*(1 + E^x)^(2/3)))/(4*(1 - (1 + E^x)^(1/6) + (1 + E^x)^(1/3))) + (E^x/(6*(1 + E^x)^(5/6)) + E^x/(3*(1 + `
- `cM-3-158` (unknown): larger 64 -> 523 — `1/(x*Sqrt[1 + a + Log[x]]*(1 + Log[x]^2))`
- `cM-3-163` (unknown): larger 157 -> 1608 — `(2*(-1 + Log[1 + 2*x]))/((1 + 2*x)*(2 + a + Log[1 + 2*x])^(3/2)*(2 + Log[1 + 2*x] + Log[1 + 2*x]^2))`
- `cM-3-173` (unknown): larger 257 -> 517 — `2/(3 + a + E^(4*x))^(1/3)`
- `cM-3-222` (unknown): larger 128 -> 529 — `(2*x)/(3 + a + E^(4*x^2))^(1/3)`
- `cM-3-261` (unknown): larger 35 -> 187 — `E^x/((3 + E^x)^(1/3)*(2 + a + E^x))`
- `cM-3-266` (unknown): lost ELEM -> TIMEOUT — `(-2*E^x^2*(2 + E^x^2)*Sqrt[2 + E^x^2 + E^(2*x^2)]*x)/(2 + (2 + a)*E^x^2 - E^(2*x^2))`
- `cM-3-278` (unknown): larger 586 -> 935, slower 0.1s -> 20.6s — `((2 + E^x)*(1 + a + E^x)^(2/3))/(-2 + E^x)`
- `cM-3-347` (unknown): larger 37 -> 214 — `(2*E^x^2*x)/((1 + E^x^2)^(1/3)*(-2 + a + E^x^2))`
- `cM-3-369` (unknown): larger 47 -> 381 — `E^x/(Sqrt[1 + E^x]*(3 + a + E^(2*x)))`
- `cM-4-68` (unknown): larger 122 -> 517 — `1/(3 + a + E^(4*x))^(1/3)`
- `cM-4-114` (unknown): larger 80 -> 196 — `(3 + Log[x])^(1/3)/(x*(2 + a + Log[x]))`
- `cM-4-167` (unknown): larger 86 -> 851 — `1/(x*(3 + a + Log[x]^2)*Sqrt[1 + Log[x] + Log[x]^2])`
- `cM-4-181` (elem): larger 53 -> 870 — `E^x/(3*(2 + a + E^x)^(1/4) + a^2*(2 + a + E^x)^(1/4) + 2*a*E^x*(2 + a + E^x)^(1/4) + E^(2*x)*(2 + a + E^x)^(1/4))`
- `cM-4-200` (unknown): larger 139 -> 1374, slower 0.1s -> 7.0s — `(2*E^(2*x))/((2 + a + E^(2*x))^(3/2)*(2 + E^(4*x)))`
- `cB-7-14` (unknown): larger 304 -> 816 — `x*(-2 + (3 + a)*x)^(3/2)*Log[-1 + 3*x]`
- `cB-9-85` (elem): larger 71 -> 451 — `(672 + 372*x - 300*x^2 + 496*(3 + 3*x) - 800*x*(3 + 3*x))/(2025*(3 + 3*x)^(3/4)) + (-4*Log[2 + 2*x] + x*Log[2 + 2*x] + 5*x^2*Log[2 + 2*x])/(15*(3 + a + 3*x)^(3/`
- `cB-10-152` (elem): lost ELEM -> TIMEOUT — `-3/(2*(-3 + 3*x)^(1/3)) + (-3 + 3*x)^(2/3)/(-1 + 2*x) + (27 - 27*a - 54*x + 27*a*x + 27*x^2)/(-36*(-3 + a + 3*x)^(1/3) + 108*x*(-3 + a + 3*x)^(1/3) - 108*x^2*(-`
- `cB-10-224` (unknown): larger 193 -> 408 — `(x*Log[1 + 3*x])/Sqrt[2 + (3 + a)*x]`
- `cB-11-127` (unknown): larger 237 -> 438 — `x*Sqrt[2 + (3 + a)*x]*Log[2 + 2*x]`
- `cB-11-170` (unknown): larger 237 -> 438 — `x*Sqrt[-1 + (3 + a)*x]*Log[2 + 2*x]`

## Better

- `cB-7-107` (unknown): TIMEOUT -> ELEM — `x*(-2 + (2 + a)*x)^(1/4)*Log[2 + 3*x]`
- `cB-7-178` (unknown): TIMEOUT -> ELEM — `x*(2 + (2 + a)*x)^(1/3)*Log[1 + 3*x]`
- `cB-7-190` (unknown): TIMEOUT -> ELEM — `(-1 + (2 + a)*x)^(2/3)*Log[-1 + 3*x]`
- `cB-7-230` (unknown): TIMEOUT -> ELEM — `x*(-2 + (3 + a)*x)^(2/3)*Log[1 + x]`
- `cB-7-240` (unknown): UNEVAL -> ELEM — `(-1 + (2 + a)*x)^(2/3)*Log[-1 + x]`
- `cB-7-269` (unknown): UNEVAL -> ELEM — `x*(-1 + (2 + a)*x)^(2/3)*Log[-2 + x]`
- `cB-8-4` (unknown): TIMEOUT -> ELEM — `(-2 + (2 + a)*x)^(1/3)*Log[3 + 2*x]`
- `cB-8-24` (unknown): TIMEOUT -> ELEM — `x*(-1 + (2 + a)*x)^(1/4)*Log[2 + 3*x]`
- `cB-8-58` (unknown): TIMEOUT -> ELEM — `x*(3 + (2 + a)*x)^(2/3)*Log[-1 + 3*x]`
- `cB-8-65` (elem): TIMEOUT -> ELEM — `(18*x^2*(-1 + a + x)^(1/6))/(-26 + 39*x) - (108*(-1 + x)^(1/6))/(-182 + 273*x) - (18*(-1 + x)^(1/6)*x)/(-182 + 273*x) + (7438 - 4392*(-1 + x) - 732*x - 10584*(-`
- `cB-8-156` (unknown): TIMEOUT -> ELEM — `(-2 + (2 + a)*x)^(2/3)*Log[-1 + 3*x]`
- `cB-8-159` (unknown): TIMEOUT -> ELEM — `(-1 + (2 + a)*x)^(1/6)*Log[-1 + x]`
- `cB-8-258` (unknown): TIMEOUT -> ELEM — `x*(-1 + (3 + a)*x)^(1/4)*Log[1 + 2*x]`
- `cB-8-261` (unknown): TIMEOUT -> ELEM — `x*(-1 + (3 + a)*x)^(1/4)*Log[3 + 2*x]`
- `cB-9-71` (unknown): UNEVAL -> ELEM — `x*(-1 + a + x)^(1/4)*Log[-3 + x]`
- `cB-9-102` (unknown): TIMEOUT -> ELEM — `x*(-1 + (2 + a)*x)^(2/3)*Log[2 + 3*x]`
- `cB-9-125` (unknown): TIMEOUT -> ELEM — `x*(2 + (3 + a)*x)^(2/3)*Log[1 + 3*x]`
- `cB-9-273` (unknown): TIMEOUT -> ELEM — `(-2 + (3 + a)*x)^(2/3)*Log[-1 + x]`
- `cB-9-294` (unknown): TIMEOUT -> ELEM — `x*(-3 + (3 + a)*x)^(1/4)*Log[-3 + 2*x]`
- `cB-10-2` (elem): TIMEOUT -> ELEM — `(6*x^2*(1 + 3*x)^(1/3))/(7 + 14*x) - (1 + 3*x)^(1/3)/(14 + 28*x) + (x*(1 + 3*x + a*x)^(1/3))/(14 + 28*x) + (-79 + 93*x - 144*x^2 + 93*(1 + 3*x) - 288*x*(1 + 3*x`
- `cB-10-40` (elem): TIMEOUT -> ELEM — `(9*x^2*(-1 + a + x)^(1/3))/(-7 + 21*x) - (27*(-1 + x)^(1/3))/(-28 + 84*x) - (9*(-1 + x)^(1/3)*x)/(-28 + 84*x) + (857 + 45*(-1 + x) + 15*x - 864*(-1 + x)*x - 144`
- `cB-10-198` (unknown): TIMEOUT -> ELEM — `(3 + (3 + a)*x)^(2/3)*Log[-2 + 2*x]`
- `cB-10-234` (elem): TIMEOUT -> ELEM — `(3*x^2*(2 + a + x)^(2/3))/(-16 + 8*x) - (9*(2 + x)^(2/3))/(-20 + 10*x) + (3*x*(2 + x)^(2/3))/(-20 + 10*x) + (-1452 - 396*x - 75*x^2 - 594*(2 + x) - 225*x*(2 + x`
- `cB-10-240` (unknown): TIMEOUT -> ELEM — `x*(1 + (2 + a)*x)^(2/3)*Log[-2 + x]`
- `cB-11-9` (elem): TIMEOUT -> ELEM — `(9*x^2*(-1 + a + x)^(2/3))/(16 + 24*x) - (9*(-1 + x)^(2/3)*x)/(40 + 60*x) - (27*(-1 + x)^(2/3))/(80 + 120*x) + (-349 + 711*(-1 + x) + 474*x - 675*(-1 + x)*x - 2`
- `cB-11-166` (unknown): TIMEOUT -> ELEM — `(-3 + (3 + a)*x)^(2/3)*Log[3 + 3*x]`
- `cB-11-172` (unknown): TIMEOUT -> ELEM — `(-3 + (3 + a)*x)^(1/6)*Log[-2 + 3*x]`
