# 258 — parametric A/B of the Trager and Bronstein cascade stages

11 cases; A = `Integrate`$CascadeParameters = False`, B = `{"Trager"}`; timeout 30 s, check at a = 1/7.

- hard failures in B: **1**
- B worse than A: **3**
- B better than A: **0**

| A | B | cases |
|---|---|---|
| ELEM | ELEM | 6 |
| TIMEOUT | TIMEOUT | 4 |
| WRONG | WRONG | 1 |

## Hard failures

- `trager-cheb-002` (unknown): WRONG — `x^2*Sqrt[-1 + a + 2*x^4]`

## Worse

- `trager-evensub-005` (unknown): larger 107 -> 280 — `x/(Sqrt[1 + a + x^4]*(1 + x^2 + x^4))`
- `cA-1-263` (elem): slower 0.0s -> 10.0s — `(2*a*x^2 - 2*a*x^2*(1 + x^3)^(1/6) - 2*a*x^2*(1 + x^3)^(1/3) - a^2*x^2*(1 + x^3)^(1/3) + 4*a*x^2*Sqrt[1 + x^3] + a^2*x^2*Sqrt[1 + x^3] - 2*a*x^2*(1 + x^3)^(2/3)`
- `cM-2-67` (unknown): larger 78 -> 179 — `-((E^x*(-2 + E^x))/(Sqrt[2 + a + E^x]*(-1 + E^x - 2*E^(2*x))))`

## Better

