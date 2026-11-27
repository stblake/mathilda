# 265-goursat-stress (2026-11-27): correctness pass; timings pending a quiet machine

**Correctness pass; timings pending a quiet machine.** This run asks only whether an answer is right: wrong answers and verdicts do not depend on machine load, times and timeouts do. The times below are NOT a baseline, and every timeout is provisional (it may close on a quiet machine, or within the timed run's own limit). The timed baseline is `python3 tools/goursat/goursat_stress.py` on a quiet machine.

1183 cases. Mathilda 0.6764; per-case timeout 60 s for Mathilda and sympy, 4 parallel processes; wall 156.7 min.
Machine load (1-min average, 16 CPUs): start 6.2, mean 17.1, max 162.4, end 4.9 -- other agents were building and testing on the same machine, so times are noisy upper bounds. The src/internal .m stamp changed during the run (a package was edited mid-run).

Every Mathilda and sympy answer is judged by PossibleZeroQ[D[F, x] - f] in Mathilda: zero is `verified`, nonzero is `wrong`, anything else is `unverified`. `cert`: certified non-elementary. goursat.py is read from tests/data/generalised_goursat/sympy_*.tsv (bench_sympy.py, 120 s limit, its own numerical verification; no answer text, so not PossibleZeroQ-judged).

Times: `median s` over decided runs (verified answer or certificate), `max s` over all runs (a timeout counts its full limit). `x goursat.py`: median per-case time ratio method / goursat.py over the cases both decided (in brackets: total-time ratio, and the number of such cases); below 1 means faster than goursat.py.

## All cases

| run | n | verified | wrong | unverified | cert | false cert | uneval | timeout | error | median s | max s | median leaves | x goursat.py |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 1183 | 804 (68 special) | 0 | 0 | 137 | 0 | 94 | 148 | 0 | 0.576 | 143 | 136 | 0.296 (0.738, 896) |
| Mathilda:GoursatAlgebraic | 1183 | 275 | 0 | 0 | 0 | 0 | 883 | 25 | 0 | 0.0669 | 60.4 | 103 | 0.108 (0.308, 274) |
| Mathilda:GeneralisedGoursat | 1183 | 874 (68 special) | 1 | 116 | 141 | 0 | 19 | 32 | 0 | 0.125 | 123 | 112 | 0.0782 (0.368, 981) |
| goursat.py | 1183 | 913 | 1 | 0 | 211 | 0 | 5 | 35 | 18 | 2.35 | 120 | - | - |

goursat.py `verified` = its own numerical check passed (defect <= 1e-20).

## By family

Families pool the corpora they occur in (the stress families pool stress_2026, stress_7 and stress_11).

| family | expect | run | n | verified | wrong | unverified | cert | false cert | uneval | timeout | error | median s | max s | median leaves | x goursat.py |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev | elem 56 nonelem 48 | Mathilda:Automatic | 104 | 104 (48 special) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.00447 | 2.23 | 40 | 0.0252 (0.844, 104) |
| chebyshev | elem 56 nonelem 48 | Mathilda:GoursatAlgebraic | 104 | 38 | 0 | 0 | 0 | 0 | 66 | 0 | 0 | 0.0212 | 0.0579 | 36 | 0.135 (0.13, 38) |
| chebyshev | elem 56 nonelem 48 | Mathilda:GeneralisedGoursat | 104 | 92 (36 special) | 0 | 0 | 12 | 0 | 0 | 0 | 0 | 0.0455 | 2.67 | 54 | 0.251 (1.07, 104) |
| chebyshev | elem 56 nonelem 48 | goursat.py | 104 | 56 | 0 | 0 | 48 | 0 | 0 | 0 | 0 | 0.204 | 4.15 | - | - |
| unit | elem 17 nonelem 4 unknown 2 | Mathilda:Automatic | 23 | 19 (2 special) | 0 | 0 | 4 | 0 | 0 | 0 | 0 | 0.219 | 13.4 | 50 | 0.374 (1.4, 21) |
| unit | elem 17 nonelem 4 unknown 2 | Mathilda:GoursatAlgebraic | 23 | 5 | 0 | 0 | 0 | 0 | 18 | 0 | 0 | 0.0402 | 6.07 | 71 | 0.0471 (0.131, 5) |
| unit | elem 17 nonelem 4 unknown 2 | Mathilda:GeneralisedGoursat | 23 | 17 | 0 | 0 | 4 | 0 | 2 | 0 | 0 | 0.0269 | 19.7 | 58 | 0.135 (0.852, 21) |
| unit | elem 17 nonelem 4 unknown 2 | goursat.py | 23 | 17 | 0 | 0 | 4 | 0 | 2 | 0 | 0 | 0.57 | 7.58 | - | - |
| cubic-field | elem 1 | Mathilda:Automatic | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.706 | 0.706 | 168 | 1.16 (1.16, 1) |
| cubic-field | elem 1 | Mathilda:GoursatAlgebraic | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | - | 6.14 | - | - |
| cubic-field | elem 1 | Mathilda:GeneralisedGoursat | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | - | 0.184 | - | - |
| cubic-field | elem 1 | goursat.py | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.61 | 0.61 | - | - |
| examples | unknown 29 | Mathilda:Automatic | 29 | 24 (5 special) | 0 | 0 | 5 | 0 | 0 | 0 | 0 | 0.654 | 15.5 | 86.5 | 0.277 (1.91, 28) |
| examples | unknown 29 | Mathilda:GoursatAlgebraic | 29 | 10 | 0 | 0 | 0 | 0 | 19 | 0 | 0 | 0.0398 | 6.13 | 109 | 0.165 (0.357, 10) |
| examples | unknown 29 | Mathilda:GeneralisedGoursat | 29 | 19 (1 special) | 0 | 0 | 8 | 0 | 2 | 0 | 0 | 0.0515 | 19.3 | 125 | 0.127 (0.58, 27) |
| examples | unknown 29 | goursat.py | 29 | 19 | 0 | 0 | 9 | 0 | 1 | 0 | 0 | 0.578 | 4.99 | - | - |
| literature | elem 32 unknown 11 nonelem 8 | Mathilda:Automatic | 51 | 42 (8 special) | 0 | 0 | 9 | 0 | 0 | 0 | 0 | 0.211 | 10.3 | 57.5 | 0.235 (0.634, 49) |
| literature | elem 32 unknown 11 nonelem 8 | Mathilda:GoursatAlgebraic | 51 | 20 | 0 | 0 | 0 | 0 | 31 | 0 | 0 | 0.0363 | 1.73 | 78.5 | 0.0805 (0.0536, 20) |
| literature | elem 32 unknown 11 nonelem 8 | Mathilda:GeneralisedGoursat | 51 | 41 (7 special) | 0 | 0 | 8 | 0 | 2 | 0 | 0 | 0.0555 | 21.1 | 100 | 0.117 (0.877, 49) |
| literature | elem 32 unknown 11 nonelem 8 | goursat.py | 51 | 34 | 0 | 0 | 15 | 0 | 2 | 0 | 0 | 0.464 | 11.4 | - | - |
| welz | elem 30 | Mathilda:Automatic | 30 | 30 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 10.1 | 10.8 | 144 | 9.83 (11.1, 30) |
| welz | elem 30 | Mathilda:GoursatAlgebraic | 30 | 30 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0796 | 0.188 | 144 | 0.069 (0.0957, 30) |
| welz | elem 30 | Mathilda:GeneralisedGoursat | 30 | 30 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.293 | 0.653 | 141 | 0.291 (0.31, 30) |
| welz | elem 30 | goursat.py | 30 | 30 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.02 | 1.44 | - | - |
| hermite-exact | elem 75 | Mathilda:Automatic | 75 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0206 | 0.0696 | 26 | 0.17 (0.174, 60) |
| hermite-exact | elem 75 | Mathilda:GoursatAlgebraic | 75 | 4 | 0 | 0 | 0 | 0 | 71 | 0 | 0 | 0.221 | 6.07 | 19 | 6.05 (5.63, 4) |
| hermite-exact | elem 75 | Mathilda:GeneralisedGoursat | 75 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0134 | 0.0518 | 23 | 0.0927 (0.121, 60) |
| hermite-exact | elem 75 | goursat.py | 75 | 60 | 0 | 0 | 0 | 0 | 0 | 0 | 15 | 0.135 | 0.388 | - | - |
| cbrt-symmetric | elem 75 | Mathilda:Automatic | 75 | 21 | 0 | 0 | 0 | 0 | 0 | 54 | 0 | 5.37 | 139 | 3.07e+03 | 1.6 (2.46, 21) |
| cbrt-symmetric | elem 75 | Mathilda:GoursatAlgebraic | 75 | 8 | 0 | 0 | 0 | 0 | 65 | 2 | 0 | 3.01 | 60.2 | 1.24e+03 | 0.48 (1.24, 8) |
| cbrt-symmetric | elem 75 | Mathilda:GeneralisedGoursat | 75 | 61 | 0 | 6 | 0 | 0 | 0 | 8 | 0 | 1.35 | 91.4 | 405 | 0.175 (0.441, 59) |
| cbrt-symmetric | elem 75 | goursat.py | 75 | 72 | 1 | 0 | 0 | 0 | 0 | 2 | 0 | 9.45 | 120 | - | - |
| cbrt-mixed | elem 75 | Mathilda:Automatic | 75 | 16 | 0 | 0 | 0 | 0 | 0 | 59 | 0 | 9.54 | 143 | 1.02e+04 | 3.08 (2.04, 16) |
| cbrt-mixed | elem 75 | Mathilda:GoursatAlgebraic | 75 | 1 | 0 | 0 | 0 | 0 | 74 | 0 | 0 | 0.303 | 7.66 | 254 | 0.153 (0.153, 1) |
| cbrt-mixed | elem 75 | Mathilda:GeneralisedGoursat | 75 | 69 | 0 | 3 | 0 | 0 | 0 | 3 | 0 | 4.56 | 123 | 334 | 0.214 (0.506, 66) |
| cbrt-mixed | elem 75 | goursat.py | 75 | 72 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 14.9 | 120 | - | - |
| cbrt-nonelementary | nonelem 75 | Mathilda:Automatic | 75 | 0 | 0 | 0 | 51 | 0 | 0 | 24 | 0 | 10.6 | 72.6 | - | 0.516 (0.485, 51) |
| cbrt-nonelementary | nonelem 75 | Mathilda:GoursatAlgebraic | 75 | 0 | 0 | 0 | 0 | 0 | 75 | 0 | 0 | - | 8.25 | - | - |
| cbrt-nonelementary | nonelem 75 | Mathilda:GeneralisedGoursat | 75 | 7 (7 special) | 0 | 0 | 63 | 0 | 0 | 5 | 0 | 10.2 | 88.7 | 277 | 0.532 (0.532, 67) |
| cbrt-nonelementary | nonelem 75 | goursat.py | 75 | 0 | 0 | 0 | 72 | 0 | 0 | 3 | 0 | 17.9 | 120 | - | - |
| sqrt-klein | elem 75 | Mathilda:Automatic | 75 | 73 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0.578 | 60.8 | 273 | 0.107 (0.148, 65) |
| sqrt-klein | elem 75 | Mathilda:GoursatAlgebraic | 75 | 16 | 0 | 0 | 0 | 0 | 37 | 22 | 0 | 0.151 | 60.4 | 82.5 | 0.127 (1.33, 15) |
| sqrt-klein | elem 75 | Mathilda:GeneralisedGoursat | 75 | 53 | 0 | 14 | 0 | 0 | 4 | 4 | 0 | 0.0833 | 60.4 | 235 | 0.031 (0.0627, 50) |
| sqrt-klein | elem 75 | goursat.py | 75 | 65 | 0 | 0 | 0 | 0 | 0 | 8 | 2 | 5.33 | 120 | - | - |
| sqrt-mixed | elem 75 | Mathilda:Automatic | 75 | 72 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0.515 | 60.1 | 221 | 0.14 (0.352, 66) |
| sqrt-mixed | elem 75 | Mathilda:GoursatAlgebraic | 75 | 6 | 0 | 0 | 0 | 0 | 68 | 1 | 0 | 0.321 | 60.2 | 494 | 0.254 (0.466, 6) |
| sqrt-mixed | elem 75 | Mathilda:GeneralisedGoursat | 75 | 58 | 0 | 12 | 0 | 0 | 1 | 4 | 0 | 0.0338 | 60.6 | 92 | 0.0406 (0.0732, 57) |
| sqrt-mixed | elem 75 | goursat.py | 75 | 67 | 0 | 0 | 0 | 0 | 0 | 8 | 0 | 2.08 | 120 | - | - |
| sqrt-nonelementary | nonelem 75 | Mathilda:Automatic | 75 | 5 (5 special) | 0 | 0 | 68 | 0 | 0 | 2 | 0 | 10.5 | 60.1 | 268 | 1.29 (0.98, 62) |
| sqrt-nonelementary | nonelem 75 | Mathilda:GoursatAlgebraic | 75 | 0 | 0 | 0 | 0 | 0 | 75 | 0 | 0 | - | 4.33 | - | - |
| sqrt-nonelementary | nonelem 75 | Mathilda:GeneralisedGoursat | 75 | 17 (17 special) | 0 | 4 | 46 | 0 | 2 | 6 | 0 | 0.299 | 60.8 | 271 | 0.128 (0.142, 56) |
| sqrt-nonelementary | nonelem 75 | goursat.py | 75 | 0 | 0 | 0 | 63 | 0 | 0 | 11 | 1 | 3.75 | 120 | - | - |
| sqrt-order3 | elem 75 | Mathilda:Automatic | 75 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.253 | 43.3 | 140 | 0.239 (0.398, 75) |
| sqrt-order3 | elem 75 | Mathilda:GoursatAlgebraic | 75 | 8 | 0 | 0 | 0 | 0 | 67 | 0 | 0 | 0.23 | 7.03 | 66 | 0.262 (0.307, 8) |
| sqrt-order3 | elem 75 | Mathilda:GeneralisedGoursat | 75 | 67 | 0 | 4 | 0 | 0 | 4 | 0 | 0 | 0.0271 | 10.2 | 58 | 0.0288 (0.034, 67) |
| sqrt-order3 | elem 75 | goursat.py | 75 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.07 | 47.8 | - | - |
| cbrt-second-level | elem 75 | Mathilda:Automatic | 75 | 71 | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 38.5 | 60 | 7.36e+03 | 0.879 (0.914, 71) |
| cbrt-second-level | elem 75 | Mathilda:GoursatAlgebraic | 75 | 0 | 0 | 0 | 0 | 0 | 75 | 0 | 0 | - | 0.0549 | - | - |
| cbrt-second-level | elem 75 | Mathilda:GeneralisedGoursat | 75 | 0 | 0 | 73 | 0 | 0 | 0 | 2 | 0 | - | 60.4 | - | - |
| cbrt-second-level | elem 75 | goursat.py | 75 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 41.6 | 50.9 | - | - |
| pseudoelliptic-sqrt | elem 152 | Mathilda:Automatic | 152 | 139 | 0 | 0 | 0 | 0 | 13 | 0 | 0 | 0.177 | 16.1 | 91 | 0.198 (0.432, 139) |
| pseudoelliptic-sqrt | elem 152 | Mathilda:GoursatAlgebraic | 152 | 93 | 0 | 0 | 0 | 0 | 59 | 0 | 0 | 0.0714 | 6.07 | 111 | 0.102 (0.124, 93) |
| pseudoelliptic-sqrt | elem 152 | Mathilda:GeneralisedGoursat | 152 | 150 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0.0321 | 11.6 | 64 | 0.0335 (0.0905, 150) |
| pseudoelliptic-sqrt | elem 152 | goursat.py | 152 | 152 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.961 | 7.46 | - | - |
| pseudoelliptic-cbrt | elem 118 | Mathilda:Automatic | 118 | 37 | 0 | 0 | 0 | 0 | 81 | 0 | 0 | 0.346 | 56.4 | 176 | 0.0598 (0.277, 37) |
| pseudoelliptic-cbrt | elem 118 | Mathilda:GoursatAlgebraic | 118 | 36 | 0 | 0 | 0 | 0 | 82 | 0 | 0 | 0.188 | 6.06 | 172 | 0.0863 (0.095, 36) |
| pseudoelliptic-cbrt | elem 118 | Mathilda:GeneralisedGoursat | 118 | 118 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.308 | 8.94 | 390 | 0.0602 (0.123, 118) |
| pseudoelliptic-cbrt | elem 118 | goursat.py | 118 | 118 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.42 | 17.7 | - | - |

## By corpus (GeneralisedGoursat vs goursat.py)

| corpus | n | GG verified | GG cert | GG timeout | goursat.py elementary | goursat.py non-elementary | goursat.py timeout | GG median s | goursat.py median s | x goursat.py |
|---|---|---|---|---|---|---|---|---|---|---|
| unit | 128 | 109 | 16 | 0 | 74 | 52 | 0 | 0.0422 | 0.241 | 0.208 (1, 125) |
| examples | 29 | 19 | 8 | 0 | 19 | 9 | 0 | 0.0515 | 0.578 | 0.127 (0.58, 27) |
| literature | 51 | 41 | 8 | 0 | 34 | 15 | 0 | 0.0555 | 0.464 | 0.117 (0.877, 49) |
| welz | 30 | 30 | 0 | 0 | 30 | 0 | 0 | 0.293 | 1.02 | 0.291 (0.31, 30) |
| stress_2026 | 225 | 149 | 34 | 2 | 163 | 44 | 12 | 0.227 | 5.33 | 0.13 (0.417, 169) |
| stress_7 | 225 | 131 | 38 | 10 | 159 | 46 | 14 | 0.22 | 7.51 | 0.0761 (0.317, 158) |
| stress_11 | 225 | 127 | 37 | 20 | 164 | 45 | 9 | 0.237 | 7.4 | 0.0907 (0.461, 155) |
| pseudoelliptic_sqrt | 152 | 150 | 0 | 0 | 152 | 0 | 0 | 0.0321 | 0.961 | 0.0335 (0.0905, 150) |
| pseudoelliptic_cbrt | 118 | 118 | 0 | 0 | 118 | 0 | 0 | 0.308 | 4.42 | 0.0602 (0.123, 118) |

## Wrong answers from Mathilda (1)

An answer PossibleZeroQ[D[F, x] - f] rejects (or a garbage answer carrying $Failed or a package-private symbol). Each is a wrong answer.

- `G-psqrt-pseudoelliptic-sqrt-001` (pseudoelliptic-sqrt #1: (x/(x^2 + 2*Sqrt[2] + 3))/Sqrt[x^3 + x]) Mathilda:GeneralisedGoursat, 11.65 s: `x/((3 + 2*Sqrt[2] + x^2)*Sqrt[x + x^3])` -> `((-1 + Sqrt[2])*ArcTan[-1 + Sqrt[x + x^3]/x])/4 + ((-1 + Sqrt[2])*ArcTan[1 + Sqrt[x + x^3]/x])/4 + (1/8 + (Sqrt[-1/48 - 13/(768*((-35/2048 + ((9*I)/1024)*Sqrt[3])/2)^(1/3)) - ((-35/2048 + ((9*I)/1024)*Sqrt[3])/2)^(1/3)/3] - Sqrt[-1/24 + 13/`

## False certificates (0)

A non-elementary certificate on a case the corpus constructs as elementary.


## Unverified Mathilda answers (116)

PossibleZeroQ did not decide (timed out, or D[F, x] - f does not evaluate numerically).

- `G-s2026-cbrt-symmetric-000` Mathilda:GeneralisedGoursat (undecided-nonnumeric): `(13608 - 835280*x + 15080112*x^2 - 68724738*x^3 + 165022002*x^4 - 313835700*x^5 + 468801960*x^6 - 492763974*x^7 + 433252686*x^8 - 262055814*x^9 + 65842898*x^10)`
- `G-s2026-cbrt-symmetric-001` Mathilda:GeneralisedGoursat (undecided): `(3*x + 163*x^2 - 2432*x^3 + 4446*x^4 + 58752*x^5 - 239598*x^6 - 130734*x^7 + 1829061*x^8 - 2332071*x^9 - 1233468*x^10 + 4133430*x^11 - 2125764*x^12)/((-x + 3*x^`
- `G-s2026-cbrt-symmetric-006` Mathilda:GeneralisedGoursat (undecided): `(5180 - 29512*x + 61612*x^2 - 58511*x^3 + 26292*x^4 - 5162*x^5 + 356*x^6 + x^7)/((2 - x - 2*x^2 + x^3)^(1/3)*(-1924 + 8686*x - 13050*x^2 + 4510*x^3 + 8702*x^4 -`
- `G-s2026-cbrt-mixed-012` Mathilda:GeneralisedGoursat (undecided): `(7584813 + 40801293*x - 19179513*x^2 - 243182481*x^3 - 489361824*x^4 - 518452864*x^5 - 230259700*x^6 + 217588308*x^7 + 418793296*x^8 + 207985664*x^9 - 48712000*`
- `G-s2026-sqrt-klein-005` Mathilda:GeneralisedGoursat (undecided): `(-721884960 - 3774719232*x - 8432856864*x^2 - 11122859088*x^3 - 9913088880*x^4 - 6412943988*x^5 - 3128526576*x^6 - 1168160364*x^7 - 332404011*x^8 - 70362216*x^9`
- `G-s2026-sqrt-klein-011` Mathilda:GeneralisedGoursat (undecided): `(66 + 420*x + 840*x^2 + 876*x^3 + 780*x^4 + 876*x^5 + 840*x^6 + 420*x^7 + 66*x^8)/(Sqrt[-9 - 9*x + 9*x^3 + 9*x^4]*(65 + 144*x + 164*x^2 + 176*x^3 + 198*x^4 + 17`
- `G-s2026-sqrt-mixed-005` Mathilda:GeneralisedGoursat (undecided): `(-92630952 + 156368628*x + 53693424*x^2 - 335291592*x^3 + 448732347*x^4 - 375517836*x^5 + 231250450*x^6 - 104693295*x^7 + 33953811*x^8 - 7757100*x^9 + 1203480*x`
- `G-s2026-sqrt-mixed-017` Mathilda:GeneralisedGoursat (undecided): `(-506536 + 595983*x + 10537746*x^2 - 58046913*x^3 + 137022390*x^4 - 143106531*x^5 - 13370482*x^6 + 203959401*x^7 - 243568518*x^8 + 161272620*x^9 - 75805560*x^10`
- `G-s2026-sqrt-nonelementary-006` Mathilda:GeneralisedGoursat (undecided): `(164 + 2192*x - 3336*x^2 - 1142*x^3 + 2866*x^4 + 1056*x^5 - 1117*x^6 + 674*x^7 + 990*x^8 + 189*x^9)/(Sqrt[1 + x^3]*(192 - 336*x + 384*x^3 - 294*x^4 + 192*x^6 + `
- `G-s2026-sqrt-nonelementary-019` Mathilda:GeneralisedGoursat (undecided): `(-636979 + 1885356*x - 3231430*x^2 + 3921552*x^3 - 3396249*x^4 + 2358372*x^5 - 1456359*x^6 + 736200*x^7 - 261651*x^8 + 58438*x^9 - 7281*x^10 + 385*x^11)/(Sqrt[-`
- `G-s2026-sqrt-order3-015` Mathilda:GeneralisedGoursat (undecided): `(48492 + 146598*x + 192471*x^2 + 194108*x^3 + 192669*x^4 + 148716*x^5 + 75890*x^6 + 25272*x^7 + 5589*x^8 + 882*x^9 + 81*x^10)/(Sqrt[-9 - 3*x - 3*x^2 - x^3]*(-89`
- `G-s2026-cbrt-second-level-000` Mathilda:GeneralisedGoursat (undecided): `(2837234832666842365952*x - 2040647904553557032960*x^4 + 449249098763385962496*x^7 - 47498262429924589568*x^10 + 2968641401870286848*x^13 - 109679955752779776*x`
- `G-s2026-cbrt-second-level-001` Mathilda:GeneralisedGoursat (undecided): `(-55831680773263369854517248 + 34899007672401088563118080*x^3 - 4144783073051274570104832*x^6 + 279514838216553711796224*x^9 - 17468855790792609890304*x^12 + 10`
- `G-s2026-cbrt-second-level-002` Mathilda:GeneralisedGoursat (undecided): `(-1395317740800*x + 7751763499776*x^4 - 13520841615360*x^7 + 11333568430752*x^10 - 5666784215376*x^13 + 1690105201920*x^16 - 242242609368*x^19 + 10900919850*x^2`
- `G-s2026-cbrt-second-level-003` Mathilda:GeneralisedGoursat (undecided): `(-113490905713920*x - 638911731942144*x^4 - 1118473211492352*x^7 - 940759885920672*x^10 - 470379942960336*x^13 - 139809151436544*x^16 - 19965991623192*x^19 - 88`
- `G-s2026-cbrt-second-level-004` Mathilda:GeneralisedGoursat (undecided): `(781365958933184671664217414132864*x - 166410096301273264453930637537280*x^4 + 10852978868366899853722988693952*x^7 - 339935178011661054724804550064*x^10 + 6295`
- `G-s2026-cbrt-second-level-005` Mathilda:GeneralisedGoursat (undecided): `(6203835880306979761029120*x - 4459343618608256035848192*x^4 + 981560887542835788644352*x^7 - 103762670732626633949184*x^10 + 6485166920789164621824*x^13 - 2396`
- `G-s2026-cbrt-second-level-006` Mathilda:GeneralisedGoursat (undecided): `(-55828090782832820666499072 + 34906861550158238372069376*x^3 - 4146978064360304120168448*x^6 + 279698072583276173721600*x^9 - 17478337002346534404096*x^12 + 10`
- `G-s2026-cbrt-second-level-007` Mathilda:GeneralisedGoursat (undecided): `(-1962443366013835363221504*x + 1409512455230955769036800*x^4 - 310185226002052386127872*x^7 + 32783767802822042320896*x^10 - 2048985487676377645056*x^13 + 7572`
- `G-s2026-cbrt-second-level-008` Mathilda:GeneralisedGoursat (undecided): `(112036264704*x - 652507512480*x^4 + 1153082418048*x^7 - 978301038942*x^10 + 489150519471*x^13 - 144135302256*x^16 + 20390859765*x^19 - 875283318*x^22)/((-2 + x`
- `G-s2026-cbrt-second-level-009` Mathilda:GeneralisedGoursat (undecided): `(-12621091852800 - 64774581023232*x^3 - 63214147961856*x^6 - 34489664832384*x^9 - 17083162657728*x^12 - 8622416208096*x^15 - 3950884247616*x^18 - 1012102828488*`
- `G-s2026-cbrt-second-level-010` Mathilda:GeneralisedGoursat (undecided): `(-2613753608316561131381494272 + 484056111568281313384502784*x^3 - 17032628589098571158473728*x^6 + 340329016337118838922880*x^9 - 6302248645341451432896*x^12 +`
- `G-s2026-cbrt-second-level-011` Mathilda:GeneralisedGoursat (undecided): `(-113843540180736*x + 647221602781440*x^4 - 1136096985480192*x^7 + 957998368901088*x^10 - 478999184450544*x^13 + 142012123185024*x^16 - 20225675086920*x^19 + 88`
- `G-s2026-cbrt-second-level-012` Mathilda:GeneralisedGoursat (undecided): `(-97667298956307328795253086622208*x + 20798230922548363130476805494560*x^4 - 1356383060477426800680576120144*x^7 + 42483209437087572250704723408*x^10 - 7867261`
- `G-s2026-cbrt-second-level-013` Mathilda:GeneralisedGoursat (undecided): `(387749091728124228403200*x - 278736733148749253050368*x^4 + 61354951507518538383360*x^7 - 6486073177014991847424*x^10 + 405379573563436990464*x^13 - 1497923620`
- `G-s2026-cbrt-second-level-014` Mathilda:GeneralisedGoursat (undecided): `(55837065975326067331694592 - 34887227964901841537335296*x^3 + 4141491127129914482884608*x^6 - 279240085470072023285760*x^9 + 17454642526313634594816*x^12 - 109`
- `G-s2026-cbrt-second-level-015` Mathilda:GeneralisedGoursat (undecided): `(112036264704*x - 652507512480*x^4 + 1153082418048*x^7 - 978301038942*x^10 + 489150519471*x^13 - 144135302256*x^16 + 20390859765*x^19 - 875283318*x^22)/((-2 + x`
- `G-s2026-cbrt-second-level-016` Mathilda:GeneralisedGoursat (undecided): `(-1389066286457864676466710869488128 + 257240487432694165382883753063936*x^3 - 9051276038179527709356302945280*x^6 + 180850915230036760470445688448*x^9 - 334906`
- `G-s2026-cbrt-second-level-017` Mathilda:GeneralisedGoursat (undecided): `(112036264704*x - 652507512480*x^4 + 1153082418048*x^7 - 978301038942*x^10 + 489150519471*x^13 - 144135302256*x^16 + 20390859765*x^19 - 875283318*x^22)/((-2 + x`
- `G-s2026-cbrt-second-level-018` Mathilda:GeneralisedGoursat (undecided): `(502599561186062688659177472*x - 361470666030756681228484608*x^4 + 79576737175298503958593536*x^7 - 8413391296701038714683392*x^10 + 525836956043814919667712*x^`
- `G-s2026-cbrt-second-level-019` Mathilda:GeneralisedGoursat (undecided): `(-689287396496585313484800 + 430835787040563916701696*x^3 - 51165645299018215981056*x^6 + 3450423509612136235008*x^9 - 215645383115409457152*x^12 + 134782168344`
- `G-s2026-cbrt-second-level-020` Mathilda:GeneralisedGoursat (undecided): `(-7060122942336*x - 39157523688960*x^4 - 68268417682752*x^7 - 57199751793648*x^10 - 28599875896824*x^13 - 8533552210344*x^16 - 1223672615280*x^19 - 55157210487*`
- `G-s2026-cbrt-second-level-021` Mathilda:GeneralisedGoursat (undecided): `(-502429929139783347512279040*x + 360967046349235311913992192*x^4 - 79442544838580696128684032*x^7 + 8396948375471622422790144*x^10 - 524809273466976401424384*x`
- `G-s2026-cbrt-second-level-022` Mathilda:GeneralisedGoursat (undecided): `(48835165678703603379772929763200*x - 10400449147643432842668762352992*x^4 + 678296820834075327937890194160*x^7 - 21245427333385140554392258296*x^10 + 393433839`
- `G-s2026-cbrt-second-level-023` Mathilda:GeneralisedGoursat (undecided): `(-12673236091392 - 63856292405760*x^3 - 61154383352832*x^6 - 33101789173632*x^9 - 16502336888256*x^12 - 8275447293408*x^15 - 3822148959552*x^18 - 997754568840*x`
- `G-s2026-cbrt-second-level-024` Mathilda:GeneralisedGoursat (undecided): `(-218093671657588352090112 + 136321692223786822139904*x^3 - 16189844414585722896384*x^6 + 1091795199471088828416*x^9 - 68234632405202239488*x^12 + 4264824997933`
- `G-s7-cbrt-symmetric-000` Mathilda:GeneralisedGoursat (undecided-nonnumeric): `(-27648 + 169344*x - 220384*x^2 - 369168*x^3 + 1119744*x^4 - 953964*x^5 - 931338*x^6 + 1850202*x^7 - 708588*x^8)/((2*x - 3*x^2)^(2/3)*(-1152*x + 6192*x^2 - 9720`
- `G-s7-cbrt-symmetric-022` Mathilda:GeneralisedGoursat (undecided): `(25007781883412480 + 143720397474365440*x + 347982563769647104*x^2 + 467487673184993280*x^3 + 388498948720148480*x^4 + 211449781407517696*x^5 + 7746740202830668`
- `G-s7-sqrt-klein-000` Mathilda:GeneralisedGoursat (undecided): `(254085 - 1175772*x + 2268644*x^2 - 2465428*x^3 + 1695009*x^4 - 785400*x^5 + 255024*x^6 - 56712*x^7 + 5175*x^8 + 2196*x^9 - 1044*x^10 + 92*x^11 + 35*x^12)/(Sqrt`
- `G-s7-sqrt-klein-005` Mathilda:GeneralisedGoursat (undecided): `(-144 + 288*x + 72*x^2 + 3*x^4)/(Sqrt[12*x + 6*x^2 + x^3]*(2304 + 192*x - 288*x^2 + 16*x^3 + 16*x^4))`

## Verdicts against the corpus labels (0)

A certificate on an `elem` case, or a verified elementary answer on a `nonelem` case (a special-function answer to a `nonelem` case is consistent and not listed). For goursat.py `verified` is its own numerical check.


How each run ends on the `nonelem` cases (a certificate or a special-function answer is a decision; uneval/timeout is not):

| run | CERT | SPECIAL | ELEM | UNEVAL | TIMEOUT | UNVERIFIED | WRONG | ERROR |
|---|---|---|---|---|---|---|---|---|
| Mathilda:Automatic | 123 | 61 | 0 | 0 | 26 | 0 | 0 | 0 |
| Mathilda:GoursatAlgebraic | 0 | 0 | 0 | 210 | 0 | 0 | 0 | 0 |
| Mathilda:GeneralisedGoursat | 131 | 62 | 0 | 2 | 11 | 4 | 0 | 0 |
| goursat.py | 195 | 0 | 0 | 0 | 14 | 0 | 0 | 1 |

GeneralisedGoursat and goursat.py reach opposite verdicts on 0 cases.

Provisional timeouts (at 60 s, under load): Mathilda:Automatic 148, Mathilda:GeneralisedGoursat 32, Mathilda:GoursatAlgebraic 25.


## Worst offenders (times under load: indicative only)

### Mathilda:Automatic: slowest decided

| case | label | s | class | goursat.py s |
|---|---|---|---|---|
| `G-s7-cbrt-second-level-005` | cbrt-second-level #5 | 58.75 | ELEM | 35.53 (ELEM) |
| `G-s7-cbrt-second-level-004` | cbrt-second-level #4 | 57.94 | ELEM | 36.68 (ELEM) |
| `G-s2026-cbrt-second-level-004` | cbrt-second-level #4 | 56.29 | ELEM | 42.61 (ELEM) |
| `G-s2026-cbrt-second-level-024` | cbrt-second-level #24 | 55.98 | ELEM | 25.84 (ELEM) |
| `G-s2026-cbrt-second-level-012` | cbrt-second-level #12 | 53.67 | ELEM | 41.56 (ELEM) |
| `G-s7-cbrt-second-level-012` | cbrt-second-level #12 | 52.58 | ELEM | 36.06 (ELEM) |
| `G-s7-sqrt-mixed-019` | sqrt-mixed #19 | 51.90 | ELEM | 13.59 (ELEM) |
| `G-s11-cbrt-symmetric-018` | cbrt-symmetric #18 | 50.65 | ELEM | 3.51 (ELEM) |
| `G-s2026-cbrt-second-level-001` | cbrt-second-level #1 | 50.64 | ELEM | 45.18 (ELEM) |
| `G-s2026-cbrt-second-level-018` | cbrt-second-level #18 | 49.77 | ELEM | 45.17 (ELEM) |

Mathilda:Automatic: 228 cases goursat.py decides and this run does not (TIMEOUT 134, UNEVAL 94).

Mathilda:Automatic: slower than goursat.py on 199 decided cases (over 1 s); the worst:

- `G-unit-unit-022`: 10.16 s vs 0.03 s (366.8x)
- `G-ex-examples-025`: 10.36 s vs 0.04 s (266.3x)
- `G-welz-welz-016`: 10.27 s vs 0.32 s (32.0x)
- `G-s7-cbrt-symmetric-016`: 26.01 s vs 0.82 s (31.8x)
- `G-welz-welz-015`: 10.07 s vs 0.38 s (26.2x)
- `G-welz-welz-026`: 10.18 s vs 0.43 s (23.8x)
- `G-welz-welz-011`: 10.22 s vs 0.48 s (21.2x)
- `G-unit-chebyshev-018`: 2.19 s vs 0.10 s (21.2x)
- `G-welz-welz-014`: 10.08 s vs 0.48 s (21.0x)
- `G-welz-welz-008`: 10.33 s vs 0.49 s (20.9x)

### Mathilda:GoursatAlgebraic: slowest decided

| case | label | s | class | goursat.py s |
|---|---|---|---|---|
| `G-s7-cbrt-symmetric-003` | cbrt-symmetric #3 | 22.74 | ELEM | 5.69 (ELEM) |
| `G-s7-sqrt-klein-008` | sqrt-klein #8 | 18.67 | ELEM | 2.16 (ELEM) |
| `G-s7-cbrt-symmetric-002` | cbrt-symmetric #2 | 14.12 | ELEM | 5.19 (ELEM) |
| `G-s7-cbrt-symmetric-004` | cbrt-symmetric #4 | 8.22 | ELEM | 4.69 (ELEM) |
| `G-s11-cbrt-symmetric-017` | cbrt-symmetric #17 | 5.14 | ELEM | 10.70 (ELEM) |
| `G-ex-examples-023` | irreducible cubic: 1/(H (t^3+t+1)^(1/3)) | 3.16 | ELEM | 3.10 (ELEM) |
| `G-s7-sqrt-klein-007` | sqrt-klein #7 | 3.10 | ELEM | 2.19 (ELEM) |
| `G-pcbrt-pseudoelliptic-cbrt-050` | pseudoelliptic-cbrt #50: (2/(2*x - 1 - Sqrt[3]*I))/(x^3 - 3* | 2.06 | ELEM | 15.36 (ELEM) |
| `G-s2026-sqrt-mixed-011` | sqrt-mixed #11 | 1.85 | ELEM | 2.82 (ELEM) |
| `G-pcbrt-pseudoelliptic-cbrt-051` | pseudoelliptic-cbrt #51: (2/(2*x - 1 + Sqrt[3]*I))/(x^3 - 3* | 1.72 | ELEM | 16.52 (ELEM) |

Mathilda:GoursatAlgebraic: 850 cases goursat.py decides and this run does not (UNEVAL 829, TIMEOUT 21).

Mathilda:GoursatAlgebraic: slower than goursat.py on 6 decided cases (over 1 s); the worst:

- `G-s7-sqrt-klein-008`: 18.67 s vs 2.16 s (8.7x)
- `G-s7-cbrt-symmetric-003`: 22.74 s vs 5.69 s (4.0x)
- `G-s7-cbrt-symmetric-002`: 14.12 s vs 5.19 s (2.7x)
- `G-s7-cbrt-symmetric-004`: 8.22 s vs 4.69 s (1.8x)
- `G-s7-sqrt-klein-007`: 3.10 s vs 2.19 s (1.4x)
- `G-ex-examples-023`: 3.16 s vs 3.10 s (1.0x)

### Mathilda:GeneralisedGoursat: slowest decided

| case | label | s | class | goursat.py s |
|---|---|---|---|---|
| `G-s11-cbrt-nonelementary-010` | cbrt-nonelementary #10 | 57.78 | CERT | 24.93 (CERT) |
| `G-s2026-cbrt-symmetric-009` | cbrt-symmetric #9 | 55.41 | ELEM | 24.92 (ELEM) |
| `G-s7-cbrt-mixed-011` | cbrt-mixed #11 | 55.13 | ELEM | 59.47 (ELEM) |
| `G-s11-cbrt-nonelementary-008` | cbrt-nonelementary #8 | 53.78 | CERT | 64.35 (CERT) |
| `G-s7-cbrt-mixed-013` | cbrt-mixed #13 | 52.50 | ELEM | 73.60 (ELEM) |
| `G-s7-cbrt-mixed-020` | cbrt-mixed #20 | 38.82 | ELEM | 24.83 (ELEM) |
| `G-s11-cbrt-mixed-013` | cbrt-mixed #13 | 36.54 | ELEM | 42.35 (ELEM) |
| `G-s2026-cbrt-mixed-007` | cbrt-mixed #7 | 34.04 | ELEM | 1.27 (ELEM) |
| `G-s11-cbrt-mixed-010` | cbrt-mixed #10 | 33.51 | ELEM | 36.89 (ELEM) |
| `G-s11-cbrt-nonelementary-012` | cbrt-nonelementary #12 | 33.28 | CERT | 58.87 (CERT) |

Mathilda:GeneralisedGoursat: 143 cases goursat.py decides and this run does not (UNVERIFIED 107, TIMEOUT 21, UNEVAL 14, WRONG 1).

- `G-unit-cubic-field-000` (irreducible resolvent: g - F0(g), g=1/(t-2), R=t^4+t+1): UNEVAL 0.2 s; goursat.py ELEM 0.6 s
- `G-ex-examples-024` (irreducible resolvent: (F - F0), F=1/(t-2), R=t^4+t+1): UNEVAL 0.2 s; goursat.py ELEM 0.8 s
- `G-s2026-cbrt-symmetric-000` (cbrt-symmetric #0): UNVERIFIED 30.7 s; goursat.py ELEM 21.4 s
- `G-s2026-cbrt-symmetric-006` (cbrt-symmetric #6): UNVERIFIED 37.1 s; goursat.py ELEM 24.7 s
- `G-s2026-cbrt-symmetric-017` (cbrt-symmetric #17): TIMEOUT 91.4 s; goursat.py ELEM 45.9 s
- `G-s2026-cbrt-mixed-012` (cbrt-mixed #12): UNVERIFIED 31.2 s; goursat.py ELEM 37.4 s
- `G-s2026-sqrt-klein-011` (sqrt-klein #11): UNVERIFIED 10.5 s; goursat.py ELEM 37.3 s
- `G-s2026-sqrt-klein-019` (sqrt-klein #19): UNEVAL 0.3 s; goursat.py ELEM 13.9 s
- `G-s2026-sqrt-mixed-005` (sqrt-mixed #5): UNVERIFIED 5.1 s; goursat.py ELEM 18.1 s
- `G-s2026-sqrt-mixed-008` (sqrt-mixed #8): UNEVAL 6.6 s; goursat.py ELEM 104.4 s
- `G-s2026-sqrt-mixed-017` (sqrt-mixed #17): UNVERIFIED 17.3 s; goursat.py ELEM 19.6 s
- `G-s2026-sqrt-nonelementary-006` (sqrt-nonelementary #6): UNVERIFIED 7.1 s; goursat.py CERT 30.8 s
- `G-s2026-sqrt-nonelementary-021` (sqrt-nonelementary #21): UNEVAL 0.4 s; goursat.py CERT 14.6 s
- `G-s2026-sqrt-order3-015` (sqrt-order3 #15): UNVERIFIED 7.1 s; goursat.py ELEM 47.8 s
- `G-s2026-sqrt-order3-021` (sqrt-order3 #21): UNEVAL 10.2 s; goursat.py ELEM 21.5 s
- `G-s2026-cbrt-second-level-000` (cbrt-second-level #0): UNVERIFIED 9.5 s; goursat.py ELEM 45.6 s
- `G-s2026-cbrt-second-level-001` (cbrt-second-level #1): UNVERIFIED 12.1 s; goursat.py ELEM 45.2 s
- `G-s2026-cbrt-second-level-002` (cbrt-second-level #2): UNVERIFIED 13.1 s; goursat.py ELEM 45.3 s
- `G-s2026-cbrt-second-level-003` (cbrt-second-level #3): UNVERIFIED 10.6 s; goursat.py ELEM 45.2 s
- `G-s2026-cbrt-second-level-004` (cbrt-second-level #4): UNVERIFIED 22.7 s; goursat.py ELEM 42.6 s
- `G-s2026-cbrt-second-level-005` (cbrt-second-level #5): UNVERIFIED 15.5 s; goursat.py ELEM 39.7 s
- `G-s2026-cbrt-second-level-006` (cbrt-second-level #6): UNVERIFIED 14.1 s; goursat.py ELEM 39.3 s
- `G-s2026-cbrt-second-level-007` (cbrt-second-level #7): UNVERIFIED 14.0 s; goursat.py ELEM 39.5 s
- `G-s2026-cbrt-second-level-008` (cbrt-second-level #8): UNVERIFIED 9.4 s; goursat.py ELEM 39.3 s
- `G-s2026-cbrt-second-level-009` (cbrt-second-level #9): UNVERIFIED 14.1 s; goursat.py ELEM 41.2 s
- `G-s2026-cbrt-second-level-010` (cbrt-second-level #10): UNVERIFIED 8.9 s; goursat.py ELEM 39.3 s
- `G-s2026-cbrt-second-level-011` (cbrt-second-level #11): UNVERIFIED 10.4 s; goursat.py ELEM 41.6 s
- `G-s2026-cbrt-second-level-012` (cbrt-second-level #12): UNVERIFIED 13.4 s; goursat.py ELEM 41.6 s
- `G-s2026-cbrt-second-level-013` (cbrt-second-level #13): UNVERIFIED 13.5 s; goursat.py ELEM 41.2 s
- `G-s2026-cbrt-second-level-014` (cbrt-second-level #14): UNVERIFIED 11.2 s; goursat.py ELEM 41.4 s
- `G-s2026-cbrt-second-level-015` (cbrt-second-level #15): UNVERIFIED 10.1 s; goursat.py ELEM 43.3 s
- `G-s2026-cbrt-second-level-016` (cbrt-second-level #16): UNVERIFIED 11.8 s; goursat.py ELEM 41.8 s
- `G-s2026-cbrt-second-level-017` (cbrt-second-level #17): UNVERIFIED 9.4 s; goursat.py ELEM 44.9 s
- `G-s2026-cbrt-second-level-018` (cbrt-second-level #18): UNVERIFIED 11.3 s; goursat.py ELEM 45.2 s
- `G-s2026-cbrt-second-level-019` (cbrt-second-level #19): UNVERIFIED 14.2 s; goursat.py ELEM 44.2 s
- `G-s2026-cbrt-second-level-020` (cbrt-second-level #20): UNVERIFIED 18.1 s; goursat.py ELEM 45.2 s
- `G-s2026-cbrt-second-level-021` (cbrt-second-level #21): UNVERIFIED 11.9 s; goursat.py ELEM 43.9 s
- `G-s2026-cbrt-second-level-022` (cbrt-second-level #22): UNVERIFIED 14.3 s; goursat.py ELEM 38.5 s
- `G-s2026-cbrt-second-level-023` (cbrt-second-level #23): UNVERIFIED 13.1 s; goursat.py ELEM 27.7 s
- `G-s2026-cbrt-second-level-024` (cbrt-second-level #24): UNVERIFIED 13.8 s; goursat.py ELEM 25.8 s
- `G-s7-cbrt-symmetric-000` (cbrt-symmetric #0): UNVERIFIED 4.7 s; goursat.py ELEM 17.3 s
- `G-s7-cbrt-symmetric-015` (cbrt-symmetric #15): TIMEOUT 60.2 s; goursat.py ELEM 59.5 s
- `G-s7-cbrt-symmetric-018` (cbrt-symmetric #18): TIMEOUT 60.1 s; goursat.py ELEM 54.0 s
- `G-s7-cbrt-symmetric-022` (cbrt-symmetric #22): UNVERIFIED 27.2 s; goursat.py ELEM 50.5 s
- `G-s7-cbrt-nonelementary-024` (cbrt-nonelementary #24): TIMEOUT 60.2 s; goursat.py CERT 85.2 s
- `G-s7-sqrt-klein-005` (sqrt-klein #5): UNVERIFIED 6.0 s; goursat.py ELEM 35.4 s
- `G-s7-sqrt-klein-010` (sqrt-klein #10): UNVERIFIED 6.5 s; goursat.py ELEM 45.2 s
- `G-s7-sqrt-klein-013` (sqrt-klein #13): UNEVAL 0.3 s; goursat.py ELEM 15.1 s
- `G-s7-sqrt-klein-014` (sqrt-klein #14): UNVERIFIED 7.0 s; goursat.py ELEM 38.8 s
- `G-s7-sqrt-klein-017` (sqrt-klein #17): UNEVAL 6.9 s; goursat.py ELEM 91.7 s
- `G-s7-sqrt-klein-024` (sqrt-klein #24): UNVERIFIED 7.3 s; goursat.py ELEM 34.7 s
- `G-s7-sqrt-mixed-013` (sqrt-mixed #13): UNVERIFIED 6.3 s; goursat.py ELEM 33.0 s
- `G-s7-sqrt-mixed-019` (sqrt-mixed #19): UNVERIFIED 10.8 s; goursat.py ELEM 13.6 s
- `G-s7-sqrt-nonelementary-001` (sqrt-nonelementary #1): UNVERIFIED 7.7 s; goursat.py CERT 27.2 s
- `G-s7-sqrt-nonelementary-012` (sqrt-nonelementary #12): TIMEOUT 60.4 s; goursat.py CERT 31.2 s
- `G-s7-sqrt-nonelementary-013` (sqrt-nonelementary #13): TIMEOUT 60.8 s; goursat.py CERT 14.7 s
- `G-s7-sqrt-nonelementary-019` (sqrt-nonelementary #19): UNVERIFIED 5.5 s; goursat.py CERT 16.0 s
- `G-s7-sqrt-order3-004` (sqrt-order3 #4): UNVERIFIED 6.4 s; goursat.py ELEM 29.7 s
- `G-s7-sqrt-order3-016` (sqrt-order3 #16): UNEVAL 0.4 s; goursat.py ELEM 12.0 s
- `G-s7-sqrt-order3-018` (sqrt-order3 #18): UNVERIFIED 5.2 s; goursat.py ELEM 32.6 s
- `G-s7-sqrt-order3-021` (sqrt-order3 #21): UNVERIFIED 7.4 s; goursat.py ELEM 37.0 s
- `G-s7-sqrt-order3-023` (sqrt-order3 #23): UNEVAL 0.6 s; goursat.py ELEM 16.0 s
- `G-s7-cbrt-second-level-000` (cbrt-second-level #0): UNVERIFIED 19.1 s; goursat.py ELEM 35.9 s
- `G-s7-cbrt-second-level-001` (cbrt-second-level #1): UNVERIFIED 13.8 s; goursat.py ELEM 37.0 s
- `G-s7-cbrt-second-level-002` (cbrt-second-level #2): UNVERIFIED 14.6 s; goursat.py ELEM 35.7 s
- `G-s7-cbrt-second-level-003` (cbrt-second-level #3): UNVERIFIED 22.4 s; goursat.py ELEM 35.9 s
- `G-s7-cbrt-second-level-004` (cbrt-second-level #4): UNVERIFIED 33.7 s; goursat.py ELEM 36.7 s
- `G-s7-cbrt-second-level-005` (cbrt-second-level #5): UNVERIFIED 18.6 s; goursat.py ELEM 35.5 s
- `G-s7-cbrt-second-level-006` (cbrt-second-level #6): UNVERIFIED 19.4 s; goursat.py ELEM 36.5 s
- `G-s7-cbrt-second-level-007` (cbrt-second-level #7): UNVERIFIED 31.9 s; goursat.py ELEM 36.1 s
- `G-s7-cbrt-second-level-008` (cbrt-second-level #8): UNVERIFIED 18.2 s; goursat.py ELEM 35.8 s
- `G-s7-cbrt-second-level-009` (cbrt-second-level #9): UNVERIFIED 13.8 s; goursat.py ELEM 36.5 s
- `G-s7-cbrt-second-level-010` (cbrt-second-level #10): UNVERIFIED 15.7 s; goursat.py ELEM 35.4 s
- `G-s7-cbrt-second-level-011` (cbrt-second-level #11): UNVERIFIED 28.8 s; goursat.py ELEM 35.4 s
- `G-s7-cbrt-second-level-012` (cbrt-second-level #12): UNVERIFIED 19.4 s; goursat.py ELEM 36.1 s
- `G-s7-cbrt-second-level-013` (cbrt-second-level #13): UNVERIFIED 17.1 s; goursat.py ELEM 35.4 s
- `G-s7-cbrt-second-level-014` (cbrt-second-level #14): UNVERIFIED 23.2 s; goursat.py ELEM 35.8 s
- `G-s7-cbrt-second-level-015` (cbrt-second-level #15): UNVERIFIED 29.3 s; goursat.py ELEM 36.1 s
- `G-s7-cbrt-second-level-016` (cbrt-second-level #16): UNVERIFIED 36.0 s; goursat.py ELEM 35.5 s
- `G-s7-cbrt-second-level-017` (cbrt-second-level #17): UNVERIFIED 19.8 s; goursat.py ELEM 34.2 s

Mathilda:GeneralisedGoursat: slower than goursat.py on 65 decided cases (over 1 s); the worst:

- `G-s2026-cbrt-mixed-007`: 34.04 s vs 1.27 s (26.9x)
- `G-lit-literature-032`: 2.52 s vs 0.14 s (18.0x)
- `G-unit-chebyshev-021`: 2.42 s vs 0.15 s (16.4x)
- `G-unit-chebyshev-092`: 2.09 s vs 0.13 s (16.1x)
- `G-unit-chebyshev-015`: 2.09 s vs 0.14 s (15.4x)
- `G-unit-chebyshev-100`: 2.38 s vs 0.16 s (15.2x)
- `G-unit-chebyshev-089`: 2.11 s vs 0.15 s (14.5x)
- `G-unit-chebyshev-081`: 2.11 s vs 0.16 s (13.5x)
- `G-unit-chebyshev-093`: 1.42 s vs 0.11 s (13.1x)
- `G-unit-chebyshev-024`: 2.43 s vs 0.19 s (13.0x)

GeneralisedGoursat decides 34 cases goursat.py does not (ERROR 17, TIMEOUT 16, WRONG 1).

