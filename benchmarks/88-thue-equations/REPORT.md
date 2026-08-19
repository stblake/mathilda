# Benchmark 88 — Thue-equation stress tests

Mathilda's `Solve[F(x,y)==m && Element[{x,y},Integers], {x,y}, Integers]` vs **PARI/GP `thue()`** (gold-standard oracle) over an adversarial corpus. A `WRONG` verdict (finite sets differ) is a completeness bug.

| verdict | n | meaning |
|---|---|---|
| CORRECT | 98 | matches PARI |
| DECLINE | 6 | unevaluated (honest gap); PARI solved |

## Bottlenecks — slowest CORRECT solves

Both times are COMPUTE time (Mathilda `Timing[]`, PARI `gettime()`), excluding each process's ~30–50 ms startup, so the ratio is fair.

| label | form == m | Mathilda | PARI | ratio | #sol |
|---|---|---:|---:|---:|---:|
| adv-big-coef-2 | `131*y^3 - 211*x*y^2 + 97*x^2*y + x^3 == 1` | 383.6ms | 11ms | 34.9x | 1 |
| binom-sextic-d3 | `-3*y^6 + x^6 == 1` | 249.6ms | 11ms | 22.7x | 2 |
| binom-cubic-d2-m100 | `-2*y^3 + x^3 == 100` | 244.6ms | 8ms | 30.6x | 0 |
| septic-x7-2 | `-2*y^7 + x^7 == 1` | 180.4ms | 20ms | 9.0x | 2 |
| binom-cubic-d41-m1 | `-41*y^3 + x^3 == -1` | 146.0ms | 8ms | 18.2x | 1 |
| binom-cubic-d41-p1 | `-41*y^3 + x^3 == 1` | 142.2ms | 8ms | 17.8x | 1 |
| adv-big-coef-cubic | `-97*y^3 + x^3 == 1` | 132.7ms | 9ms | 14.7x | 1 |
| binom-sextic-d2 | `-2*y^6 + x^6 == 1` | 132.3ms | 16ms | 8.3x | 2 |
| quartic-x4-x-1 | `-y^4 - x*y^3 + x^4 == 1` | 86.1ms | 12ms | 7.2x | 4 |
| nosol-cyclic | `y^3 - 3*x*y^2 + x^3 == 4` | 75.0ms | 8ms | 9.4x | 0 |
| binom-quintic-d3 | `-3*y^5 + x^5 == 1` | 48.5ms | 8ms | 6.1x | 1 |
| quartic-biquad-1 | `y^4 - 4*x^2*y^2 + x^4 == 1` | 45.4ms | 11ms | 4.1x | 12 |
| cyclic-cubic-m2 | `y^3 - 3*x*y^2 + x^3 == 2` | 42.3ms | 8ms | 5.3x | 0 |
| binom-quartic-d6-p1 | `-6*y^4 + x^4 == 1` | 38.4ms | 8ms | 4.8x | 2 |
| binom-quartic-d6-m1 | `-6*y^4 + x^4 == -1` | 38.1ms | 8ms | 4.8x | 0 |

## Coverage by family

| family | CORRECT | DECLINE | WRONG | CRASH | TIMEOUT |
|---|---:|---:|---:|---:|---:|
| adversarial-precision | 5 | 0 | 0 | 0 | 0 |
| binomial-cubic | 32 | 0 | 0 | 0 | 0 |
| binomial-cubic-large | 2 | 0 | 0 | 0 | 0 |
| binomial-quartic | 15 | 2 | 0 | 0 | 0 |
| general-cubic | 4 | 1 | 0 | 0 | 0 |
| general-quartic | 4 | 2 | 0 | 0 | 0 |
| high-m | 9 | 0 | 0 | 0 | 0 |
| large-solutions | 1 | 0 | 0 | 0 | 0 |
| many-solutions | 1 | 0 | 0 | 0 | 0 |
| no-solution | 2 | 0 | 0 | 0 | 0 |
| non-monogenic | 3 | 0 | 0 | 0 | 0 |
| quintic | 3 | 1 | 0 | 0 | 0 |
| reducible | 3 | 0 | 0 | 0 | 0 |
| septic | 1 | 0 | 0 | 0 | 0 |
| sextic | 2 | 0 | 0 | 0 | 0 |
| simplest-cubic | 5 | 0 | 0 | 0 | 0 |
| simplest-cubic-family | 6 | 0 | 0 | 0 | 0 |

## All cases

| label | form == m | verdict | M#/P# | M time | note |
|---|---|---|---|---:|---|
| binom-cubic-d2-p1 | `-2*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.2ms | x^3 - 2 y^3 = 1 |
| binom-cubic-d2-m1 | `-2*y^3 + x^3 == -1` | CORRECT | 2/2 | 6.0ms | x^3 - 2 y^3 = -1 |
| binom-cubic-d3-p1 | `-3*y^3 + x^3 == 1` | CORRECT | 1/1 | 2.4ms | x^3 - 3 y^3 = 1 |
| binom-cubic-d3-m1 | `-3*y^3 + x^3 == -1` | CORRECT | 1/1 | 2.4ms | x^3 - 3 y^3 = -1 |
| binom-cubic-d5-p1 | `-5*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.9ms | x^3 - 5 y^3 = 1 |
| binom-cubic-d5-m1 | `-5*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.5ms | x^3 - 5 y^3 = -1 |
| binom-cubic-d6-p1 | `-6*y^3 + x^3 == 1` | CORRECT | 1/1 | 7.5ms | x^3 - 6 y^3 = 1 |
| binom-cubic-d6-m1 | `-6*y^3 + x^3 == -1` | CORRECT | 1/1 | 8.1ms | x^3 - 6 y^3 = -1 |
| binom-cubic-d7-p1 | `-7*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.0ms | x^3 - 7 y^3 = 1 |
| binom-cubic-d7-m1 | `-7*y^3 + x^3 == -1` | CORRECT | 2/2 | 5.8ms | x^3 - 7 y^3 = -1 |
| binom-cubic-d10-p1 | `-10*y^3 + x^3 == 1` | CORRECT | 1/1 | 7.1ms | x^3 - 10 y^3 = 1 |
| binom-cubic-d10-m1 | `-10*y^3 + x^3 == -1` | CORRECT | 1/1 | 7.0ms | x^3 - 10 y^3 = -1 |
| binom-cubic-d11-p1 | `-11*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.7ms | x^3 - 11 y^3 = 1 |
| binom-cubic-d11-m1 | `-11*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.7ms | x^3 - 11 y^3 = -1 |
| binom-cubic-d12-p1 | `-12*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.5ms | x^3 - 12 y^3 = 1 |
| binom-cubic-d12-m1 | `-12*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.6ms | x^3 - 12 y^3 = -1 |
| binom-cubic-d13-p1 | `-13*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.5ms | x^3 - 13 y^3 = 1 |
| binom-cubic-d13-m1 | `-13*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.3ms | x^3 - 13 y^3 = -1 |
| binom-cubic-d15-p1 | `-15*y^3 + x^3 == 1` | CORRECT | 1/1 | 27.9ms | x^3 - 15 y^3 = 1 |
| binom-cubic-d15-m1 | `-15*y^3 + x^3 == -1` | CORRECT | 1/1 | 32.1ms | x^3 - 15 y^3 = -1 |
| binom-cubic-d20-p1 | `-20*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.3ms | x^3 - 20 y^3 = 1 |
| binom-cubic-d20-m1 | `-20*y^3 + x^3 == -1` | CORRECT | 2/2 | 6.2ms | x^3 - 20 y^3 = -1 |
| binom-cubic-d30-p1 | `-30*y^3 + x^3 == 1` | CORRECT | 1/1 | 12.0ms | x^3 - 30 y^3 = 1 |
| binom-cubic-d30-m1 | `-30*y^3 + x^3 == -1` | CORRECT | 1/1 | 11.8ms | x^3 - 30 y^3 = -1 |
| binom-cubic-d41-p1 | `-41*y^3 + x^3 == 1` | CORRECT | 1/1 | 142.2ms | x^3 - 41 y^3 = 1 |
| binom-cubic-d41-m1 | `-41*y^3 + x^3 == -1` | CORRECT | 1/1 | 146.0ms | x^3 - 41 y^3 = -1 |
| binom-cubic-d42-p1 | `-42*y^3 + x^3 == 1` | CORRECT | 1/1 | 31.6ms | x^3 - 42 y^3 = 1 |
| binom-cubic-d42-m1 | `-42*y^3 + x^3 == -1` | CORRECT | 1/1 | 31.0ms | x^3 - 42 y^3 = -1 |
| binom-cubic-d43-p1 | `-43*y^3 + x^3 == 1` | CORRECT | 2/2 | 9.4ms | x^3 - 43 y^3 = 1 |
| binom-cubic-d43-m1 | `-43*y^3 + x^3 == -1` | CORRECT | 2/2 | 8.7ms | x^3 - 43 y^3 = -1 |
| binom-cubic-d17-p1 | `-17*y^3 + x^3 == 1` | CORRECT | 2/2 | 28.8ms | x^3-17y^3=1 -> (18,7); Q(cbrt17) non-monogenic (index 3) |
| binom-cubic-d20-p1 | `-20*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.1ms | x^3-20y^3=1 |
| binom-cubic-plus-d2 | `2*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.9ms | x^3+2y^3=1 |
| binom-cubic-plus-d7 | `7*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.0ms | x^3+7y^3=1 |
| binom-cubic-d2-m2 | `-2*y^3 + x^3 == 2` | CORRECT | 1/1 | 10.7ms | x^3-2y^3=2 |
| binom-cubic-d2-m3 | `-2*y^3 + x^3 == 3` | CORRECT | 2/2 | 12.0ms | x^3-2y^3=3 |
| binom-cubic-d2-m4 | `-2*y^3 + x^3 == 4` | CORRECT | 0/0 | 12.0ms | x^3-2y^3=4 |
| binom-cubic-d2-m5 | `-2*y^3 + x^3 == 5` | CORRECT | 0/0 | 11.3ms | x^3-2y^3=5 |
| binom-cubic-d2-m9 | `-2*y^3 + x^3 == 9` | CORRECT | 0/0 | 14.7ms | x^3-2y^3=9 |
| binom-cubic-d2-m10 | `-2*y^3 + x^3 == 10` | CORRECT | 2/2 | 13.6ms | x^3-2y^3=10 |
| binom-cubic-d2-m100 | `-2*y^3 + x^3 == 100` | CORRECT | 0/0 | 244.6ms | x^3-2y^3=100 |
| binom-cubic-d2-m73 | `-2*y^3 + x^3 == 73` | CORRECT | 0/0 | 11.6ms | x^3-2y^3=73 (has solutions) |
| cyclic-cubic-m2 | `y^3 - 3*x*y^2 + x^3 == 2` | CORRECT | 0/0 | 42.3ms | x^3-3xy^2+y^3=2 |
| cyclic-cond9-p1 | `y^3 - 3*x*y^2 + x^3 == 1` | CORRECT | 6/6 | 8.4ms | cond 9, x^3-3xy^2+y^3=1 |
| cyclic-cond9-m1 | `y^3 - 3*x*y^2 + x^3 == -1` | CORRECT | 6/6 | 8.2ms | cond 9, =-1 |
| thomas-cond7-p1 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | CORRECT | 9/9 | 38.1ms | x^3+x^2y-2xy^2-y^3=1 (cond 7, 9 solutions) |
| thomas-cond7-m1 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == -1` | CORRECT | 9/9 | 37.6ms | cond 7, =-1 |
| thomas-t1 | `-y^3 - 3*x*y^2 - x^2*y + x^3 == 1` | CORRECT | 3/3 | 0.8ms | x^3-1x^2y-3xy^2-y^3=1 |
| thomas-t2 | `-y^3 - 4*x*y^2 - 2*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-2x^2y-4xy^2-y^3=1 |
| thomas-t3 | `-y^3 - 5*x*y^2 - 3*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-3x^2y-5xy^2-y^3=1 |
| thomas-t4 | `-y^3 - 6*x*y^2 - 4*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-4x^2y-6xy^2-y^3=1 |
| thomas-t5 | `-y^3 - 7*x*y^2 - 5*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-5x^2y-7xy^2-y^3=1 |
| thomas-t10 | `-y^3 - 12*x*y^2 - 10*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-10x^2y-12xy^2-y^3=1 |
| shanks-simplest-n3 | `y^3 - 3*x*y^2 - 6*x^2*y + x^3 == 1` | CORRECT | 4/4 | 11.8ms | Shanks x^3-3x^2y-6xy^2-... variant |
| cubic-irr-1 | `y^3 + x*y^2 - 2*x^2*y + x^3 == 1` | CORRECT | 4/4 | 2.9ms | x^3-2x^2y+xy^2+y^3 |
| cubic-irr-2 | `-y^3 + 3*x*y^2 - 3*x^2*y + x^3 == 1` | DECLINE | None/None | 1.1ms | totally real-ish |
| cubic-irr-3 | `y^3 - x^2*y + x^3 == 1` | CORRECT | 5/5 | 2.8ms | x^3-xy^2+y^3? small disc |
| cubic-neg-disc | `y^3 + x*y^2 + x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.5ms | x^3+x^2y+xy^2+y^3 (reducible? check) |
| cubic-big-coef | `-5*y^3 + 11*x*y^2 - 7*x^2*y + x^3 == 1` | CORRECT | 1/1 | 0.5ms | bigger coeffs |
| nosol-d2-m4 | `-2*y^3 + x^3 == 4` | CORRECT | 0/0 | 12.3ms | x^3-2y^3=4 (|m|!=1 -> decline; PARI {}) |
| nosol-cyclic | `y^3 - 3*x*y^2 + x^3 == 4` | CORRECT | 0/0 | 75.0ms | cyclic, =4 |
| binom-quartic-d2-p1 | `-2*y^4 + x^4 == 1` | CORRECT | 2/2 | 6.0ms | x^4 - 2 y^4 = 1 |
| binom-quartic-d2-m1 | `-2*y^4 + x^4 == -1` | CORRECT | 4/4 | 6.2ms | x^4 - 2 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d3-p1 | `-3*y^4 + x^4 == 1` | CORRECT | 2/2 | 9.3ms | x^4 - 3 y^4 = 1 |
| binom-quartic-d3-m1 | `-3*y^4 + x^4 == -1` | CORRECT | 0/0 | 9.5ms | x^4 - 3 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d5-p1 | `-5*y^4 + x^4 == 1` | CORRECT | 6/6 | 14.3ms | x^4 - 5 y^4 = 1 |
| binom-quartic-d5-m1 | `-5*y^4 + x^4 == -1` | CORRECT | 0/0 | 14.1ms | x^4 - 5 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d6-p1 | `-6*y^4 + x^4 == 1` | CORRECT | 2/2 | 38.4ms | x^4 - 6 y^4 = 1 |
| binom-quartic-d6-m1 | `-6*y^4 + x^4 == -1` | CORRECT | 0/0 | 38.1ms | x^4 - 6 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d7-p1 | `-7*y^4 + x^4 == 1` | CORRECT | 2/2 | 19.3ms | x^4 - 7 y^4 = 1 |
| binom-quartic-d7-m1 | `-7*y^4 + x^4 == -1` | CORRECT | 0/0 | 18.3ms | x^4 - 7 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d8-p1 | `-8*y^4 + x^4 == 1` | CORRECT | 2/2 | 6.7ms | x^4 - 8 y^4 = 1 |
| binom-quartic-d8-m1 | `-8*y^4 + x^4 == -1` | CORRECT | 0/0 | 6.5ms | x^4 - 8 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d10-p1 | `-10*y^4 + x^4 == 1` | DECLINE | None/2 | 64.4ms | x^4 - 10 y^4 = 1 |
| binom-quartic-d10-m1 | `-10*y^4 + x^4 == -1` | DECLINE | None/0 | 64.6ms | x^4 - 10 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d17-p1 | `-17*y^4 + x^4 == 1` | CORRECT | 2/2 | 10.6ms | x^4 - 17 y^4 = 1 |
| binom-quartic-d17-m1 | `-17*y^4 + x^4 == -1` | CORRECT | 4/4 | 11.5ms | x^4 - 17 y^4 = -1 (subfield unit -> case iii) |
| ljunggren-x4-2y4 | `-2*y^4 + x^4 == -1` | CORRECT | 4/4 | 6.3ms | Ljunggren-type x^4-2y^4=-1 |
| quartic-biquad-1 | `y^4 - 4*x^2*y^2 + x^4 == 1` | CORRECT | 12/12 | 45.4ms | x^4-4x^2y^2+y^4 (biquadratic, quadratic subfield) |
| quartic-tr-1 | `y^4 - 4*x*y^3 + 6*x^2*y^2 - 4*x^3*y + x^4 == 1` | DECLINE | None/None | 1.4ms | (x-y)^4+... totally real? |
| quartic-mixed-1 | `-y^4 + 2*x^3*y + x^4 == 1` | CORRECT | 2/2 | 33.2ms | x^4+2x^3y-y^4 |
| quartic-cyclotomic | `y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4 == 1` | DECLINE | None/6 | 2.0ms | x^4+x^3y+x^2y^2+xy^3+y^4 (cyclotomic, reducible over cyclo) |
| quartic-x4-x-1 | `-y^4 - x*y^3 + x^4 == 1` | CORRECT | 4/4 | 86.1ms | x^4-xy^3-y^4 |
| quartic-totallyimag | `y^4 + x^2*y^2 + x^4 == 1` | CORRECT | 4/4 | 0.5ms | x^4+x^2y^2+y^4 (may be reducible) |
| binom-quintic-d2 | `-2*y^5 + x^5 == 1` | CORRECT | 2/2 | 17.8ms | x^5-2y^5=1 |
| binom-quintic-d3 | `-3*y^5 + x^5 == 1` | CORRECT | 1/1 | 48.5ms | x^5-3y^5=1 |
| binom-quintic-d5 | `-5*y^5 + x^5 == 1` | DECLINE | None/1 | 148.7ms | x^5-5y^5=1 |
| quintic-general | `y^5 - x*y^4 - x^4*y + x^5 == 1` | CORRECT | 2/2 | 0.6ms | general quintic form |
| binom-sextic-d2 | `-2*y^6 + x^6 == 1` | CORRECT | 2/2 | 132.3ms | x^6-2y^6=1 |
| binom-sextic-d3 | `-3*y^6 + x^6 == 1` | CORRECT | 2/2 | 249.6ms | x^6-3y^6=1 |
| septic-x7-2 | `-2*y^7 + x^7 == 1` | CORRECT | 2/2 | 180.4ms | x^7-2y^7=1 |
| nonmono-d17 | `-17*y^3 + x^3 == 1` | CORRECT | 2/2 | 26.7ms | Q(cbrt17) index 3 at p=3 |
| nonmono-d19 | `-19*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.3ms | x^3-19y^3=1 |
| nonmono-cubic-idx | `-8*y^3 - 2*x*y^2 - x^2*y + x^3 == 1` | CORRECT | 1/1 | 28.6ms | Dedekind cubic t^3-t^2-2t-8 (index 2 at 2) |
| reducible-x3-y3 | `-y^3 + x^3 == 7` | CORRECT | 2/2 | 0.4ms | x^3 - y^3 = 7 -> (x-y)(x^2+xy+y^2); reducible |
| reducible-x4-y4 | `-y^4 + x^4 == 15` | CORRECT | 4/4 | 0.7ms | x^4 - y^4 = 15 reducible |
| reducible-biquad | `4*y^4 - 5*x^2*y^2 + x^4 == 3` | CORRECT | 0/0 | 0.5ms | x^4-5x^2y^2+4y^4=(x^2-y^2)(x^2-4y^2) |
| adv-close-roots-1 | `y^4 - 7*x*y^3 + 14*x^2*y^2 - 7*x^3*y + x^4 == 1` | CORRECT | 4/4 | 0.8ms | near-symmetric quartic, clustered roots |
| adv-big-coef-cubic | `-97*y^3 + x^3 == 1` | CORRECT | 1/1 | 132.7ms | x^3-97y^3=1 |
| adv-big-coef-2 | `131*y^3 - 211*x*y^2 + 97*x^2*y + x^3 == 1` | CORRECT | 1/1 | 383.6ms | large mixed coeffs |
| adv-scaled-1 | `-999*y^3 + x^3 == 1` | CORRECT | 2/2 | 23.6ms | x^3-999y^3=1 (999=3^3*37) |
| adv-large-reg | `-y^3 - x*y^2 - x^2*y + x^3 == 1` | CORRECT | 4/4 | 6.4ms | x^3-x^2y-xy^2-y^3=1 (tribonacci field, larger regulator) |
| many-sol-cond7 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | CORRECT | 9/9 | 36.8ms | Thomas cond7 (9 solutions) |
| large-sol-1 | `-3*y^3 + x^3 == 1` | CORRECT | 1/1 | 2.2ms | x^3-3y^3=1 |
