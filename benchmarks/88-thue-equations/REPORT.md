# Benchmark 88 — Thue-equation stress tests

Mathilda's `Solve[F(x,y)==m && Element[{x,y},Integers], {x,y}, Integers]` vs **PARI/GP `thue()`** (gold-standard oracle) over an adversarial corpus. A `WRONG` verdict (finite sets differ) is a completeness bug.

| verdict | n | meaning |
|---|---|---|
| CORRECT | 93 | matches PARI |
| DECLINE | 11 | unevaluated (honest gap); PARI solved |

## Bottlenecks — slowest CORRECT solves

Both times are COMPUTE time (Mathilda `Timing[]`, PARI `gettime()`), excluding each process's ~30–50 ms startup, so the ratio is fair.

| label | form == m | Mathilda | PARI | ratio | #sol |
|---|---|---:|---:|---:|---:|
| adv-big-coef-2 | `131*y^3 - 211*x*y^2 + 97*x^2*y + x^3 == 1` | 379.2ms | 11ms | 34.5x | 1 |
| binom-sextic-d3 | `-3*y^6 + x^6 == 1` | 242.2ms | 11ms | 22.0x | 2 |
| binom-cubic-d2-m100 | `-2*y^3 + x^3 == 100` | 234.5ms | 7ms | 33.5x | 0 |
| septic-x7-2 | `-2*y^7 + x^7 == 1` | 176.6ms | 22ms | 8.0x | 2 |
| binom-cubic-d41-m1 | `-41*y^3 + x^3 == -1` | 138.8ms | 8ms | 17.4x | 1 |
| binom-cubic-d41-p1 | `-41*y^3 + x^3 == 1` | 137.9ms | 8ms | 17.2x | 1 |
| binom-sextic-d2 | `-2*y^6 + x^6 == 1` | 128.4ms | 17ms | 7.6x | 2 |
| adv-big-coef-cubic | `-97*y^3 + x^3 == 1` | 121.5ms | 9ms | 13.5x | 1 |
| quartic-x4-x-1 | `-y^4 - x*y^3 + x^4 == 1` | 87.0ms | 13ms | 6.7x | 4 |
| binom-quintic-d3 | `-3*y^5 + x^5 == 1` | 50.1ms | 9ms | 5.6x | 1 |
| quartic-biquad-1 | `y^4 - 4*x^2*y^2 + x^4 == 1` | 44.7ms | 12ms | 3.7x | 12 |
| binom-quartic-d6-m1 | `-6*y^4 + x^4 == -1` | 40.1ms | 8ms | 5.0x | 0 |
| binom-quartic-d6-p1 | `-6*y^4 + x^4 == 1` | 39.8ms | 7ms | 5.7x | 2 |
| many-sol-cond7 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | 37.0ms | 10ms | 3.7x | 9 |
| thomas-cond7-p1 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | 35.6ms | 9ms | 4.0x | 9 |

## Coverage by family

| family | CORRECT | DECLINE | WRONG | CRASH | TIMEOUT |
|---|---:|---:|---:|---:|---:|
| adversarial-precision | 4 | 1 | 0 | 0 | 0 |
| binomial-cubic | 32 | 0 | 0 | 0 | 0 |
| binomial-cubic-large | 2 | 0 | 0 | 0 | 0 |
| binomial-quartic | 13 | 4 | 0 | 0 | 0 |
| general-cubic | 4 | 1 | 0 | 0 | 0 |
| general-quartic | 4 | 2 | 0 | 0 | 0 |
| high-m | 8 | 1 | 0 | 0 | 0 |
| large-solutions | 1 | 0 | 0 | 0 | 0 |
| many-solutions | 1 | 0 | 0 | 0 | 0 |
| no-solution | 1 | 1 | 0 | 0 | 0 |
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
| binom-cubic-d2-p1 | `-2*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.8ms | x^3 - 2 y^3 = 1 |
| binom-cubic-d2-m1 | `-2*y^3 + x^3 == -1` | CORRECT | 2/2 | 5.6ms | x^3 - 2 y^3 = -1 |
| binom-cubic-d3-p1 | `-3*y^3 + x^3 == 1` | CORRECT | 1/1 | 2.2ms | x^3 - 3 y^3 = 1 |
| binom-cubic-d3-m1 | `-3*y^3 + x^3 == -1` | CORRECT | 1/1 | 2.2ms | x^3 - 3 y^3 = -1 |
| binom-cubic-d5-p1 | `-5*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.5ms | x^3 - 5 y^3 = 1 |
| binom-cubic-d5-m1 | `-5*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.0ms | x^3 - 5 y^3 = -1 |
| binom-cubic-d6-p1 | `-6*y^3 + x^3 == 1` | CORRECT | 1/1 | 7.3ms | x^3 - 6 y^3 = 1 |
| binom-cubic-d6-m1 | `-6*y^3 + x^3 == -1` | CORRECT | 1/1 | 7.2ms | x^3 - 6 y^3 = -1 |
| binom-cubic-d7-p1 | `-7*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.5ms | x^3 - 7 y^3 = 1 |
| binom-cubic-d7-m1 | `-7*y^3 + x^3 == -1` | CORRECT | 2/2 | 5.5ms | x^3 - 7 y^3 = -1 |
| binom-cubic-d10-p1 | `-10*y^3 + x^3 == 1` | CORRECT | 1/1 | 9.6ms | x^3 - 10 y^3 = 1 |
| binom-cubic-d10-m1 | `-10*y^3 + x^3 == -1` | CORRECT | 1/1 | 9.7ms | x^3 - 10 y^3 = -1 |
| binom-cubic-d11-p1 | `-11*y^3 + x^3 == 1` | CORRECT | 1/1 | 5.9ms | x^3 - 11 y^3 = 1 |
| binom-cubic-d11-m1 | `-11*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.1ms | x^3 - 11 y^3 = -1 |
| binom-cubic-d12-p1 | `-12*y^3 + x^3 == 1` | CORRECT | 1/1 | 5.8ms | x^3 - 12 y^3 = 1 |
| binom-cubic-d12-m1 | `-12*y^3 + x^3 == -1` | CORRECT | 1/1 | 5.9ms | x^3 - 12 y^3 = -1 |
| binom-cubic-d13-p1 | `-13*y^3 + x^3 == 1` | CORRECT | 1/1 | 6.0ms | x^3 - 13 y^3 = 1 |
| binom-cubic-d13-m1 | `-13*y^3 + x^3 == -1` | CORRECT | 1/1 | 6.2ms | x^3 - 13 y^3 = -1 |
| binom-cubic-d15-p1 | `-15*y^3 + x^3 == 1` | CORRECT | 1/1 | 27.5ms | x^3 - 15 y^3 = 1 |
| binom-cubic-d15-m1 | `-15*y^3 + x^3 == -1` | CORRECT | 1/1 | 27.2ms | x^3 - 15 y^3 = -1 |
| binom-cubic-d20-p1 | `-20*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.5ms | x^3 - 20 y^3 = 1 |
| binom-cubic-d20-m1 | `-20*y^3 + x^3 == -1` | CORRECT | 2/2 | 5.7ms | x^3 - 20 y^3 = -1 |
| binom-cubic-d30-p1 | `-30*y^3 + x^3 == 1` | CORRECT | 1/1 | 12.2ms | x^3 - 30 y^3 = 1 |
| binom-cubic-d30-m1 | `-30*y^3 + x^3 == -1` | CORRECT | 1/1 | 11.1ms | x^3 - 30 y^3 = -1 |
| binom-cubic-d41-p1 | `-41*y^3 + x^3 == 1` | CORRECT | 1/1 | 137.9ms | x^3 - 41 y^3 = 1 |
| binom-cubic-d41-m1 | `-41*y^3 + x^3 == -1` | CORRECT | 1/1 | 138.8ms | x^3 - 41 y^3 = -1 |
| binom-cubic-d42-p1 | `-42*y^3 + x^3 == 1` | CORRECT | 1/1 | 30.9ms | x^3 - 42 y^3 = 1 |
| binom-cubic-d42-m1 | `-42*y^3 + x^3 == -1` | CORRECT | 1/1 | 30.0ms | x^3 - 42 y^3 = -1 |
| binom-cubic-d43-p1 | `-43*y^3 + x^3 == 1` | CORRECT | 2/2 | 8.1ms | x^3 - 43 y^3 = 1 |
| binom-cubic-d43-m1 | `-43*y^3 + x^3 == -1` | CORRECT | 2/2 | 9.0ms | x^3 - 43 y^3 = -1 |
| binom-cubic-d17-p1 | `-17*y^3 + x^3 == 1` | CORRECT | 2/2 | 28.5ms | x^3-17y^3=1 -> (18,7); Q(cbrt17) non-monogenic (index 3) |
| binom-cubic-d20-p1 | `-20*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.7ms | x^3-20y^3=1 |
| binom-cubic-plus-d2 | `2*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.6ms | x^3+2y^3=1 |
| binom-cubic-plus-d7 | `7*y^3 + x^3 == 1` | CORRECT | 2/2 | 5.6ms | x^3+7y^3=1 |
| binom-cubic-d2-m2 | `-2*y^3 + x^3 == 2` | CORRECT | 1/1 | 9.5ms | x^3-2y^3=2 |
| binom-cubic-d2-m3 | `-2*y^3 + x^3 == 3` | CORRECT | 2/2 | 11.1ms | x^3-2y^3=3 |
| binom-cubic-d2-m4 | `-2*y^3 + x^3 == 4` | CORRECT | 0/0 | 11.4ms | x^3-2y^3=4 |
| binom-cubic-d2-m5 | `-2*y^3 + x^3 == 5` | CORRECT | 0/0 | 10.8ms | x^3-2y^3=5 |
| binom-cubic-d2-m9 | `-2*y^3 + x^3 == 9` | CORRECT | 0/0 | 14.2ms | x^3-2y^3=9 |
| binom-cubic-d2-m10 | `-2*y^3 + x^3 == 10` | CORRECT | 2/2 | 12.9ms | x^3-2y^3=10 |
| binom-cubic-d2-m100 | `-2*y^3 + x^3 == 100` | CORRECT | 0/0 | 234.5ms | x^3-2y^3=100 |
| binom-cubic-d2-m73 | `-2*y^3 + x^3 == 73` | CORRECT | 0/0 | 10.7ms | x^3-2y^3=73 (has solutions) |
| cyclic-cubic-m2 | `y^3 - 3*x*y^2 + x^3 == 2` | DECLINE | None/0 | 1.8ms | x^3-3xy^2+y^3=2 |
| cyclic-cond9-p1 | `y^3 - 3*x*y^2 + x^3 == 1` | CORRECT | 6/6 | 7.7ms | cond 9, x^3-3xy^2+y^3=1 |
| cyclic-cond9-m1 | `y^3 - 3*x*y^2 + x^3 == -1` | CORRECT | 6/6 | 7.7ms | cond 9, =-1 |
| thomas-cond7-p1 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | CORRECT | 9/9 | 35.6ms | x^3+x^2y-2xy^2-y^3=1 (cond 7, 9 solutions) |
| thomas-cond7-m1 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == -1` | CORRECT | 9/9 | 34.7ms | cond 7, =-1 |
| thomas-t1 | `-y^3 - 3*x*y^2 - x^2*y + x^3 == 1` | CORRECT | 3/3 | 0.7ms | x^3-1x^2y-3xy^2-y^3=1 |
| thomas-t2 | `-y^3 - 4*x*y^2 - 2*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.7ms | x^3-2x^2y-4xy^2-y^3=1 |
| thomas-t3 | `-y^3 - 5*x*y^2 - 3*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.7ms | x^3-3x^2y-5xy^2-y^3=1 |
| thomas-t4 | `-y^3 - 6*x*y^2 - 4*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.9ms | x^3-4x^2y-6xy^2-y^3=1 |
| thomas-t5 | `-y^3 - 7*x*y^2 - 5*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.8ms | x^3-5x^2y-7xy^2-y^3=1 |
| thomas-t10 | `-y^3 - 12*x*y^2 - 10*x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.7ms | x^3-10x^2y-12xy^2-y^3=1 |
| shanks-simplest-n3 | `y^3 - 3*x*y^2 - 6*x^2*y + x^3 == 1` | CORRECT | 4/4 | 11.6ms | Shanks x^3-3x^2y-6xy^2-... variant |
| cubic-irr-1 | `y^3 + x*y^2 - 2*x^2*y + x^3 == 1` | CORRECT | 4/4 | 2.7ms | x^3-2x^2y+xy^2+y^3 |
| cubic-irr-2 | `-y^3 + 3*x*y^2 - 3*x^2*y + x^3 == 1` | DECLINE | None/None | 1.1ms | totally real-ish |
| cubic-irr-3 | `y^3 - x^2*y + x^3 == 1` | CORRECT | 5/5 | 2.6ms | x^3-xy^2+y^3? small disc |
| cubic-neg-disc | `y^3 + x*y^2 + x^2*y + x^3 == 1` | CORRECT | 2/2 | 0.5ms | x^3+x^2y+xy^2+y^3 (reducible? check) |
| cubic-big-coef | `-5*y^3 + 11*x*y^2 - 7*x^2*y + x^3 == 1` | CORRECT | 1/1 | 0.5ms | bigger coeffs |
| nosol-d2-m4 | `-2*y^3 + x^3 == 4` | CORRECT | 0/0 | 12.0ms | x^3-2y^3=4 (|m|!=1 -> decline; PARI {}) |
| nosol-cyclic | `y^3 - 3*x*y^2 + x^3 == 4` | DECLINE | None/0 | 1.9ms | cyclic, =4 |
| binom-quartic-d2-p1 | `-2*y^4 + x^4 == 1` | CORRECT | 2/2 | 5.8ms | x^4 - 2 y^4 = 1 |
| binom-quartic-d2-m1 | `-2*y^4 + x^4 == -1` | CORRECT | 4/4 | 6.3ms | x^4 - 2 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d3-p1 | `-3*y^4 + x^4 == 1` | CORRECT | 2/2 | 9.3ms | x^4 - 3 y^4 = 1 |
| binom-quartic-d3-m1 | `-3*y^4 + x^4 == -1` | CORRECT | 0/0 | 9.2ms | x^4 - 3 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d5-p1 | `-5*y^4 + x^4 == 1` | CORRECT | 6/6 | 15.5ms | x^4 - 5 y^4 = 1 |
| binom-quartic-d5-m1 | `-5*y^4 + x^4 == -1` | CORRECT | 0/0 | 15.9ms | x^4 - 5 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d6-p1 | `-6*y^4 + x^4 == 1` | CORRECT | 2/2 | 39.8ms | x^4 - 6 y^4 = 1 |
| binom-quartic-d6-m1 | `-6*y^4 + x^4 == -1` | CORRECT | 0/0 | 40.1ms | x^4 - 6 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d7-p1 | `-7*y^4 + x^4 == 1` | CORRECT | 2/2 | 17.8ms | x^4 - 7 y^4 = 1 |
| binom-quartic-d7-m1 | `-7*y^4 + x^4 == -1` | CORRECT | 0/0 | 17.8ms | x^4 - 7 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d8-p1 | `-8*y^4 + x^4 == 1` | CORRECT | 2/2 | 5.9ms | x^4 - 8 y^4 = 1 |
| binom-quartic-d8-m1 | `-8*y^4 + x^4 == -1` | CORRECT | 0/0 | 6.1ms | x^4 - 8 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d10-p1 | `-10*y^4 + x^4 == 1` | DECLINE | None/2 | 65.9ms | x^4 - 10 y^4 = 1 |
| binom-quartic-d10-m1 | `-10*y^4 + x^4 == -1` | DECLINE | None/0 | 66.9ms | x^4 - 10 y^4 = -1 (subfield unit -> case iii) |
| binom-quartic-d17-p1 | `-17*y^4 + x^4 == 1` | DECLINE | None/2 | 73.3ms | x^4 - 17 y^4 = 1 |
| binom-quartic-d17-m1 | `-17*y^4 + x^4 == -1` | DECLINE | None/4 | 69.9ms | x^4 - 17 y^4 = -1 (subfield unit -> case iii) |
| ljunggren-x4-2y4 | `-2*y^4 + x^4 == -1` | CORRECT | 4/4 | 5.8ms | Ljunggren-type x^4-2y^4=-1 |
| quartic-biquad-1 | `y^4 - 4*x^2*y^2 + x^4 == 1` | CORRECT | 12/12 | 44.7ms | x^4-4x^2y^2+y^4 (biquadratic, quadratic subfield) |
| quartic-tr-1 | `y^4 - 4*x*y^3 + 6*x^2*y^2 - 4*x^3*y + x^4 == 1` | DECLINE | None/None | 1.3ms | (x-y)^4+... totally real? |
| quartic-mixed-1 | `-y^4 + 2*x^3*y + x^4 == 1` | CORRECT | 2/2 | 33.6ms | x^4+2x^3y-y^4 |
| quartic-cyclotomic | `y^4 + x*y^3 + x^2*y^2 + x^3*y + x^4 == 1` | DECLINE | None/6 | 2.1ms | x^4+x^3y+x^2y^2+xy^3+y^4 (cyclotomic, reducible over cyclo) |
| quartic-x4-x-1 | `-y^4 - x*y^3 + x^4 == 1` | CORRECT | 4/4 | 87.0ms | x^4-xy^3-y^4 |
| quartic-totallyimag | `y^4 + x^2*y^2 + x^4 == 1` | CORRECT | 4/4 | 0.5ms | x^4+x^2y^2+y^4 (may be reducible) |
| binom-quintic-d2 | `-2*y^5 + x^5 == 1` | CORRECT | 2/2 | 17.5ms | x^5-2y^5=1 |
| binom-quintic-d3 | `-3*y^5 + x^5 == 1` | CORRECT | 1/1 | 50.1ms | x^5-3y^5=1 |
| binom-quintic-d5 | `-5*y^5 + x^5 == 1` | DECLINE | None/1 | 155.1ms | x^5-5y^5=1 |
| quintic-general | `y^5 - x*y^4 - x^4*y + x^5 == 1` | CORRECT | 2/2 | 0.6ms | general quintic form |
| binom-sextic-d2 | `-2*y^6 + x^6 == 1` | CORRECT | 2/2 | 128.4ms | x^6-2y^6=1 |
| binom-sextic-d3 | `-3*y^6 + x^6 == 1` | CORRECT | 2/2 | 242.2ms | x^6-3y^6=1 |
| septic-x7-2 | `-2*y^7 + x^7 == 1` | CORRECT | 2/2 | 176.6ms | x^7-2y^7=1 |
| nonmono-d17 | `-17*y^3 + x^3 == 1` | CORRECT | 2/2 | 26.8ms | Q(cbrt17) index 3 at p=3 |
| nonmono-d19 | `-19*y^3 + x^3 == 1` | CORRECT | 2/2 | 6.1ms | x^3-19y^3=1 |
| nonmono-cubic-idx | `-8*y^3 - 2*x*y^2 - x^2*y + x^3 == 1` | CORRECT | 1/1 | 28.0ms | Dedekind cubic t^3-t^2-2t-8 (index 2 at 2) |
| reducible-x3-y3 | `-y^3 + x^3 == 7` | CORRECT | 2/2 | 0.4ms | x^3 - y^3 = 7 -> (x-y)(x^2+xy+y^2); reducible |
| reducible-x4-y4 | `-y^4 + x^4 == 15` | CORRECT | 4/4 | 0.7ms | x^4 - y^4 = 15 reducible |
| reducible-biquad | `4*y^4 - 5*x^2*y^2 + x^4 == 3` | CORRECT | 0/0 | 0.5ms | x^4-5x^2y^2+4y^4=(x^2-y^2)(x^2-4y^2) |
| adv-close-roots-1 | `y^4 - 7*x*y^3 + 14*x^2*y^2 - 7*x^3*y + x^4 == 1` | DECLINE | None/4 | 1.6ms | near-symmetric quartic, clustered roots |
| adv-big-coef-cubic | `-97*y^3 + x^3 == 1` | CORRECT | 1/1 | 121.5ms | x^3-97y^3=1 |
| adv-big-coef-2 | `131*y^3 - 211*x*y^2 + 97*x^2*y + x^3 == 1` | CORRECT | 1/1 | 379.2ms | large mixed coeffs |
| adv-scaled-1 | `-999*y^3 + x^3 == 1` | CORRECT | 2/2 | 23.6ms | x^3-999y^3=1 (999=3^3*37) |
| adv-large-reg | `-y^3 - x*y^2 - x^2*y + x^3 == 1` | CORRECT | 4/4 | 6.1ms | x^3-x^2y-xy^2-y^3=1 (tribonacci field, larger regulator) |
| many-sol-cond7 | `-y^3 - 2*x*y^2 + x^2*y + x^3 == 1` | CORRECT | 9/9 | 37.0ms | Thomas cond7 (9 solutions) |
| large-sol-1 | `-3*y^3 + x^3 == 1` | CORRECT | 1/1 | 2.2ms | x^3-3y^3=1 |
