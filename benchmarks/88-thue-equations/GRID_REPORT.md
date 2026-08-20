# Benchmark 88 — randomized Thue stress grid

Deterministic random binary forms F(x,y) (degree 3-6, mixed m) vs **PARI/GP `thue()`**. Seed `20260820`, 400 cases. Reproducible: the same seed regenerates the identical corpus, so any `WRONG` is a stable, re-runnable completeness bug. Regenerate: `python3 grid.py --n 400 --seed 20260820`.

| verdict | n | meaning |
|---|---|---|
| CORRECT | 261 | finite set matches PARI |
| DECLINE | 138 | unevaluated (honest gap); PARI solved |
| UNVERIFIED | 1 | PARI rejected the form (perfect power / repeated factor) |

**Result: 0 bug(s).** Every form Mathilda solved matched PARI exactly.

## Solve paths exercised: 261 CORRECT

A sample (the grid's value is that these are *machine-generated*, not chosen to look good):

| label | form == m | #sol |
|---|---|---:|
| grid-0003-n4 | `y^4 + 4*x*y^3 + 4*x^2*y^2 + 2*x^3*y - x^4 == -1` | 2 |
| grid-0004-n4 | `2*y^4 + x*y^3 - x^4 == 4` | 0 |
| grid-0005-n3 | `-y^3 - 4*x*y^2 - 3*x^2*y + x^3 == 1` | 9 |
| grid-0010-n3 | `3*y^3 + 2*x*y^2 - x^2*y + x^3 == -8` | 2 |
| grid-0011-n5 | `-2*y^5 + 2*x*y^4 + x^3*y^2 - 2*x^4*y + x^5 == -1` | 1 |
| grid-0012-n5 | `2*y^5 - 3*x*y^4 - 2*x^2*y^3 - 2*x^3*y^2 - x^4*y - x^5 == -1` | 1 |
| grid-0014-n3 | `2*y^3 - 3*x*y^2 + x^3 == 1` | 1 |
| grid-0016-n3 | `-5*y^3 - 2*x*y^2 - 5*x^2*y + x^3 == -16` | 0 |
| grid-0017-n4 | `-3*y^4 - x*y^3 - x^2*y^2 - 3*x^3*y - x^4 == -1` | 4 |
| grid-0018-n6 | `-y^6 - x*y^5 + 2*x^3*y^3 + x^5*y + x^6 == 1` | 2 |
| grid-0019-n4 | `-2*y^4 + 3*x*y^3 - 2*x^3*y + x^4 == 3` | 0 |
| grid-0020-n3 | `-4*y^3 + 6*x*y^2 + x^2*y + x^3 == -1` | 1 |
| grid-0021-n3 | `-y^3 - 4*x*y^2 - 4*x^2*y + x^3 == -1` | 3 |
| grid-0022-n4 | `-y^4 - x*y^3 + 3*x^2*y^2 + x^3*y - x^4 == -1` | 4 |
| grid-0023-n3 | `-2*y^3 - x*y^2 + x^2*y - x^3 == -39` | 0 |
| grid-0025-n5 | `3*y^5 + 3*x^2*y^3 + 2*x^3*y^2 + 2*x^4*y - x^5 == -1` | 1 |
| grid-0026-n4 | `2*y^4 - 4*x*y^3 + x^2*y^2 + 4*x^3*y - x^4 == -1` | 2 |
| grid-0030-n4 | `-4*y^4 - 4*x*y^3 + 4*x^2*y^2 - 4*x^3*y + x^4 == 1` | 2 |
| grid-0031-n4 | `-2*y^4 - 2*x*y^3 - x^2*y^2 + 2*x^3*y + x^4 == 1` | 2 |
| grid-0033-n3 | `-2*y^3 + 6*x*y^2 + 5*x^2*y + x^3 == -1` | 1 |

