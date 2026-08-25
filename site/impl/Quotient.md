---
source: src/core.c
---
`builtin_quotient` floors the ratio for real arguments (`mpz_fdiv_q` on exact integers), so that `n Quotient[m, n] + Mod[m, n] == m` holds exactly; the three-argument form applies the same `d`-offset as `Mod`. For complex `m` or `n` it switches to Gaussian-integer division — it forms the exact ratio and rounds each component to the nearest integer (ties to even), which is the quotient minimising the norm of the remainder and deliberately differs from the real (floored) branch. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; non-numeric arguments stay symbolic.
