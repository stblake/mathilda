---
source: src/lucas.c
---
**Algorithm.** `builtin_lucasl` mirrors `Fibonacci`. For exact integer order it fast-doubles the Fibonacci pair `(F_m, F_{m+1})` in GMP and derives `L_m = 2 F_{m+1} - F_m` in `O(log n)`, with negative orders via `L_{-m} = (-1)^m L_m`. `LucasL[n, x]` (Lucas polynomial) iterates `L_k = x L_{k-1} + L_{k-2}` from `L_0 = 2, L_1 = x`, Expand-ing each step. For inexact/non-integer order it builds the closed form `phi^n + Cos[Pi n] phi^-n` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2`) and hands it to `numericalize`. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` integer pair; `Expr` trees through `eval_and_free` for the polynomial and closed-form branches.
