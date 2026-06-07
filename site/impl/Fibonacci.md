---
source: src/fibonacci.c
---
**Algorithm.** `builtin_fibonacci` follows the system's two-tier split. For exact integer order `n`, `Fibonacci[n]` uses GMP **fast doubling** (`fib_mpz`: `F(2k)=F(k)(2F(k+1)-F(k))`, `F(2k+1)=F(k+1)^2+F(k)^2`) in `O(log n)` big-integer operations, with negative orders via `F(-m)=(-1)^{m+1}F(m)`. `Fibonacci[n, x]` (Fibonacci polynomial) iterates the recurrence `F_k = x F_{k-1} + F_{k-2}`, evaluating each step so the partial result stays canonical. For inexact or non-integer order, it builds the generalized closed form `(phi^n - Cos[Pi n] phi^-n)/Sqrt[5]` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2` for the polynomial) as an expression tree and lets `numericalize` drive the reduction at the inputs' precision — which also yields complex results for complex order for free. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` for the integer fast-doubling pair; `Expr` trees (built with `mk_fn1`/`mk_fn2` helpers, reduced via `eval_and_free`) for polynomial and closed-form paths.
