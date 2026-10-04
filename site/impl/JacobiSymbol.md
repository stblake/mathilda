---
references:
  - "K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — the Legendre, Jacobi and Kronecker symbols and quadratic reciprocity (Chapter 5)."
source: src/numbertheory/jacobisymbol.c
---
**Algorithm.** `builtin_jacobisymbol` requires two integer-like arguments and returns
`mpz_kronecker(n, m)`, i.e. the full Kronecker-symbol generalisation of the Jacobi symbol:
`m` need not be odd or positive, and `n` may be negative. For prime `m` this is the Legendre
symbol — `+1` if `n` is a non-zero quadratic residue modulo `m`, `-1` for a non-residue, `0`
when `m` divides `n`; for odd composite `m` it is the product of the Legendre symbols over
the prime factors, computed by quadratic reciprocity without factoring `m`. The result is
always `-1`, `0` or `1`.

**Data structures.** Two `mpz_t`, extracted with `expr_to_mpz`; the answer is a machine
`Integer`. There is no ND/packed/`Compile` kernel — `Listable` threading over list arguments
is applied by the evaluator before the builtin is reached.

**Complexity / limits.** `mpz_kronecker` is `O((log m)^2)`, so machine integers and
arbitrary-precision bigints are handled uniformly (e.g. `JacobiSymbol[10^10 + 1,
Prime[1000]]`). Symbolic arguments leave the call unevaluated, and a wrong argument count
emits `JacobiSymbol::argrx`.
