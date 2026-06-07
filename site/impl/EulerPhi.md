---
source: src/facint.c
---
`builtin_eulerphi` computes Euler's totient. It takes `|n|` (since `phi(-n)=phi(n)`), factors a working copy via the shared `factorize_mpz` cascade (trial division → Pollard rho → ECM), then applies `phi(n) = n * prod (1 - 1/p_i)` per distinct prime as `phi <- (phi / p) * (p - 1)` with GMP `mpz_divexact`/`mpz_mul`, keeping intermediates exact. `phi(0) = 0`, `phi(1) = 1`. Non-integer arguments return `NULL`. Its cost is dominated by the factorisation of `n`.
