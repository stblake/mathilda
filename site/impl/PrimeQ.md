---
source: src/facint.c
---
`builtin_primeq` is a `*Q` predicate: it always returns `True` or `False`, never unevaluated. For an `EXPR_INTEGER`/`EXPR_BIGINT` it takes `|n|` and runs GMP's `mpz_probab_prime_p(n, 25)` (Baillie–PSW plus 25 Miller–Rabin rounds). With `GaussianIntegers -> True` (parsed by `primeq_parse_options`; a malformed option list yields `False`), a rational integer is a Gaussian prime iff `|n|` is prime and `|n| ≡ 3 (mod 4)`, and a `Complex[a, b]` with integer parts is tested by `gaussian_prime_test` — pure-real/pure-imaginary need the `≡ 3 mod 4` condition, mixed needs `a^2 + b^2` prime. Reals, rationals, strings, symbols, and symbolic functions are all `False`.
