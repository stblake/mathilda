---
source: src/numbertheory.c
---
**Algorithm.** `builtin_multiplicativeorder` computes the multiplicative order of `k` modulo `n` — the least `m` with `k^m ≡ 1 (mod n)`. `mo_order_mpz` reduces `k mod n`, checks `gcd(k, n) = 1` (returning unevaluated otherwise), computes Euler's totient `φ(n)` (`mo_eulerphi_mpz`), and then deflates the order down from `φ(n)`: for each distinct prime `q | φ(n)` it repeatedly divides the running order by `q` while `k^(order/q) ≡ 1 (mod n)` (using `mpz_powm`). This yields the order without enumerating all exponents. The 3-argument form `MultiplicativeOrder[k, n, {r1, ...}]` instead searches for the least `m ≤ order` with `k^m` congruent to one of the residues `r_i` (`mo_search_residues`, capped at `MO_SEARCH_CAP = 10^8` iterations and requiring the order to fit in `unsigned long`). Wrong arg counts emit `MultiplicativeOrder::argt`.

**Data structures.** GMP `mpz_t` throughout; the distinct prime divisors of `φ(n)` are collected into a fixed `mpz_t primes[]` array (`pr_collect_distinct_primes`).

**Complexity / limits.** Order computation is `O(ω(φ(n)) · log φ(n))` modular exponentiations after factoring `φ(n)`; the residue-search form is bounded by `MO_SEARCH_CAP`.
