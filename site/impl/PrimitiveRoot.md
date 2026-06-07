---
source: src/numbertheory.c
---
**Algorithm.** `builtin_primitiveroot` returns the smallest primitive root of `n` ≥ an optional second-argument start (`PrimitiveRoot[n]` / `PrimitiveRoot[n, k]`). It first classifies `n` with `pr_classify` to confirm the unit group `(Z/nZ)*` is cyclic (i.e. `n ∈ {1, 2, 4, p^e, 2p^e}` for odd prime p), then computes `φ(n)` and its distinct prime divisors. `pr_smallest_primitive_root` scans candidates `g`, testing each with `pr_is_primitive_root`: `g` is a primitive root iff `gcd(g, n) = 1` and `g^(φ(n)/q) ≢ 1 (mod n)` for every prime `q | φ(n)` (via `mpz_powm`). Non-integer numeric input emits `PrimitiveRoot::intg`; `n < 2` likewise; wrong arg count emits `PrimitiveRoot::argt`; symbolic input returns unevaluated.

**Data structures.** GMP `mpz_t`; distinct primes of `φ(n)` in a fixed `mpz_t[]` array.

**Complexity / limits.** Primitive-root density is `φ(φ(n))/φ(n)`, so the scan finds one in roughly `O(log log p)` candidates on average; each test is `ω(φ(n))` modular exponentiations.
