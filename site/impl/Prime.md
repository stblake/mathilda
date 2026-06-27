---
references:
  - "M. Cipolla, \"La determinazione assintotica dell'n-esimo numero primo\", Rend. Accad. Sci. Fis. Mat. Napoli 8 (1902), 132–166."
source: src/numbertheory/prime.c
---
**Algorithm.** `builtin_prime` (`src/numbertheory/prime.c`) computes `p_n`, the
`n`-th prime, as the functional inverse of `PrimePi`. The index `n` is resolved
to a positive `int64_t` (Integer or Bigint); a non-positive integer emits
`Prime::intpp`, any other numeric argument also emits `Prime::intpp`, and a
symbolic argument is left unevaluated. A wrong argument count emits `Prime::argx`.

**Small `n` — direct table.** `primecount_init` builds (once) a sieve table of
the primes below `10^6`; if `n` is within that table, `Prime` returns
`primecount_small_prime(n)` directly — a single indexed read.

**Large `n` — estimate, refine, walk.** Beyond the table:

1. **Estimate.** `prime_estimate` evaluates Cipolla's asymptotic expansion
   `p_n ≈ n(ln n + ln ln n − 1 + (ln ln n − 2)/ln n − (ln²ln n − 6 ln ln n +
   11)/(2 ln²n))`, accurate to a relative error that shrinks with `n`.
2. **Newton-refine.** Starting from that anchor `x` with `c = π(x)` (computed by
   `prime_count(x, PC_AUTOMATIC)`), it takes up to `PRIME_NEWTON_MAX = 12` Newton
   steps `x ← x + (n − c)·ln x` against the *exact* counter until `|n − c|` falls
   below `PRIME_WALK_LIMIT = 10^6`.
3. **Walk.** From the refined anchor, GMP's `mpz_nextprime` / `mpz_prevprime`
   step exactly `|n − c|` primes in the right direction to land on `p_n`. The
   result is returned as a normalized bigint.

**Range.** The Newton search relies on `PrimePi`, so `Prime` is exact while the
estimate stays within `PI_COUNT_MAX = 5×10^13` — i.e. `n` up to about
`1.4×10^12` (`Prime[10^10] = 252097800623`). For `n` whose estimated prime would
exceed that bound, the call is left unevaluated.

**Attributes.** `Listable` (threads over a list of indices), `Protected`.
