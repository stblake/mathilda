---
references:
  - "M. Deléglise and J. Rivat, \"Computing π(x): the Meissel, Lehmer, Lagarias, Miller, Odlyzko method\", Math. Comp. 65 (1996), 235–245."
  - "J. C. Lagarias, V. S. Miller and A. M. Odlyzko, \"Computing π(x): the Meissel–Lehmer method\", Math. Comp. 44 (1985), 537–560."
  - "D. H. Lehmer, \"On the exact number of primes less than a given limit\", Illinois J. Math. 3 (1959), 381–388."
source: src/numbertheory/primecount.c
---
**Algorithm.** `builtin_primepi` (`src/numbertheory/prime.c`) floors the argument
to an `int64_t` (Integer, Bigint, Real and Rational accepted; `x < 2` returns
`0`), parses any trailing `Method -> _` options, and delegates the count to
`prime_count(x, method)` in `src/numbertheory/primecount.c`. A count of `-1`
(argument beyond the chosen method's range) leaves the call unevaluated.

**Shared infrastructure (`primecount.c`).** All combinatorial methods build on a
common toolkit:

- a once-built sieve table of the primes below `10^6` (`primecount_init`, also
  serving `Prime[n]`), with `pi_small` answering `π(x)` for `x ≤ 10^6` by binary
  search;
- `gen_primes(bound)` — a plain sieve returning every prime up to a modest bound
  (`≤ √x`);
- `PiTable` — a bit-packed sieve up to `B` with per-word prefix popcounts, giving
  `O(1)` `π(v)` lookups for `v ≤ B`;
- `phi(v, b)` — the partial-sieve function `#{n ≤ v : n coprime to the first b
  primes}`, evaluated via a *PhiTiny* wheel over the first 7 primes (`phi_wheel`)
  plus **Lehmer's `p_b² > v` prune** (`phi_rec`): once the `b`-th prime squared
  exceeds `v`, an `n ≤ v` coprime to the first `b` primes is either `1` or a
  single prime in `(p_b, v]`, so `φ(v, b) = π(v) − b + 1` (answered from the
  `PiTable`). This prune is what turns the otherwise-exponential Legendre
  recursion into a bounded one.

**Methods (`prime_count` dispatch).** Each is a recognised setting for `Method`;
they overlap in range so the test suite cross-validates them against one another
and against the known `π(10^k)`:

| Method | Identity / technique | Range cap | Notes |
|--------|----------------------|-----------|-------|
| `"Sieve"` | segmented sieve of Eratosthenes — counts every prime directly | `10^10` | ground-truth oracle; `O(x log log x)` time, `O(√x)` space |
| `"LucyHedgehog"` | Lucy_Hedgehog dynamic program over the `O(√x)` distinct values of `⌊x/i⌋` | `5×10^13` | `O(x^{3/4})` time, `O(√x)` space; tiny setup, wins for small `x` |
| `"Legendre"` | `π(x) = φ(x, π(√x)) + π(√x) − 1` | `10^9` | `PiTable` spans `[0, x]`; pedagogical, slowest combinatorial |
| `"Meissel"` | `π(x) = φ(x, a) + a − 1 − P2`, `a = π(x^{1/3})` | `10^13` | `PiTable` spans `[0, x^{2/3}]`; `P2` corrects two-prime products |
| `"Lehmer"` | Meissel refined to `a = π(x^{1/4})` with an added `P3` term | `10^12` | `PiTable` spans `[0, x^{3/4}]` |
| `"LMO"` | Lagarias–Miller–Odlyzko: `φ(x, a)` split into ordinary leaves `S1` and special leaves `S2`, the latter swept over `[1, x/y]` in segments carrying a Fenwick tree | `5×10^13` | `O(x^{2/3} / log x)`-ish time, `O(√x)` space |
| `"DelegliseRivat"` | LMO with the special leaves partitioned by quotient size (see below) | `5×10^13` | fastest combinatorial method here |

**Deléglise–Rivat refinement.** `count_dr` keeps the Meissel skeleton and LMO's
`S1`/`P2`, but classifies each special leaf `v = ⌊x / (p_b · m)⌋` by size:

- *trivial* (`v < p_b`): `φ(v, b−1) = 1`;
- *easy* (`p_{b−1}² > v`): `φ(v, b−1) = π(v) − (b−1) + 1`, an `O(1)` `PiTable`
  lookup (the Lehmer prune again);
- *hard* (`p_{b−1}² ≤ v`): the only leaves that pay for an incremental Fenwick
  query.

Trivial and easy leaves dominate the leaf count, so only a minority touch the
segmented sieve — strictly cheaper than LMO, which answers every special leaf
with a Fenwick query. Because `S1`/`P2`/`S2` are otherwise identical to LMO, the
two methods agree exactly.

**`Automatic` dispatch.** `x ≤ 10^6` reads the direct sieve table; `x ≤ 10^9`
uses Lucy_Hedgehog (smallest setup); above `10^9`, Deléglise–Rivat. The
combinatorial methods exact-count up to `PI_COUNT_MAX = 5×10^13`; larger `x` is
reported out of range and the call is left unevaluated.

**Relative performance.** Measured wall-clock (single core): Deléglise–Rivat
≈ 0.09 s at `10^10`, 1.6 s at `10^12`, 9.9 s at `10^13`; Lucy_Hedgehog 0.12 s /
2.5 s / 13.5 s at the same points. At or below `10^9` both finish in a few
hundredths of a second, where Lucy's lighter setup edges ahead — hence the
`Automatic` crossover at `10^9`.
