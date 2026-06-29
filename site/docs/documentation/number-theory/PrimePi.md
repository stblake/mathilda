# PrimePi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimePi[x] gives the number of primes less than or equal to x, exact for x up to 5*10^13 (larger x is left unevaluated). The option Method -> m selects the algorithm: Automatic (default), "Sieve", "Legendre", "Meissel", "Lehmer", "LMO" (Lagarias-Miller-Odlyzko), "DelegliseRivat", or "LucyHedgehog".
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimePi[10^9]
Out[1]= 50847534

In[2]:= PrimePi[10^9, Method -> "LMO"]
Out[2]= 50847534

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## Implementation notes

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

- `Listable`, `Protected`.  Default `Method -> Automatic`.
- Several recognised prime-counting algorithms are available via `Method`:

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- M. Deléglise and J. Rivat, "Computing π(x): the Meissel, Lehmer, Lagarias, Miller, Odlyzko method", Math. Comp. 65 (1996), 235–245.
- J. C. Lagarias, V. S. Miller and A. M. Odlyzko, "Computing π(x): the Meissel–Lehmer method", Math. Comp. 44 (1985), 537–560.
- D. H. Lehmer, "On the exact number of primes less than a given limit", Illinois J. Math. 3 (1959), 381–388.
- Source: [`src/numbertheory/primecount.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/primecount.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PrimePi[100]
Out[1]= 25
```

The count scales to large bounds; there are 78498 primes below one million:

```mathematica
In[1]:= PrimePi[10^6]
Out[1]= 78498
```

`PrimePi` threads over lists (it is `Listable`):

```mathematica
In[1]:= PrimePi[{10, 100}]
Out[1]= {4, 25}
```

### Choosing a method

`PrimePi[x, Method -> m]` counts with a specific algorithm. Every method returns
the same `π(x)` (they cross-validate one another in the test suite); they differ
only in speed and in the largest `x` they support. The default
`Method -> Automatic` picks the best one for the size of `x`.

```mathematica
In[1]:= PrimePi[10^9, Method -> "DelegliseRivat"]
Out[1]= 50847534

In[2]:= PrimePi[10^9, Method -> "LMO"]
Out[2]= 50847534
```

| `Method` | What it does | Best for | Max `x` |
|----------|--------------|----------|---------|
| `Automatic` | table ≤ `10^6`, Lucy_Hedgehog ≤ `10^9`, Deléglise–Rivat above | everything | `5×10^13` |
| `"Sieve"` | segmented sieve of Eratosthenes; counts every prime directly. The slow-but-obvious ground truth used to check the others. | small `x`, validation | `10^10` |
| `"LucyHedgehog"` | Lucy_Hedgehog dynamic programming over the `O(√x)` distinct values of `⌊x/i⌋`. `O(x^{3/4})` time, `O(√x)` memory, almost no setup. | `x ≲ 10^9` | `5×10^13` |
| `"Legendre"` | Legendre's formula `π(x) = φ(x, π(√x)) + π(√x) − 1`. Simplest combinatorial method, but the partial sieve `φ` needs a table up to `x`, so it is the slowest and most memory-hungry here. | teaching / small `x` | `10^9` |
| `"Meissel"` | Meissel's `π(x) = φ(x, a) + a − 1 − P2` with `a = π(x^{1/3})`; the `P2` term corrects for products of two large primes. Sieve table only up to `x^{2/3}`. | mid-range `x` | `10^13` |
| `"Lehmer"` | Lehmer's refinement of Meissel using `a = π(x^{1/4})` and an extra `P3` correction. | mid-range `x` | `10^12` |
| `"LMO"` | Lagarias–Miller–Odlyzko. Splits `φ(x, a)` into *ordinary* and *special* leaves and sweeps the special leaves through a segmented sieve with a Fenwick tree — `O(√x)` memory, far fewer operations than a sieve to `x`. | large `x` | `5×10^13` |
| `"DelegliseRivat"` | Refines LMO by sorting the special leaves into *trivial*, *easy* and *hard*, answering the first two with `O(1)` prime-count lookups and sending only the *hard* leaves through the sieve. The fastest combinatorial method here. | large `x` | `5×10^13` |

**Relative performance.** Below `10^6` every method just reads the cached sieve
table (instant). For larger `x` the practical ranking is, fastest first,
Deléglise–Rivat ≳ LMO > Lehmer > Meissel > Legendre, with Lucy_Hedgehog
competitive up to about `10^9` thanks to its negligible setup cost (which is why
`Automatic` prefers it there). On a single core, Deléglise–Rivat counts to
`10^12` in ≈ 1.6 s and `10^13` in ≈ 9.9 s, versus ≈ 2.5 s and 13.5 s for
Lucy_Hedgehog; both finish in a few hundredths of a second at `10^9` and below.
The combinatorial methods are exact up to `5×10^13`; beyond that, and for the
range-capped methods past their limit, `PrimePi` is left unevaluated.

### Notes

`PrimePi[x]` gives the prime-counting function `π(x)`, the number of primes
less than or equal to `x`. For `x = 100` the answer is `25`; for `x = 10^6` it
is `78498`, consistent with the asymptotic `π(x) ~ x / Log[x]`. An unrecognised
`Method` setting emits `PrimePi::method` and leaves the call unevaluated.
