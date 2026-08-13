# Prime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Prime[n] gives the nth prime number. Listable. Small n is read from a sieve table; large n inverts PrimePi via an asymptotic estimate refined against the exact prime counter. Defined for positive integers up to n ~ 1.4*10^12; non-positive-integer arguments give Prime::intpp.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Prime[100]
Out[1]= 541

In[2]:= Prime[{1, 3, 4, 10}]
Out[2]= {2, 5, 7, 29}

In[3]:= Prime[10^10]
Out[3]= 252097800623
```

### Applications (3)

```mathematica
In[4]:= Prime[100]
Out[4]= 541

In[5]:= Prime[{1, 3, 4, 10}]
Out[5]= {2, 5, 7, 29}

In[6]:= Prime[10^10]
Out[6]= 252097800623
```

## Options & behaviour

> **Packed arrays.** `Prime[list]` over an `int64` buffer is one **sieve** to
> the largest requested index followed by a gather, rather than a separate
> prime count per element — so the whole array costs about what its largest
> single index used to.

## Algorithm

prime.c -- Prime[n] (the nth prime) and PrimePi[x] (the prime-counting function).

The prime-counting algorithms live in primecount.c; this file holds the two

```text
builtins.  PrimePi[x] accepts a Method option selecting the algorithm:
  Automatic (default), "Sieve", "Legendre", "Meissel", "Lehmer", "LMO",
  "DelegliseRivat", "LucyHedgehog".
```

Prime[n] is the functional inverse of PrimePi: small n are read straight from the sieve table; large n are found by seeding Cipolla's asymptotic estimate, refining it with a Newton step driven by the exact counter, then walking with GMP's nextprime/prevprime to land exactly on p_n.

## Implementation notes

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

- `Listable`, `Protected`.
- Small `n` is read directly from a sieve table of the primes below $10^6$.
- Large `n` inverts `PrimePi`: a Cipolla asymptotic estimate for $p_n$ is refined
  by a Newton step against the exact prime counter, then `NextPrime`/`PrevPrime`
  steps land exactly on $p_n$.  Exact for `n` up to about $1.4 \times 10^{12}$
  ($p_n \le 5 \times 10^{13}$); beyond that the call is left unevaluated.
- A non-positive-integer argument emits `Prime::intpp`; a wrong argument count
  emits `Prime::argx`; both leave the call unevaluated.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [PrimePi](../../number-theory/PrimePi/), [NextPrime](../../number-theory/NextPrime/)

- M. Cipolla, "La determinazione assintotica dell'n-esimo numero primo", Rend. Accad. Sci. Fis. Mat. Napoli 8 (1902), 132–166.
- Source: [`src/numbertheory/prime.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/prime.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_jacobisymbol.c`](https://github.com/stblake/mathilda/blob/main/tests/test_jacobisymbol.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_prime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_prime.c)
- Tests: [`tests/test_primitive_root.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primitive_root.c)

## Notes & additional examples

### Notes

`Prime[n]` gives the `n`-th prime `p_n`, so `Prime[1] = 2`, `Prime[2] = 3`,
`Prime[3] = 5`, and so on. It is the functional inverse of [`PrimePi`](PrimePi.md):
`PrimePi[Prime[n]] == n` for every positive integer `n`.

Small `n` is read straight from the sieve table of the primes below `10^6`. For
larger `n`, `Prime` seeds Cipolla's asymptotic estimate of `p_n`, refines it with
a Newton step driven by the exact prime counter, then walks with
`NextPrime`/`PrevPrime` to land exactly on `p_n`. This is exact for `n` up to
about `1.4×10^12` (`p_n ≤ 5×10^13`); beyond that the call is left unevaluated.

A non-positive-integer argument emits `Prime::intpp` and a wrong argument count
emits `Prime::argx`; in both cases the call is returned unevaluated.
