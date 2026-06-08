# PrimePi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimePi[x] gives the number of primes less than or equal to x.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimePi[10]
Out[1]= 4

In[2]:= PrimePi[100]
Out[2]= 25

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## Implementation notes

**Algorithm.** `builtin_primepi` computes `π(x)`, the number of primes `≤ x`. The argument is floored to an integer (Real and Rational accepted). `init_primepi` (run once, cached) sieves all primes up to `MAX_PRIME_LIMIT = 10^6` with a sieve of Eratosthenes into `pi_primes[]`, and allocates the `phi_cache` (Legendre partial-sieve memo table, `CACHE_A = 100` rows × `CACHE_X = 20000` columns).

For `x ≤ 10^6`, `pi_base` answers directly by binary-searching the cached prime list. For larger `x`, `count_primes` uses the Meissel–Legendre identity `π(x) = φ(x, a) + a - 1 - P2(x, a)` with `a = π(x^{1/3})`: `phi_rec` evaluates the Legendre partial-sieve function `φ(x, a)` by the recurrence `φ(x, a) = φ(x, a-1) - φ(x/p_a, a-1)` (base cases `φ(x,0)=x`, `φ(x,1)=x - ⌊x/2⌋`), memoised in `phi_cache` for small `(a, x)`; `p2_calc` computes the two-prime-factor correction `P2(x, a) = Σ_{a<i≤π(√x)} (π(x/p_i) - i + 1)`.

**Data structures.** `pi_primes` is a flat `uint32_t[]` of all primes ≤ 10^6 (`total_pi_primes` of them); `phi_cache` is a 2-D `int32_t**` memo for the Legendre φ. All counting is in `uint64_t`; result returned as `EXPR_INTEGER`.

**Complexity / limits.** The Meissel–Legendre method is roughly `O(x^{2/3})`-ish in practice for the φ recursion plus the `P2` correction, far cheaper than a full sieve to `x`. The sieve floor is fixed at 10^6; the φ memo only covers `a < 100, x < 20000`.

- `Listable`, `Protected`.
- Uses Meissel-Lehmer algorithm for efficient counting.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- M. Deléglise and J. Rivat, "Computing π(x): the Meissel, Lehmer, Lagarias, Miller, Odlyzko method", Math. Comp. 65, 1996.
- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
