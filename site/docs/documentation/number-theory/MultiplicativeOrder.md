# MultiplicativeOrder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MultiplicativeOrder[k, n] gives the multiplicative order of k modulo n, the smallest positive integer m such that k^m is congruent to 1 modulo n; reversing it -- recovering the exponent x with k^x congruent to a given target -- is the discrete logarithm problem, whose presumed hardness underlies Diffie-Hellman key exchange and the ElGamal and DSA schemes.`**

**`MultiplicativeOrder[k, n, {r1, r2, ...}]`**

gives the smallest positive integer m such that k^m is congruent to one of the ri modulo n -- a multi-target discrete logarithm.

<details>
<summary>Notes</summary>

Computing the order is easy (it divides EulerPhi\[n\]); the discrete logarithm it inverts is not, so the three-argument form is a teaching tool for small moduli rather than a practical logarithm at cryptographic sizes.  Returns unevaluated when gcd(k, n) is not 1, when no power of k lands in the residue set, or when n is zero.  All arithmetic is exact via GMP, so k and n may be arbitrary-precision integers.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= MultiplicativeOrder[5, 8]
Out[1]= 2

In[2]:= MultiplicativeOrder[5, 7]
Out[2]= 6

In[3]:= MultiplicativeOrder[-5, 7]
Out[3]= 3

In[4]:= MultiplicativeOrder[5, 7, {3, 11}]
Out[4]= 2

In[5]:= MultiplicativeOrder[10^10000, 7919]
Out[5]= 3959

In[6]:= Select[Range[43], MultiplicativeOrder[#, 43] == EulerPhi[43] &]
Out[6]= {3, 5, 12, 18, 19, 20, 26, 28, 29, 30, 33, 34}

In[7]:= MultiplicativeOrder[10, 22]
Out[7]= MultiplicativeOrder[10, 22]
```

### Applications (5)

```mathematica
In[8]:= MultiplicativeOrder[2, 7]
Out[8]= 3

In[9]:= MultiplicativeOrder[10, 7]
Out[9]= 6

In[10]:= MultiplicativeOrder[7, 1000000007]
Out[10]= 500000003

In[11]:= MultiplicativeOrder[3, 998244353]
Out[11]= 998244352

In[12]:= MultiplicativeOrder[2, 11, {1, 10}]
Out[12]= 5
```

## Implementation notes

**Algorithm.** `builtin_multiplicativeorder` computes the multiplicative order of `k` modulo `n` — the least `m` with `k^m ≡ 1 (mod n)`. `mo_order_mpz` reduces `k mod n`, checks `gcd(k, n) = 1` (returning unevaluated otherwise), computes Euler's totient `φ(n)` (`mo_eulerphi_mpz`), and then deflates the order down from `φ(n)`: for each distinct prime `q | φ(n)` it repeatedly divides the running order by `q` while `k^(order/q) ≡ 1 (mod n)` (using `mpz_powm`). This yields the order without enumerating all exponents. The 3-argument form `MultiplicativeOrder[k, n, {r1, ...}]` instead searches for the least `m ≤ order` with `k^m` congruent to one of the residues `r_i` (`mo_search_residues`, capped at `MO_SEARCH_CAP = 10^8` iterations and requiring the order to fit in `unsigned long`). Wrong arg counts emit `MultiplicativeOrder::argt`.

**Data structures.** GMP `mpz_t` throughout; the distinct prime divisors of `φ(n)` are collected into a fixed `mpz_t primes[]` array (`pr_collect_distinct_primes`).

**Complexity / limits.** Order computation is `O(ω(φ(n)) · log φ(n))` modular exponentiations after factoring `φ(n)`; the residue-search form is bounded by `MO_SEARCH_CAP`.

- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be
  arbitrary-precision bignums.
- Negative `n` is treated as `|n|`; `k` (and each `r_i`) is reduced modulo
  `n` before the search, so negative or out-of-range inputs work
  transparently.
- Returns unevaluated when `gcd(k, n) != 1`, when `n` is zero, or when
  no power of `k` lands in the residue set (3-arg form).
- The 2-arg form factors `phi(n)` and successively strips prime factors
  from `phi(n)` whose corresponding exponent still maps `k` to 1 — so
  the work is dominated by factoring `phi(n)`, not by walking the orbit.
- The 3-arg form walks `k^m mod n` for `m = 1, ..., order(k, n)`. To
  guard against pathological group sizes, the call returns unevaluated
  if the order exceeds `10^8` or does not fit in an unsigned long.
- Non-integer numeric inputs (`Real`, `Complex`, `Rational`) and symbolic
  arguments flow through unevaluated with no diagnostic.
- Diagnostic: `MultiplicativeOrder::argt` when called with anything other
  than 2 or 3 arguments.

**Attributes:** `Protected`.

## References

**See also:** [Complex](../../arithmetic/Complex/), [Rational](../../arithmetic/Rational/)

- A. J. Menezes, P. C. van Oorschot and S. A. Vanstone, *Handbook of Applied Cryptography*, CRC Press, 1996 — §3.6 covers the discrete logarithm problem and its algorithms (baby-step giant-step, Pohlig–Hellman, index calculus).
- R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, 2nd ed., Springer, 2005 — §5.2, discrete logarithms.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_multiplicative_order.c`](https://github.com/stblake/mathilda/blob/main/tests/test_multiplicative_order.c)

## Notes & additional examples

### The discrete logarithm problem

The multiplicative order of `k` modulo `n` is the size of the cyclic subgroup that `k`
generates in the units mod `n`. It is intimately tied to the **discrete logarithm problem**
(DLP): given a base `k` and a target `b`, find the exponent `x` with `k^x ≡ b (mod n)`. The
three-argument form `MultiplicativeOrder[k, n, {r1, ...}]` returns the least exponent whose
power lands in the residue set, so it *is* a (multi-target) discrete logarithm:

```mathematica
In[1]:= MultiplicativeOrder[3, 7, {5}]
Out[1]= 5
```

Here `3` is a primitive root of `7`, and the least `x` with `3^x ≡ 5 (mod 7)` is `5` — the
discrete logarithm of `5` to base `3`. Computing the *order* is easy: it divides `EulerPhi[n]`
and is recovered by stripping prime factors from it. Inverting it — the DLP itself — is
believed **hard**: no polynomial-time algorithm is known over a general prime field, and that
presumed hardness is the foundation of Diffie–Hellman key exchange and the ElGamal and DSA
signature schemes. The standard algorithms (baby-step giant-step, Pollard's rho for
logarithms, Pohlig–Hellman, and index calculus) are surveyed in the references below. Mathilda's
three-argument form simply walks the orbit `k^1, k^2, …`, so it is a teaching tool for small
moduli, not a practical logarithm at cryptographic sizes.

### Notes

`MultiplicativeOrder[k, n]` is the smallest `m > 0` with `k^m ≡ 1 (mod n)`. The
order `6` for `10` modulo `7` reflects that `1/7 = 0.142857...` has a repeating
block of length `6`. The two large-modulus cases use prime moduli: `3` is a
primitive root of the NTT prime `998244353`, so its order equals `n - 1`. The
three-argument form `MultiplicativeOrder[k, n, {r1, ...}]` finds the least `m`
with `k^m` congruent to one of the listed residues. All arithmetic is exact via
GMP. The result is unevaluated when `gcd(k, n) ≠ 1`.
