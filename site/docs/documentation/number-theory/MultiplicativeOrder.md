# MultiplicativeOrder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MultiplicativeOrder[k, n]
    gives the multiplicative order of k modulo n, the smallest positive integer m such that k^m is congruent to 1 modulo n.
MultiplicativeOrder[k, n, {r1, r2, ...}]
    gives the smallest positive integer m such that k^m is congruent to one of the ri modulo n.

Returns unevaluated when gcd(k, n) is not 1, when no power of k lands in the residue set, or when n is zero.  All arithmetic is exact via GMP, so k and n may be arbitrary-precision integers.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Algorithm.** `builtin_multiplicativeorder` computes the multiplicative order of `k` modulo `n` — the least `m` with `k^m ≡ 1 (mod n)`. `mo_order_mpz` reduces `k mod n`, checks `gcd(k, n) = 1` (returning unevaluated otherwise), computes Euler's totient `φ(n)` (`mo_eulerphi_mpz`), and then deflates the order down from `φ(n)`: for each distinct prime `q | φ(n)` it repeatedly divides the running order by `q` while `k^(order/q) ≡ 1 (mod n)` (using `mpz_powm`). This yields the order without enumerating all exponents. The 3-argument form `MultiplicativeOrder[k, n, {r1, ...}]` instead searches for the least `m ≤ order` with `k^m` congruent to one of the residues `r_i` (`mo_search_residues`, capped at `MO_SEARCH_CAP = 10^8` iterations and requiring the order to fit in `unsigned long`). Wrong arg counts emit `MultiplicativeOrder::argt`.

**Data structures.** GMP `mpz_t` throughout; the distinct prime divisors of `φ(n)` are collected into a fixed `mpz_t primes[]` array (`pr_collect_distinct_primes`).

**Complexity / limits.** Order computation is `O(ω(φ(n)) · log φ(n))` modular exponentiations after factoring `φ(n)`; the residue-search form is bounded by `MO_SEARCH_CAP`.

- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
