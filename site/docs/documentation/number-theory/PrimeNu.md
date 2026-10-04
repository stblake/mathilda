# PrimeNu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimeNu[n] gives the number of distinct prime factors of n, nu(n). PrimeNu[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts distinct Gaussian prime factors over Z[i]. PrimeNu[1] is 0; PrimeNu[0] is left unevaluated.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= PrimeNu[24]
Out[1]= 2

In[2]:= PrimeNu[105]
Out[2]= 3

In[3]:= PrimeNu[{4, 28, 180}]
Out[3]= {1, 2, 3}

In[4]:= PrimeNu[50!]
Out[4]= 15

In[5]:= PrimeNu[3 + I]
Out[5]= 2
```

### Options (1)

```mathematica
In[6]:= PrimeNu[105, GaussianIntegers -> True]
Out[6]= 4
```

### Applications (5)

Distinct prime factors {2, 3, 5}, regardless of exponents

```mathematica
In[7]:= PrimeNu[360]
Out[7]= 3
```

A prime power has a single distinct factor

```mathematica
In[8]:= PrimeNu[2^10]
Out[8]= 1
```

No prime factors

```mathematica
In[9]:= PrimeNu[1]
Out[9]= 0
```

Threaded element-wise over a range by Listable

```mathematica
In[10]:= PrimeNu[Range[10]]
Out[10]= {0, 1, 1, 1, 1, 2, 1, 1, 1, 2}
```

Equal to PrimeOmega only when n is square-free; 72 = 2^3 * 3^2 is not

```mathematica
In[11]:= PrimeNu[72]
Out[11]= 2
```

## Algorithm

primenu.c -- PrimeNu[]. Split from numbertheory.c; see numbertheory.h and numbertheory_internal.h for the subsystem layout.

```text
PrimeNu[n] = nu(n), the number of DISTINCT prime factors of n.  It is the
```

additive companion to PrimeOmega (which counts prime factors with multiplicity): for n = u p_1^k_1 ... p_m^k_m with u a unit and p_i distinct

```text
primes, PrimeNu[n] returns m.  nu and Omega coincide exactly when n is
square-free.  PrimeNu shares all factoring machinery and argument handling
```

with PrimeOmega/LiouvilleLambda; it simply returns the count of factors rather than the sum of the exponents.

## Implementation notes

**Algorithm.** `builtin_primenu` factors `|n|` and returns `nu(n)`, the number of distinct
prime factors — the count of factors, independent of their exponents (`primenu_from_count`).
It is the additive companion to `PrimeOmega` (which sums the exponents); `nu` and `Omega`
coincide exactly when `n` is square-free. It takes one positional argument plus an optional
`GaussianIntegers` rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`,
the count runs over the distinct Gaussian prime factors. `nu(1) = nu(-1) = 0`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; only the factor count is used, and the result is a machine
`Integer`. Unlike `MoebiusMu`, `PrimeNu` has no dedicated ND/packed/`Compile` kernel;
`Listable` threading is handled by the evaluator.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger`;
`df_factor_mpz` confirms each base prime with 40 Miller–Rabin rounds and declines — leaving
the call unevaluated — on an unfactored composite cofactor rather than returning a wrong
count. `n == 0` is left unevaluated, and a wrong argument count emits `PrimeNu::argt`.

- `Listable`, `Protected`.
- Additive on coprime arguments: `nu(m n) = nu(m) + nu(n)` when
  `GCD[m, n] == 1`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeNu[1]` (and `PrimeNu[-1]`) is `0`; the sign of `n` is ignored
  (`nu(-n) = nu(n)`).
- Gaussian integers: `PrimeNu[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the distinct Gaussian prime factors. Because a rational prime `p ≡ 1 (mod 4)`
  splits into two conjugate Gaussian primes, e.g.
  `PrimeNu[105, GaussianIntegers -> True]` is `4` (from `3`, the split pair over
  `5`, and `7`) while `PrimeNu[105]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeNu::argt` message.
- Relations: for a square-free `n`, `MoebiusMu[n] == (-1)^PrimeNu[n]` and
  `LiouvilleLambda[n] == (-1)^PrimeNu[n]`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [PrimeOmega](../../number-theory/PrimeOmega/)

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the prime-factor counting functions nu and Omega (Chapter 22).
- Source: [`src/numbertheory/primenu.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/primenu.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_moebiusmu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_moebiusmu.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)

## Notes & additional examples

### Notes

`PrimeNu[n]` is the number of *distinct* prime factors of `n`, independent of their
multiplicities. Its companion [`PrimeOmega`](PrimeOmega.md) counts prime factors *with*
multiplicity, so `PrimeNu[n] <= PrimeOmega[n]`, with equality exactly when `n` is
square-free (the condition [`MoebiusMu`](MoebiusMu.md) tests). The sign of `n` is ignored
and `PrimeNu[1] = 0`.

A non-real Gaussian-integer argument, or `GaussianIntegers -> True`, counts distinct Gaussian
prime factors over `Z[i]`. `PrimeNu[0]` is left unevaluated.
