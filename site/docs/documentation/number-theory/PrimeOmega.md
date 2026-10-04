# PrimeOmega

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimeOmega[n] gives the number of prime factors of n counted with multiplicity, Omega(n). PrimeOmega[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts Gaussian prime factors over Z[i]. PrimeOmega[1] is 0; PrimeOmega[0] is left unevaluated.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= PrimeOmega[30]
Out[1]= 3

In[2]:= PrimeOmega[12]
Out[2]= 3

In[3]:= PrimeOmega[{4, 12, 24}]
Out[3]= {2, 3, 4}

In[4]:= PrimeOmega[30!]
Out[4]= 59

In[5]:= PrimeOmega[5 + 9 I]
Out[5]= 2
```

### Options (1)

```mathematica
In[6]:= PrimeOmega[12, GaussianIntegers -> True]
Out[6]= 5
```

### Applications (5)

Counted with multiplicity: 360 = 2^3 * 3^2 * 5 gives 3 + 2 + 1

```mathematica
In[7]:= PrimeOmega[360]
Out[7]= 6
```

A prime power contributes its full exponent

```mathematica
In[8]:= PrimeOmega[2^10]
Out[8]= 10
```

The empty factorisation gives 0

```mathematica
In[9]:= PrimeOmega[1]
Out[9]= 0
```

Threaded element-wise over a range by Listable

```mathematica
In[10]:= PrimeOmega[Range[10]]
Out[10]= {0, 1, 1, 2, 1, 2, 1, 3, 2, 2}
```

Compare with PrimeNu[72] = 2; 72 = 2^3 * 3^2 is not square-free, so Omega is larger

```mathematica
In[11]:= PrimeOmega[72]
Out[11]= 5
```

## Algorithm

primeomega.c -- PrimeOmega[]. Split from numbertheory.c; see numbertheory.h and numbertheory_internal.h for the subsystem layout.

PrimeOmega[n] = Omega(n), the number of prime factors of n counted with

```text
multiplicity (the sum of the exponents in the prime factorization).  This is
```

the quantity LiouvilleLambda computes internally before taking (-1)^Omega, so the two share the same factoring machinery and argument handling; PrimeOmega simply returns Omega itself.

## Implementation notes

**Algorithm.** `builtin_primeomega` factors `|n|` and returns `Omega(n)`, the number of
prime factors counted with multiplicity — the sum of the exponents in the prime
factorisation (`primeomega_from_exps`). This is the quantity `LiouvilleLambda` forms before
taking `(-1)^Omega`. It takes one positional argument plus an optional `GaussianIntegers`
rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`, the count runs over
the Gaussian prime factors. `Omega(1) = Omega(-1) = 0`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; the exponent sum gives a machine `Integer`. `PrimeOmega`
has no dedicated ND/packed/`Compile` kernel; `Listable` threading is handled by the
evaluator.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger`;
`df_factor_mpz` verifies each base's primality with 40 Miller–Rabin rounds and declines —
leaving the call unevaluated — rather than trust a composite cofactor. `n == 0` is left
unevaluated, and a wrong argument count emits `PrimeOmega::argt`.

- `Listable`, `Protected`.
- Completely additive: `Omega(m n) = Omega(m) + Omega(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeOmega[1]` (and `PrimeOmega[-1]`) is `0`; the sign of `n` is ignored
  (`Omega(-n) = Omega(n)`).
- Gaussian integers: `PrimeOmega[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the Gaussian prime factors with multiplicity. Because `2` factors as
  `-i (1 + i)^2` in `Z[i]`, `PrimeOmega[12, GaussianIntegers -> True]` is `5`
  (from `(1 + i)^4 3`) while `PrimeOmega[12]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeOmega::argt` message.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [LiouvilleLambda](../../number-theory/LiouvilleLambda/)

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the prime-factor counting functions nu and Omega (Chapter 22).
- Source: [`src/numbertheory/primeomega.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/primeomega.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_moebiusmu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_moebiusmu.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
- Tests: [`tests/test_primeomega.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primeomega.c)

## Notes & additional examples

### Notes

`PrimeOmega[n]` is the total number of prime factors of `n` counted with multiplicity — the
sum of the exponents in the prime factorisation. Its additive companion
[`PrimeNu`](PrimeNu.md) counts *distinct* primes, and the two agree exactly when `n` is
square-free; in general `PrimeOmega[n] >= PrimeNu[n]`. `PrimeOmega` is the exponent
[`LiouvilleLambda`](LiouvilleLambda.md) raises `-1` to, so `LiouvilleLambda[n] ==
(-1)^PrimeOmega[n]`.

The sign of `n` is ignored and `PrimeOmega[1] = 0`. A non-real Gaussian-integer argument, or
`GaussianIntegers -> True`, counts Gaussian prime factors over `Z[i]`. `PrimeOmega[0]` is
left unevaluated.
