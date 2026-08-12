# PrimeNu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimeNu[n] gives the number of distinct prime factors of n, nu(n). PrimeNu[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts distinct Gaussian prime factors over Z[i]. PrimeNu[1] is 0; PrimeNu[0] is left unevaluated.`**

## Examples (6)

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

## See also

[PrimeOmega](../../number-theory/PrimeOmega/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
