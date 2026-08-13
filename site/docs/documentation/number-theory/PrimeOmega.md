# PrimeOmega

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimeOmega[n] gives the number of prime factors of n counted with multiplicity, Omega(n). PrimeOmega[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts Gaussian prime factors over Z[i]. PrimeOmega[1] is 0; PrimeOmega[0] is left unevaluated.`**

## Examples (6)

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

## Algorithm

primeomega.c -- PrimeOmega[]. Split from numbertheory.c; see numbertheory.h and numbertheory_internal.h for the subsystem layout.

PrimeOmega[n] = Omega(n), the number of prime factors of n counted with

```text
multiplicity (the sum of the exponents in the prime factorization).  This is
```

the quantity LiouvilleLambda computes internally before taking (-1)^Omega, so the two share the same factoring machinery and argument handling; PrimeOmega simply returns Omega itself.

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
- Tests: [`tests/test_primeomega.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primeomega.c)
