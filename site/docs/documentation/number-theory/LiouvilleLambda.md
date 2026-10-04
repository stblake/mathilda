# LiouvilleLambda

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LiouvilleLambda[n] gives the Liouville function lambda(n) = (-1)^Omega(n), where Omega(n) counts the prime factors of n with multiplicity. Completely multiplicative. A non-real Gaussian-integer argument, or GaussianIntegers -> True, is handled over Z[i].`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= LiouvilleLambda[8]
Out[1]= -1

In[2]:= LiouvilleLambda[9]
Out[2]= 1

In[3]:= LiouvilleLambda[{1, 2, 3, 4, 5, 6}]
Out[3]= {1, -1, -1, 1, -1, 1}

In[4]:= LiouvilleLambda[10^30 + 1]
Out[4]= -1

In[5]:= LiouvilleLambda[2 + I]
Out[5]= -1
```

### Options (1)

```mathematica
In[6]:= LiouvilleLambda[8, GaussianIntegers -> True]
Out[6]= 1
```

### Applications (6)

Lambda(n) = (-1)^Omega(n); 12 = 2^2 * 3 has Omega = 3

```mathematica
In[7]:= LiouvilleLambda[12]
Out[7]= -1
```

With 210 = 2 * 3 * 5 * 7, Omega = 4, so lambda = 1

```mathematica
In[8]:= LiouvilleLambda[210]
Out[8]= 1
```

The empty product is 1

```mathematica
In[9]:= LiouvilleLambda[1]
Out[9]= 1
```

Threaded element-wise over a range by Listable

```mathematica
In[10]:= LiouvilleLambda[Range[10]]
Out[10]= {1, -1, -1, 1, -1, 1, -1, -1, 1, 1}
```

Complete multiplicativity, with no coprimality needed

```mathematica
In[11]:= LiouvilleLambda[6 * 35] == LiouvilleLambda[6] LiouvilleLambda[35]
Out[11]= True
```

Summing lambda over the divisors of n is 1 exactly when n is a perfect square

```mathematica
In[12]:= Table[Total[LiouvilleLambda[Divisors[n]]], {n, 1, 10}]
Out[12]= {1, 0, 0, 1, 0, 0, 0, 0, 1, 0}
```

## Implementation notes

**Algorithm.** `builtin_liouvillelambda` factors `|n|` and returns
`lambda(n) = (-1)^Omega(n)`, where `Omega(n)` is the sum of the prime-factor exponents
(`liouville_from_exps`). It takes one positional argument plus an optional `GaussianIntegers`
rule; under `GaussianIntegers -> True`, or for a non-real Gaussian `n`, the count runs over
the Gaussian prime factorisation. lambda is completely multiplicative, with
`lambda(-n) = lambda(n)` (sign ignored) and `lambda(1) = 1` (the empty product).

**Data structures.** The factorisation is a `mpz_t*` / `unsigned long*` pair from
`df_factor_mpz` (ordinary) or `df_gaussian_prime_factor` (Gaussian); only the exponent sum
is needed, and the result is a machine `Integer` of `+1` or `-1`. Unlike its sibling
`MoebiusMu`, `LiouvilleLambda` has no dedicated ND/packed/`Compile` kernel; `Listable`
threading is handled by the evaluator, so `LiouvilleLambda[Range[n]]` maps element-wise over
the boxed list.

**Complexity / limits.** Dominated by factoring `|n|` through `FactorInteger` (trial
division + Pollard rho + ECM). `df_factor_mpz` confirms each base prime with 40 Miller–Rabin
rounds and declines — leaving the call unevaluated — rather than return a confident wrong
value on an unfactored composite cofactor (the failure that once made these counting
functions wrong on an 82-digit semiprime). `n == 0` is left unevaluated, and a wrong
argument count emits `LiouvilleLambda::argt`.

- `Listable`, `Protected`.
- Completely multiplicative: `lambda(m n) = lambda(m) lambda(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `1` or `-1`.
- The sign of `n` is ignored (`lambda(-n) = lambda(n)`).
- Gaussian integers: `LiouvilleLambda[n, GaussianIntegers -> True]`, or a
  non-real Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]`
  and counts the Gaussian prime factors with multiplicity. Because `2` factors
  as `-i (1 + i)^2` in `Z[i]`, e.g. `LiouvilleLambda[2, GaussianIntegers -> True]`
  is `1` while `LiouvilleLambda[2]` is `-1`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `LiouvilleLambda::argt` message.

**Attributes:** `Listable`, `Protected`.

## References

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Liouville function lambda and its summatory divisor identity (Chapter 17, 22).
- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Liouville function lambda and the prime-factor count Omega (Chapter 17).
- Source: [`src/numbertheory/liouvillelambda.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/liouvillelambda.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_liouvillelambda.c`](https://github.com/stblake/mathilda/blob/main/tests/test_liouvillelambda.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)

## Notes & additional examples

### Notes

The Liouville function `LiouvilleLambda[n] = (-1)^Omega(n)` is `+1` when `n` has an even
number of prime factors counted with multiplicity and `-1` when odd, so it equals
`(-1)^PrimeOmega[n]` (see [`PrimeOmega`](PrimeOmega.md)). It is *completely* multiplicative —
`lambda(ab) = lambda(a) lambda(b)` for all `a, b`, unlike [`MoebiusMu`](MoebiusMu.md), which
is only multiplicative over coprime arguments — and the sign of `n` is ignored.

Its summatory divisor identity `sum_{d | n} lambda(d) = [n is a perfect square]` is the
analogue of the Moebius identity and detects squares. A non-real Gaussian-integer argument,
or `GaussianIntegers -> True`, counts Gaussian prime factors over `Z[i]`.
`LiouvilleLambda[0]` is left unevaluated.
