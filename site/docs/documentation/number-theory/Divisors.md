# Divisors

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Divisors[n] gives a list of the integers that divide n. Divisors[n, GaussianIntegers -> True] includes Gaussian-integer divisors.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Divisors[1729]
Out[1]= {1, 7, 13, 19, 91, 133, 247, 1729}

In[2]:= Divisors[6]
Out[2]= {1, 2, 3, 6}

In[3]:= Divisors[{605, 871, 824}]
Out[3]= {{1, 5, 11, 55, 121, 605}, {1, 13, 67, 871}, {1, 2, 4, 8, 103, 206, 412, 824}}

In[4]:= Divisors[6 + 4 I]
Out[4]= {1, 1 + I, 1 + 5*I, 2, 3 + 2*I, 6 + 4*I}
```

### Options (2)

```mathematica
In[5]:= Divisors[2, GaussianIntegers -> True]
Out[5]= {1, 1 + I, 2}

In[6]:= Divisors[3, GaussianIntegers -> True]
Out[6]= {1, 3}
```

### Applications (6)

A perfect number equals the sum of its proper divisors

```mathematica
In[7]:= Divisors[28]
Out[7]= {1, 2, 4, 7, 14, 28}
```

The twelve divisors, ascending

```mathematica
In[8]:= Divisors[60]
Out[8]= {1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60}
```

Summing all divisors of a perfect number gives twice the number

```mathematica
In[9]:= Total[Divisors[28]]
Out[9]= 56
```

The number of divisors is sigma_0

```mathematica
In[10]:= Length[Divisors[720]] == DivisorSigma[0, 720]
Out[10]= True
```

Threaded over a list by Listable

```mathematica
In[11]:= Divisors[{12, 15}]
Out[11]= {{1, 2, 3, 4, 6, 12}, {1, 3, 5, 15}}
```

Over Z[i]: one first-quadrant representative per associate class

```mathematica
In[12]:= Divisors[5, GaussianIntegers -> True]
Out[12]= {1, 1 + 2*I, 2 + I, 5}
```

## Implementation notes

**Algorithm.** `builtin_divisors` separates the single positional argument from an optional
`GaussianIntegers` rule, then enumerates divisors from the prime factorisation. The ordinary
path (`divisors_ordinary` in `nt_gaussian.c`) factors `|n|`, walks the divisor lattice with
a mixed-radix exponent odometer (digit `i` ranges `0..e_i`, forming each divisor as a
product of prime powers), and `qsort`s the results into ascending order; `Divisors[1]` is
`{1}`. The Gaussian path (`divisors_gaussian`), used under `GaussianIntegers -> True` or for
a non-real input, returns one first-quadrant representative per associate class, sorted by
`(Re, Im)`.

**Data structures.** Divisors are built in a `mpz_t` array via `mpz_pow_ui`/`mpz_mul` over
the prime powers, then emitted as a `List` of `Integer`/`BigInt` (ordinary) or `Complex`
(Gaussian) `Expr`. There is no ND/packed/`Compile` kernel — the result is a
variable-length list, not a machine buffer; `Listable` threading over a list of arguments is
done by the evaluator.

**Complexity / limits.** Cost is factoring plus `O(d log d)` to sort the `d = prod_i (e_i +
1)` divisors. The divisor count is computed first, and the call is left unevaluated if it
overflows `size_t` (e.g. `Divisors[100!]` has ~10^28 divisors, intractable to materialise).
`Divisors[0]`, a non-integer `n`, and a factorisation whose bases cannot be confirmed prime
are all left unevaluated; `Divisors[]` emits `Divisors::argx`. The sign of `n` is ignored.

- `Listable`, `Protected`.
- Machine integers and GMP bigints are handled uniformly; the result promotes to
  a big-integer list when needed.
- The sign of `n` is ignored (`Divisors[-12] == Divisors[12]`).
- Divisors are computed from the prime factorization (the divisor lattice), so
  cost scales with the number of divisors rather than `Sqrt[n]`.
- In Gaussian mode each divisor is the canonical first-quadrant representative of
  its associate class (`Re > 0`, `Im >= 0`), and the list is sorted by
  `(Re, Im)`. Rational primes are lifted to `Z[i]`: `2` ramifies as `1 + I`,
  primes `p ≡ 1 (mod 4)` split via a sum-of-two-squares (Cornacchia)
  decomposition, and primes `p ≡ 3 (mod 4)` stay inert.
- `Divisors[0]`, non-integer arguments, and calls whose divisor count would
  overflow (e.g. `Divisors[100!]`) are left unevaluated; `Divisors[]` issues a
  `Divisors::argx` message.

**Attributes:** `Listable`, `Protected`.

## References

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — divisors and perfect numbers (Chapter 16).
- Source: [`src/numbertheory/divisors.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/divisors.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_intersection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intersection.c)

## Notes & additional examples

### Notes

`Divisors[n]` returns the ascending list of positive integers dividing `n`, built from the
prime factorisation (the divisor lattice), so the sign of `n` is ignored and `Divisors[1]`
is `{1}`. The count of divisors is [`DivisorSigma`](DivisorSigma.md)`[0, n]` and their sum
is `DivisorSigma[1, n]`; a number is perfect when that sum is `2 n`, as with `28`.

`Divisors[n, GaussianIntegers -> True]`, or a non-real Gaussian-integer `n`, returns the
divisors in `Z[i]`, one representative per unit-associate class sorted by `(Re, Im)`.
`Divisors[0]`, a non-integer argument, and any call whose divisor count would overflow (such
as `Divisors[100!]`) are left unevaluated.
