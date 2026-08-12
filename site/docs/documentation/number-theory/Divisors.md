# Divisors

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Divisors[n] gives a list of the integers that divide n. Divisors[n, GaussianIntegers -> True] includes Gaussian-integer divisors.`**

## Examples (6)

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

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_intersection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intersection.c)
