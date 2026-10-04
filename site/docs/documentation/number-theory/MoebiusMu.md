# MoebiusMu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MoebiusMu[n] gives the Moebius function mu(n): 0 if n has a squared prime factor, otherwise (-1)^k where k is the number of distinct primes. A non-real Gaussian-integer argument is handled over Z[i].`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= MoebiusMu[11]
Out[1]= -1

In[2]:= MoebiusMu[10]
Out[2]= 1

In[3]:= MoebiusMu[1440]
Out[3]= 0

In[4]:= MoebiusMu[{4, 10, 17, 20}]
Out[4]= {0, 1, -1, 0}

In[5]:= MoebiusMu[10^50 + 1]
Out[5]= MoebiusMu[100000000000000000000000000000000000000000000000001]

In[6]:= MoebiusMu[5 + 6 I]
Out[6]= -1
```

### Applications (6)

A square-free integer with three prime factors gives (-1)^3

```mathematica
In[7]:= MoebiusMu[30]
Out[7]= -1
```

Any squared prime factor sends the value to 0

```mathematica
In[8]:= MoebiusMu[12]
Out[8]= 0
```

The empty product is 1

```mathematica
In[9]:= MoebiusMu[1]
Out[9]= 1
```

Threaded element-wise over a range by Listable

```mathematica
In[10]:= MoebiusMu[Range[12]]
Out[10]= {1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0}
```

Summing mu over the divisors of n is 1 only at n = 1 -- the identity behind Moebius inversion

```mathematica
In[11]:= Table[Total[MoebiusMu[Divisors[n]]], {n, 1, 8}]
Out[11]= {1, 0, 0, 0, 0, 0, 0, 0}
```

A non-real argument is factored over the Gaussian integers

```mathematica
In[12]:= MoebiusMu[1 + I]
Out[12]= -1
```

## Options & behaviour

> **Packed arrays.** Runs on an `int64` buffer. `MoebiusMu[0]` is undefined,
> so an array containing `0` takes the ordinary path and leaves that element
> unevaluated exactly as the scalar does.

## Implementation notes

**Algorithm.** `builtin_moebiusmu` takes exactly one argument, factors `|n|`, and returns
`mu(n)`: `0` if any exponent is `>= 2` (a squared prime factor), otherwise `(-1)^m` for `m`
distinct primes, with `mu(1) = 1` (`moebiusmu_from_exps`). A non-real Gaussian-integer
argument is auto-detected and factored over `Z[i]` (the unit factor does not count). The
sign of `n` is ignored, matching `mu(-n) = mu(n)`.

**Data structures.** Factorisation via `df_factor_mpz` / `df_gaussian_prime_factor` into
`mpz_t` / `unsigned long` arrays; the answer is a machine `Integer` in `{-1, 0, 1}`. An
exact-integer ND kernel exists: `ndk_MoebiusMu_ii` in `src/ndinteger.c` (registered with
`symtab_set_ndarray_unary_kernel`, placing `MoebiusMu` on both the `AWARE` and `INT64_OK`
lists in `src/pack.c`) factors each `int64` element in place and declines `n == 0` or a
factor past the trial-division ceiling back to the GMP List path. `Compile[]` lowers
`MoebiusMu[v]` at a rank-1 integer-array shape (`Compiled -> True`) but not at a scalar
shape — the integer-only kernel is a real path over an array and no path at all over a
scalar (the `narrowing_only` branch in `src/compile/compile.c`).

**Complexity / limits.** Dominated by factoring `|n|`. As with its siblings, `df_factor_mpz`
verifies each base's primality with 40 Miller–Rabin rounds and declines rather than trust a
composite cofactor — which is what stopped `MoebiusMu` of an 82-digit semiprime returning a
confident wrong `-1`. `MoebiusMu[0]` is left unevaluated, and a wrong argument count emits
`MoebiusMu::argx`.

- `Listable`, `Protected`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `0`, `1`, or `-1`.
- The sign of `n` is ignored (`mu(-n) = mu(n)`).
- A non-real Gaussian-integer argument `Complex[a, b]` is handled over `Z[i]`:
  the input is factored into Gaussian primes (the unit factor does not count),
  giving `0` for a repeated Gaussian prime factor and `(-1)^m` otherwise.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `MoebiusMu::argx` message.

**Attributes:** `Listable`, `Protected`.

## References

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Moebius function and Moebius inversion (Chapter 16).
- Source: [`src/numbertheory/moebiusmu.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/moebiusmu.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_moebiusmu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_moebiusmu.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)

## Notes & additional examples

### Notes

`MoebiusMu[n]` is `0` when `n` has a squared prime factor, and otherwise `(-1)^k` for `k`
distinct primes, so it detects square-freeness (`MoebiusMu[n] != 0` iff `n` is square-free,
the test [`SquareFreeQ`](SquareFreeQ.md) makes directly). The sign of `n` is ignored. The
summatory identity `sum_{d | n} mu(d) = [n == 1]` is the kernel of Moebius inversion, which
recovers an arithmetic function from its divisor sums.

`MoebiusMu` is the only member of its factor-counting family — alongside
[`PrimeNu`](PrimeNu.md), [`PrimeOmega`](PrimeOmega.md) and
[`LiouvilleLambda`](LiouvilleLambda.md) — with a machine-integer array kernel, so a packed
or `NDArray` `int64` argument is factored element-wise on the buffer, and the head also
lowers inside `Compile[]` at a rank-1 integer-array shape. A non-real Gaussian-integer
argument, or `GaussianIntegers -> True`, switches the factorisation to `Z[i]`. `MoebiusMu[0]`
is left unevaluated.
