# DivisorSigma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DivisorSigma[k, n] gives the divisor function sigma_k(n), the sum of the k-th powers of the divisors of n. DivisorSigma[k, n, GaussianIntegers -> True] sums over Gaussian-integer divisors.`**

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= DivisorSigma[1, 20]
Out[1]= 42

In[2]:= DivisorSigma[2, 20]
Out[2]= 546

In[3]:= DivisorSigma[0, 12]
Out[3]= 6

In[4]:= DivisorSigma[-2, 10]
Out[4]= 13/10

In[5]:= DivisorSigma[1/2, 12]
Out[5]= (2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))

In[6]:= DivisorSigma[k, {2, 3, 6}]
Out[6]= {(-1 + 2^(2 k))/(-1 + 2^k), (-1 + 3^(2 k))/(-1 + 3^k), ((-1 + 2^(2 k)) (-1 + 3^(2 k)))/((-1 + 2^k) (-1 + 3^k))}

In[7]:= DivisorSigma[2, {1, 2, 3, 4, 5}]
Out[7]= {1, 5, 10, 21, 26}

In[8]:= DivisorSigma[1, 3 + I]
Out[8]= 2 + 6*I
```

### Options (1)

```mathematica
In[9]:= DivisorSigma[2, 6, GaussianIntegers -> True]
Out[9]= -30 + 20*I
```

### Applications (6)

The sum of divisors; 28 is perfect, so sigma_1 = 2 * 28

```mathematica
In[10]:= DivisorSigma[1, 28]
Out[10]= 56
```

With k = 0, sigma_0 is the number of divisors

```mathematica
In[11]:= DivisorSigma[0, 60]
Out[11]= 12
```

Sum of squared divisors: 1 + 4 + 25 + 100

```mathematica
In[12]:= DivisorSigma[2, 10]
Out[12]= 130
```

Symbolic k keeps the multiplicative product in closed form

```mathematica
In[13]:= DivisorSigma[k, 12]
Out[13]= ((-1 + 2^(3 k)) (-1 + 3^(2 k)))/((-1 + 2^k) (-1 + 3^k))
```

Threaded element-wise over a range by Listable

```mathematica
In[14]:= DivisorSigma[1, Range[6]]
Out[14]= {1, 3, 4, 7, 6, 12}
```

Summed over the Gaussian-integer divisors

```mathematica
In[15]:= DivisorSigma[1, 5, GaussianIntegers -> True]
Out[15]= 4 + 8*I
```

## Options & behaviour

> **Packed arrays.** `DivisorSigma[k, list]` over an `int64` buffer factors
> each element by trial division in `int64`, with no GMP allocation per
> element. A non-negative integer `k` only: `DivisorSigma[-1, n]` is a
> `Rational`, which no buffer holds.

## Implementation notes

**Algorithm.** `builtin_divisorsigma` takes `k` and `n` (plus an optional
`GaussianIntegers` rule), factors `|n|` via `df_factor_mpz` (which delegates to
`FactorInteger`), and builds the multiplicative formula
`sigma_k(n) = prod_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)` as an `Expr` tree
(`ds_build_factor`) that it then evaluates. One path serves integer, rational, radical and
fully symbolic `k` — `DivisorSigma[k, 12]` returns the product in closed form in `k`. The
special case `k == 0` degenerates to the divisor count `prod_i (e_i + 1)`
(`ds_divisor_count`). With `GaussianIntegers -> True`, or a non-real Gaussian `n`, it
factors into Gaussian primes (`df_gaussian_prime_factor`), normalises each prime to its
first-quadrant associate, and runs the same product.

**Data structures.** Primes and exponents come back as parallel `mpz_t*` / `unsigned long*`
arrays; prime atoms become `Expr` (`Integer`/`BigInt`, or `Complex[u, v]` for a Gaussian
prime), and the multiplicative factors are assembled as `Power`/`Plus`/`Times` trees and
reduced through `eval_and_free`. An exact-integer ND kernel exists: `ndk_DivisorSigma_ii`
in `src/ndinteger.c` (registered with `symtab_set_ndarray_binary_kernel`, so `DivisorSigma`
sits on both the `AWARE` and `INT64_OK` lists in `src/pack.c`) evaluates
`prod_i (1 + p^k + ... + p^(e_i k))` entirely in `int64` for a packed or visible `int64`
array, declining `n == 0`, `k < 0` (sigma_-1 is a `Rational`, which no int64 buffer holds)
and any overflow back to the GMP List path. `Compile[]` lowers `DivisorSigma[k, v]` at a
rank-1 integer-array shape (`Compiled -> True`) but not at a scalar shape — the
integer-only kernel has no scalar opcode.

**Complexity / limits.** Dominated by factoring `|n|` (trial division + Pollard rho + ECM
through `FactorInteger`); `df_factor_mpz` confirms every base prime with 40 Miller–Rabin
rounds and declines — leaving the call unevaluated — rather than trust a composite cofactor
on a hard semiprime. Sign of `n` is ignored, `n == 0` is left unevaluated, and a wrong
argument count emits `DivisorSigma::argrx`.

- `Listable`, `NHoldAll`, `Protected`.
- Computed from the multiplicative formula
  `sigma_k(n) = Product_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)` for
  `n = Product_i p_i^e_i`, so a single path serves every exponent type: exact
  integers and rationals for integer `k`, and symbolic / radical forms for
  symbolic or rational `k`. `k == 0` returns the divisor count `sigma_0(n)`.
- The sign of `n` is ignored; machine integers and GMP bigints are handled
  uniformly.
- In Gaussian mode the product runs over the first-quadrant associates
  (`Re > 0`, `Im >= 0`) of the Gaussian prime factors of `n`. This is the
  multiplicative definition — note it differs from naively summing
  `d^k` over `Divisors[n, GaussianIntegers -> True]`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `DivisorSigma::argrx` message.

**Attributes:** `Listable`, `NHoldAll`, `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/)

- T. M. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976 — the divisor function sigma_k and its multiplicative formula (Chapter 2).
- T. M. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976 — multiplicative functions and the divisor function sigma_k (Chapter 2).
- Source: [`src/numbertheory/divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/divisorsigma.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)

## Notes & additional examples

### Notes

`DivisorSigma[k, n]` is `sum_{d | n} d^k`. Mathilda evaluates it from the multiplicative
formula `sigma_k(n) = prod_i (p_i^((e_i+1) k) - 1)/(p_i^k - 1)` built directly from the prime
factorisation, so a single path serves an integer, rational, radical or fully symbolic `k`:
`DivisorSigma[k, 12]` returns the product unevaluated in `k`. The case `k = 0` is the divisor
count [`DivisorSigma`](DivisorSigma.md)`[0, n] ==` `Length[`[`Divisors`](Divisors.md)`[n]]`,
and `k = 1` is the sum of divisors used to test perfect numbers.

For a non-negative integer `k`, `DivisorSigma[k, list]` over a packed or `NDArray` `int64`
buffer factors each element on the buffer (the head lowers inside `Compile[]` at a rank-1
array shape too); a negative `k`, such as `sigma_-1`, yields rationals and takes the general
`Expr` path instead. A non-real Gaussian-integer `n`, or `GaussianIntegers -> True`, sums
over the first-quadrant Gaussian divisors. The sign of `n` is ignored, and `n = 0` is left
unevaluated.
