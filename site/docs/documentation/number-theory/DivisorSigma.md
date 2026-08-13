# DivisorSigma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DivisorSigma[k, n] gives the divisor function sigma_k(n), the sum of the k-th powers of the divisors of n. DivisorSigma[k, n, GaussianIntegers -> True] sums over Gaussian-integer divisors.`**

## Examples (8)

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

## Options & behaviour

> **Packed arrays.** `DivisorSigma[k, list]` over an `int64` buffer factors
> each element by trial division in `int64`, with no GMP allocation per
> element. A non-negative integer `k` only: `DivisorSigma[-1, n]` is a
> `Rational`, which no buffer holds.

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)
