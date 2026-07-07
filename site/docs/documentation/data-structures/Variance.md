# Variance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Variance[data] gives the unbiased variance estimate of the elements in data.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]
Out[1]= 3

In[2]:= Variance[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[2]= 4

In[3]:= StandardDeviation[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[3]= 2
```

## Implementation notes

**Algorithm.** `builtin_variance` computes the *sample* variance (divisor `n-1`, requiring `n > 1`). It first reduces matrices column-wise via `apply_columnwise`, then requires a `List`. For real-valued data it uses Welford's online algorithm (running mean `m` and sum-of-squares `s`) and returns `expr_new_real(s/(n-1))`. For exact integer/rational data it does the computation in `int64_t` numerator/denominator pairs: it first accumulates the sum (hence the mean), then accumulates `Sum[(x_i - mean)^2]` as exact rationals (reducing by `gcd`), and returns `make_rational(sq_sum_n, sq_sum_d * (n-1))`. The symbolic fallback evaluates `Mean[data]`, forms `Sum[(x - mu) Conjugate[x - mu]]` (so complex/symbolic data gives the Hermitian variance), and divides by `n-1`.

**Limits.** The exact path uses fixed-width `int64_t`, so large rationals can overflow; `n <= 1` returns `NULL`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Variance[{1, 2, 3, 4, 5}]
Out[1]= 5/2
```

`Variance` is the *unbiased* (sample) estimator, dividing the sum of squared deviations by `n - 1`. On exact rational data the answer stays exact:

```mathematica
In[1]:= Variance[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[1]= 32/7
```

Feeding the same integers as 40-digit reals propagates the precision all the way through the deviation sum:

```mathematica
In[1]:= Variance[N[{1, 1, 2, 3, 5, 8, 13}, 40]]
Out[1]= 19.571428571428571428571428571428571428568
```

### Notes

`Variance[data]` gives the unbiased variance estimate (Bessel-corrected, `1/(n-1)` normalization) of the elements in `data`. Exact inputs yield exact rational results; arbitrary-precision inputs carry their precision through the computation.
