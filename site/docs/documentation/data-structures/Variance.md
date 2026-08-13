# Variance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Variance[data] gives the unbiased variance estimate of the elements in data.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]
Out[1]= 3

In[2]:= Variance[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[2]= 4

In[3]:= StandardDeviation[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[3]= 2
```

### Applications (3)

```mathematica
In[4]:= Variance[{1, 2, 3, 4, 5}]
Out[4]= 5/2

In[5]:= Variance[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[5]= 32/7

In[6]:= Variance[N[{1, 1, 2, 3, 5, 8, 13}, 40]]
Out[6]= 19.571428571428571428571428571428571428568
```

## Implementation notes

**Algorithm.** `builtin_variance` computes the *sample* variance (divisor `n-1`, requiring `n > 1`). It first reduces matrices column-wise via `apply_columnwise`, then requires a `List`. For real-valued data it uses Welford's online algorithm (running mean `m` and sum-of-squares `s`) and returns `expr_new_real(s/(n-1))`. For exact integer/rational data it does the computation in `int64_t` numerator/denominator pairs: it first accumulates the sum (hence the mean), then accumulates `Sum[(x_i - mean)^2]` as exact rationals (reducing by `gcd`), and returns `make_rational(sq_sum_n, sq_sum_d * (n-1))`. The symbolic fallback evaluates `Mean[data]`, forms `Sum[(x - mu) Conjugate[x - mu]]` (so complex/symbolic data gives the Hermitian variance), and divides by `n-1`.

**Limits.** The exact path uses fixed-width `int64_t`, so large rationals can overflow; `n <= 1` returns `NULL`.

**Attributes:** `Protected`.

## References

**See also:** [Median](../../data-structures/Median/), [StandardDeviation](../../data-structures/StandardDeviation/), [Mean](../../data-structures/Mean/)

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`Variance[data]` gives the unbiased variance estimate (Bessel-corrected, `1/(n-1)` normalization) of the elements in `data`. Exact inputs yield exact rational results; arbitrary-precision inputs carry their precision through the computation.
