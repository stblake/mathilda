---
references:
  - "T. M. Apostol, *Introduction to Analytic Number Theory*, Springer, 1976 — the divisor function sigma_k and its multiplicative formula (Chapter 2)."
---

### Worked examples

```mathematica
In[1]:= DivisorSigma[1, 28]  (* The sum of divisors; 28 is perfect, so sigma_1 = 2 * 28 *)
```

```mathematica
In[1]:= DivisorSigma[0, 60]  (* With k = 0, sigma_0 is the number of divisors *)
```

```mathematica
In[1]:= DivisorSigma[2, 10]  (* Sum of squared divisors: 1 + 4 + 25 + 100 *)
```

```mathematica
In[1]:= DivisorSigma[k, 12]  (* Symbolic k keeps the multiplicative product in closed form *)
```

```mathematica
In[1]:= DivisorSigma[1, Range[6]]  (* Threaded element-wise over a range by Listable *)
```

```mathematica
In[1]:= DivisorSigma[1, 5, GaussianIntegers -> True]  (* Summed over the Gaussian-integer divisors *)
```

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
