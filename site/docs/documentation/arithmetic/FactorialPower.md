# FactorialPower

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FactorialPower[n, k]`**

The falling factorial n (n - 1) (n - 2) ... (n - k + 1). For non-negative integer k, expands to a product of k linear factors. Equivalent to n! / (n - k)! when both n and k are non-negative integers. For a non-integer numeric k (under N), evaluates as Gamma\[n+1\]/Gamma\[n-k+1\].

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FactorialPower[5, 3]
Out[1]= 60

In[2]:= FactorialPower[n, 3]
Out[2]= n (-2 + n) (-1 + n)

In[3]:= FactorialPower[n, 0]
Out[3]= 1

In[4]:= D[x^n, {x, k}]
Out[4]= FactorialPower[n, k] x^(-k + n)
```

### Applications (3)

```mathematica
In[5]:= FactorialPower[10, 3]
Out[5]= 720

In[6]:= FactorialPower[x, 4]
Out[6]= x (-3 + x) (-2 + x) (-1 + x)

In[7]:= FactorialPower[10, 4] == 10!/6!
Out[7]= True
```

## Implementation notes

`builtin_factorialpower` computes the falling factorial `n^(k) = n(n-1)...(n-k+1)`. It requires `k` to be a concrete non-negative integer (`k == 0` gives `1`; negative or symbolic `k` returns `NULL`). For integer/bignum `n` it accumulates the exact product in GMP (`mpz_mul` of `n-i`), normalising via `expr_bigint_normalize`. For symbolic `n` with small `k` (`<= 32`) it expands the literal product `Times[n, n-1, ..., n-k+1]` through `eval_and_free` so `Expand`/`D` can act on it. The 1-argument step-form is not handled here.

- `Protected`, `Listable`, `NumericFunction`.
- For non-negative integer $n$ and non-negative integer $k$: exact GMP product.
- For symbolic $n$ with concrete $k \le 32$: expands to an explicit product
  `Times[n, n - 1, ..., n - k + 1]` so `Expand` and `D` can act on it.
- Equivalent to $n! / (n - k)!$ when both arguments are non-negative integers.
- For a **non-integer numeric** $k$ (under `N`), evaluates the generalised form
  $\Gamma(n+1)/\Gamma(n-k+1)$ (real, complex, arbitrary precision); exact
  non-integer $k$ stays symbolic.
- Used as the symbolic-order derivative of a power: $D[x^n, \{x, k\}] = \mathrm{FactorialPower}[n, k]\, x^{n-k}$.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Expand](../../algebra/Expand/), [D](../../calculus/D/), [N](../../arithmetic/N/)

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_deriv_symbolic_order.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv_symbolic_order.c)
- Tests: [`tests/test_numeric_domain.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_domain.c)

## Notes & additional examples

### Notes

`FactorialPower[n, k]` is the falling factorial `n (n-1) (n-2) ... (n-k+1)`, the product of `k` descending linear factors. For symbolic `n` it expands to that product, e.g. `FactorialPower[x, 4] = x (x-1) (x-2) (x-3)`, the basis polynomial that makes the forward difference operator behave like differentiation (Δ x^(k) = k x^(k-1)). For non-negative integers it equals the ratio of factorials `n! / (n-k)!`, as the identity check confirms.
