# Hyperfactorial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hyperfactorial[n]`**

gives the hyperfactorial prod\_{k=1}^{n} k^k.

<details>
<summary>Notes</summary>

Exact (GMP) for a non-negative integer n. A non-integer numeric order (under N) evaluates via Gamma\[n+1\]^n / BarnesG\[n+1\] (real, complex, arbitrary precision); symbolic orders stay unevaluated. Listable, NumericFunction.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Hyperfactorial[4]
Out[1]= 27648

In[2]:= N[Hyperfactorial[5.0]]
Out[2]= 8.64e+07

In[3]:= N[Hyperfactorial[7/2], 30]
Out[3]= 1282.1220994534574594154227123
```

### Applications (5)

The product 1^1 2^2 3^3 4^4

```mathematica
In[4]:= Hyperfactorial[4]
Out[4]= 27648
```

```mathematica
In[5]:= Hyperfactorial[0]
Out[5]= 1

In[6]:= Table[Hyperfactorial[n], {n, 0, 5}]
Out[6]= {1, 1, 4, 108, 27648, 86400000}
```

Exact non-integer orders stay symbolic

```mathematica
In[7]:= Hyperfactorial[1/2]
Out[7]= Hyperfactorial[1/2]
```

The Barnes-G continuation, to 20 digits

```mathematica
In[8]:= N[Hyperfactorial[5/2], 20]
Out[8]= 15.9842119922237168526
```

## Algorithm

Mathilda -- Hyperfactorial.

```text
  Hyperfactorial[n] = prod_{k=1}^{n} k^k   (H(0) = H(1) = 1).
```

Exact for a non-negative integer order (GMP); non-positive-integer, non-integer, or symbolic orders are left unevaluated (the analytic

```text
K-function continuation is not implemented).  N at an integer order routes
```

through the exact value and numericalize, so machine and MPFR precision come

```text
for free.  Used by Product to recognise prod k^k = Hyperfactorial[n].
```

Memory: honours the builtin ownership contract (never frees res; returns a fresh Expr* or NULL; clears every GMP temporary).

## Implementation notes

**Algorithm.** `builtin_hyperfactorial` computes `H(n) = prod_{k=1}^{n} k^k`
(with `H(0) = H(1) = 1`). A non-negative exact integer order is done exactly in
GMP — a running `mpz_t` accumulating `mpz_ui_pow_ui(k, k)` — capped at
`HYPERFACTORIAL_MAX_N = 20000` to protect memory. Negative-integer, exact
non-integer, and symbolic orders are left unevaluated (`NULL`). For an inexact
argument (under `N[...]`) it falls back to the Barnes-G analytic continuation
`H(z) = Gamma[z+1]^z / BarnesG[z+1]`, built as an `Expr` tree and evaluated so it
inherits `Gamma`/`BarnesG`'s numeric kernels (machine and MPFR); a result that
overflows a `double` is promoted to an extended-exponent MPFR real or made
`Overflow[]` by `numeric_promote_result_if_degenerate`. `Product` relies on this
head to recognise `prod k^k`.

**Data structures.** `Expr`; GMP `mpz_t` for the exact integer product; the
continuation is an `Expr` over `Gamma`/`BarnesG` evaluated through the ordinary
numeric stack (libm / MPFR). No ND kernel and not on `pack.c`'s `AWARE` list;
`Compile[]` does not lower it (`CompileDiagnostics` reports `Compiled -> False`).
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The exact path is `O(n)` growing-bignum multiplies,
hard-capped at `n = 20000`; the continuation costs one `BarnesG`/`Gamma`
evaluation. Exact non-integer and negative-integer orders stay symbolic (the
`K`-function continuation is reached only numerically, via Barnes G). No
closed-form derivative rule.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [Product](../../calculus/Product/)

- DLMF §5.17 — Barnes G-function (the analytic continuation Hyperfactorial[z] = Gamma[z+1]^z / BarnesG[z+1]).
- Source: [`src/special_functions/hyperfactorial.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/hyperfactorial.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_numeric_domain.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_domain.c)
- Tests: [`tests/test_product.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product.c)
- Tests: [`tests/test_product_special.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product_special.c)

## Notes & additional examples

### Notes

`Hyperfactorial[n]` is the product `prod_{k=1}^{n} k^k`, with
`Hyperfactorial[0] = Hyperfactorial[1] = 1`. A non-negative integer order is
computed exactly in arbitrary-precision integer arithmetic (capped at
`n = 20000`), so the values grow very fast: `1, 1, 4, 108, 27648, 86400000, ...`.

A non-integer or complex order is evaluated only numerically (under `N`),
through the Barnes-G analytic continuation `Hyperfactorial[z] = Gamma[z+1]^z /
BarnesG[z+1]`, which reproduces the integer values and extends the function to
the whole plane. Exact non-integer and negative-integer orders are left
symbolic.

The head exists partly so that `Product` can recognise `prod k^k` and return it
in closed form.
