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

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Hyperfactorial[4]
Out[1]= 27648

In[2]:= N[Hyperfactorial[5.0]]
Out[2]= 8.64e+07

In[3]:= N[Hyperfactorial[7/2], 30]
Out[3]= 1282.122099453457459415422713168
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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[N](../../arithmetic/N/), [Product](../../calculus/Product/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_numeric_domain.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_domain.c)
- Tests: [`tests/test_product.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product.c)
- Tests: [`tests/test_product_special.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product_special.c)
