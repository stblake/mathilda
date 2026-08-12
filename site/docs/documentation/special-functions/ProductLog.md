# ProductLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ProductLog[z]`**

gives the principal solution w of z == w e^w (the Lambert W function).

**`ProductLog[k, z] gives the k-th solution (k any integer, k == 0 the`**

**`ProductLog[E] = 1, ProductLog[-Pi/2] = I Pi/2 and ProductLog[k, 0] =`**

**`D[ProductLog[z], z] = ProductLog[z]/(z (1 + ProductLog[z])). Listable.`**

<details>
<summary>Notes</summary>

principal branch); branches are ordered by imaginary part. ProductLog\[z\] is real for z \>= -1/e and has a branch cut along (-Infinity, -1/e\]. Exact values include ProductLog\[0\] = 0, ProductLog\[-1/E\] = -1, -Infinity for k != 0. Inexact real or complex arguments evaluate numerically at machine or arbitrary (MPFR) precision. Satisfies

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1
```

## Algorithm

Mathilda -- ProductLog, the Lambert W function.

```text
  ProductLog[z]     principal branch W_0(z): the solution w of z = w e^w.
  ProductLog[k, z]  the k-th branch W_k(z), k any integer (k == 0 principal).
```

Evaluation is layered so each kind of argument takes the cheapest accurate route:

```text
  exact special values   ->  ProductLog[0] = 0, ProductLog[E] = 1,
                             ProductLog[-1/E] = -1, ProductLog[-Pi/2] = I Pi/2,
                             ProductLog[+-Infinity/ComplexInfinity] = Infinity,
                             ProductLog[k, 0] = -Infinity  (k != 0)
  numeric (real/complex)  ->  unified complex-MPFR Halley core; the result is
                             a Real / MPFR leaf when it is real-valued for the
                             chosen branch, otherwise Complex[..]
  everything else        ->  stays symbolic (return NULL)
```

The numeric core (pl_core) builds on the shared `ncpx` complex-MPFR toolkit (numeric_complex.h). It seeds an initial approximation by region --

```text
  - branch-point series in p = sqrt(2(e z + 1)) near z = -1/e (branches 0,-1);
  - the Maclaurin seed z(1 - z + 3/2 z^2) for the principal branch near 0;
  - otherwise the asymptotic L1 - L2 + L2/L1 with L1 = log z + 2 pi i k,
    L2 = log L1
```

-- and refines it with Halley's cubically-convergent iteration

```text
  w <- w - (w e^w - z) / (e^w (w+1) - (w+2)(w e^w - z)/(2w+2))   (Corless 1996).
```

A real seed keeps the whole iteration exactly real (every ncpx op preserves a zero imaginary part), so real-valued branches return a real leaf with no imaginary noise. Working precision carries guard bits above the requested output precision.

```text
D[ProductLog[z], z] = ProductLog[z] / (z (1 + ProductLog[z]))  (calculus/deriv.c).
```

Series at 0, at the branch point -1/E, and at Infinity live in calculus/series.c.

Attributes: Listable, NumericFunction, Protected, ReadProtected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## See also

[N](../../arithmetic/N/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
