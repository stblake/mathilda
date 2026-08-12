# QPochhammer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`QPochhammer[a, q, n]`**

gives the q-Pochhammer symbol prod\_{k=0}^{n-1} (1 - a q^k).

**`QPochhammer[a, q] gives the infinite q-Pochhammer (a;q)_Inf for |q|<1. The finite form is exact/symbolic for a non-negative integer n; the infinite form evaluates for machine-real a, q. Listable, NumericFunction.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= QPochhammer[a, q, 3]
Out[1]= (1 - a) (1 - a q) (1 - a q^2)
```

## Algorithm

Mathilda -- QPochhammer, the q-Pochhammer symbol (q-shifted factorial).

```text
  QPochhammer[a, q, n] = prod_{k=0}^{n-1} (1 - a q^k)
  QPochhammer[a, q]     = prod_{k=0}^{Infinity} (1 - a q^k)   ((a;q)_inf)
```

Finite form (3 args): for a non-negative integer n the product is built and handed to the evaluator, which reduces it exactly for exact a, q and at

```text
machine / MPFR precision for inexact a, q (so N works through it).  A
```

symbolic / non-integer n is left unevaluated -- which is exactly what Product relies on to emit QPochhammer[a, q, n] as a closed form.

Infinite form (2 args): evaluated for machine-real a, q with |q| < 1 by

```text
accumulating factors until they fall below machine epsilon.  Symbolic or
|q| >= 1 inputs stay unevaluated.
```

Memory: honours the builtin ownership contract (never frees res).

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[Product](../../calculus/Product/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_product.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product.c)
- Tests: [`tests/test_product_special.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product_special.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)
