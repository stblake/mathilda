# QPochhammer

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`QPochhammer[a, q, n]`**

gives the q-Pochhammer symbol prod\_{k=0}^{n-1} (1 - a q^k).

**`QPochhammer[a, q] gives the infinite q-Pochhammer (a;q)_Inf for |q|<1. The finite form is exact/symbolic for a non-negative integer n; the infinite form evaluates for machine-real a, q. Listable, NumericFunction.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= QPochhammer[a, q, 3]
Out[1]= (1 - a) (1 - a q) (1 - a q^2)
```

### Applications (4)

The finite q-shifted factorial

```mathematica
In[2]:= QPochhammer[a, q, 3]
Out[2]= (1 - a) (1 - a q) (1 - a q^2)
```

```mathematica
In[3]:= QPochhammer[a, q, 0]
Out[3]= 1
```

Exact rational arguments collapse to a rational

```mathematica
In[4]:= QPochhammer[2, 1/3, 4]
Out[4]= -175/729
```

The infinite form at machine precision

```mathematica
In[5]:= QPochhammer[0.5, 0.5]
Out[5]= 0.288788
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

**Algorithm.** `builtin_qpochhammer` evaluates the q-Pochhammer symbol. The
three-argument finite form `QPochhammer[a, q, n] = prod_{k=0}^{n-1} (1 - a q^k)`:
for a non-negative integer `n` the product is built as an `Expr` and handed to
the evaluator, which reduces it exactly for exact `a, q` and at machine / MPFR
precision for inexact `a, q`; `n = 0 -> 1`; a symbolic or non-integer `n` is
left unevaluated (which is exactly what `Product` relies on to emit the closed
form). The two-argument infinite form `QPochhammer[a, q] = (a; q)_inf` is
evaluated only for machine-real `a, q` with `|q| < 1` **and** at least one
inexact operand, by accumulating `(1 - a q^k)` in `double` until the factor
falls below machine epsilon; an all-exact or `|q| >= 1` input stays symbolic.

**Data structures.** `Expr`; the finite form goes through the ordinary evaluator
(GMP for exact, libm/MPFR for inexact); the infinite form is a `double`
accumulation only. ND: binary kernel
`NDK_BIN2(QPochhammer, sf_machine_qpochhammer)`, registered `REG_B`, so
`packed_aware`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The finite form is `O(n)` factors; the infinite form is
**machine `double` only** (no MPFR kernel — so `N[QPochhammer[a, q], p]` returns
a machine-precision value regardless of `p`) and requires `|q| < 1`, summing up
to 100000 factors. `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Product](../../calculus/Product/)

- DLMF §17.2 — q-series and the q-Pochhammer symbol (q-shifted factorial).
- Source: [`src/special_functions/qpochhammer.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/qpochhammer.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_product.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product.c)
- Tests: [`tests/test_product_special.c`](https://github.com/stblake/mathilda/blob/main/tests/test_product_special.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)

## Notes & additional examples

### Notes

`QPochhammer[a, q, n]` is the finite q-Pochhammer symbol (q-shifted factorial)
`prod_{k=0}^{n-1} (1 - a q^k)`, with `QPochhammer[a, q, 0] = 1`. For a
non-negative integer `n` the product is formed and evaluated — exactly for exact
`a, q`, numerically when they are inexact. A symbolic or non-integer `n` is left
unevaluated, which is what `Product` relies on to return the symbol as a closed
form.

The two-argument `QPochhammer[a, q]` is the infinite product `(a; q)_inf`,
evaluated for machine-real `a, q` with `|q| < 1`. That infinite form is a
machine-precision kernel only, so a request for extra digits via `N[..., p]`
still returns a machine-precision value.
