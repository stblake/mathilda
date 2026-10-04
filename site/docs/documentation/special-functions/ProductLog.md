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

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1
```

### Applications (7)

```mathematica
In[3]:= ProductLog[0]
Out[3]= 0
```

Since e = 1 * e^1, the principal value is 1

```mathematica
In[4]:= ProductLog[E]
Out[4]= 1
```

The branch point where the k = 0 and k = -1 branches meet

```mathematica
In[5]:= ProductLog[-1/E]
Out[5]= -1
```

```mathematica
In[6]:= ProductLog[-Pi/2]
Out[6]= (1/2*I) Pi
```

The omega constant

```mathematica
In[7]:= N[ProductLog[1], 30]
Out[7]= 0.56714329040978387299996866221
```

The k = -1 branch

```mathematica
In[8]:= ProductLog[-1, -0.2]
Out[8]= -2.54264
```

```mathematica
In[9]:= D[ProductLog[x], x]
Out[9]= ProductLog[x]/(x (1 + ProductLog[x]))
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

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_productlog` evaluates the Lambert W function:
`ProductLog[z]` is the principal branch `W_0`, `ProductLog[k, z]` the `k`-th
branch (`k` an explicit machine integer). Exact special values:
`ProductLog[0] = 0`, `ProductLog[E] = 1`, `ProductLog[-1/E] = -1`,
`ProductLog[-Pi/2] = I Pi/2`, `ProductLog[±Infinity] = Infinity`, and
`ProductLog[k, 0] = -Infinity` for `k != 0`. Numeric uses a unified complex-MPFR
core (`ncpx`): a region-chosen seed — the branch-point series in
`p = sqrt(2(e z + 1))` near `z = -1/e` (branches 0 and -1), the Maclaurin seed
`z(1 - z + 3/2 z^2)` for the principal branch near 0, else the asymptotic
`L1 - L2 + L2/L1` with `L1 = log z + 2 pi i k`, `L2 = log L1` — refined by
Halley's cubically-convergent iteration (Corless 1996). A real seed keeps the
iteration exactly real, so real-valued branches return a real leaf.
`D[ProductLog[z], z] = ProductLog[z]/(z(1 + ProductLog[z]))` lives in
`calculus/deriv.c`; Series at 0, at `-1/E`, and at Infinity in
`calculus/series.c`.

**Data structures.** `Expr`; the shared `ncpx` (`mpfr_t` re/im) complex toolkit
for the Halley core. ND: real-only unary kernel `NDKU_ProductLog = { NULL,
ndk_ProductLog_r, ... }` (the principal branch, via `sf_machine_productlog`),
registered `REG_U`, so `packed_aware`. Attributes: `Listable`, `NumericFunction`,
`Protected`.

**Complexity / limits.** Halley converges cubically (bounded at 100 iterations
plus one polishing step). A real leaf is returned only on the real-valued
domains (`W_0` for `x >= -1/e`, `W_{-1}` for `-1/e <= x < 0`); otherwise a
`Complex[..]`. Only numeric inputs evaluate — exact non-special arguments stay
symbolic. `Compile[]` lowers at both scalar and rank-1 array shapes
(`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- DLMF §4.13 — the Lambert W-function.
- R. M. Corless et al., On the Lambert W function, Adv. Comput. Math. 5 (1996) 329-359.
- Source: [`src/special_functions/productlog.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/productlog.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`ProductLog[z]` is the Lambert W function: the principal branch `W_0(z)`, the
solution `w` of `z = w e^w`. `ProductLog[k, z]` selects the `k`-th branch (`k`
an explicit integer); `k = 0` is the principal branch and `k = -1` is the other
real branch.

A handful of exact values are returned directly (`0`, `E`, `-1/E`, `-Pi/2`);
everything else is evaluated numerically by a Halley iteration in
arbitrary-precision complex arithmetic, returning a real result on the
real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}` for `-1/e <= x < 0`) and a
`Complex[..]` otherwise. `N[ProductLog[1], 30]` is the omega constant `Omega`,
the solution of `Omega e^Omega = 1`. The derivative is
`ProductLog[z]/(z (1 + ProductLog[z]))`.
