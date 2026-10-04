# BarnesG

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BarnesG[z]`**

gives the Barnes G-function.

<details>
<summary>Notes</summary>

G(z+1) = Gamma\[z\] G(z) with G(1)=G(2)=1; for a positive integer n, G(n+1) = prod\_{k=1}^{n-1} k! (exact via GMP), and G(m)=0 for non-positive integer m. A non-integer numeric order (under N) evaluates from the Barnes asymptotic expansion plus the Gamma recurrence (real, complex, arbitrary precision); symbolic orders stay unevaluated. Listable, NumericFunction.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= BarnesG[5]
Out[1]= 12

In[2]:= N[BarnesG[6.0]]
Out[2]= 288.0

In[3]:= N[BarnesG[13/2], 30]
Out[3]= 2548.745769568498989735906106363

In[4]:= N[BarnesG[2.5 + 1.0 I]]
Out[4]= 0.743798 - 0.0953168*I
```

### Worked examples (1)

```mathematica
In[5]:= Product[Gamma[i], {i, 1, n-1}]
Out[5]= BarnesG[n]
```

### Applications (5)

The superfactorial 1! 2! 3! = 12

```mathematica
In[6]:= BarnesG[5]
Out[6]= 12
```

Threads over the list: G(n) at the first six integers

```mathematica
In[7]:= BarnesG[Range[6]]
Out[7]= {1, 1, 1, 2, 12, 288}
```

The functional equation gives G(n+1)/G(n) = Gamma[n] = 7!

```mathematica
In[8]:= BarnesG[9]/BarnesG[8]
Out[8]= 5040
```

Evaluated at a half-integer through the asymptotic continuation

```mathematica
In[9]:= N[BarnesG[7/2], 20]
Out[9]= 1.25964825749519214408
```

Non-positive integers are (double) zeros

```mathematica
In[10]:= BarnesG[-2]
Out[10]= 0
```

## Algorithm

Mathilda -- BarnesG[z], the Barnes G-function.

```text
  G(1) = G(2) = 1,   G(z+1) = Gamma[z] G(z),
  integer:  G(n+1) = prod_{k=1}^{n-1} k!   (the superfactorial),
            G(m) = 0 for non-positive integer m (double zeros).
```

Exact for integer orders (GMP); non-integer orders are left unevaluated (the

```text
LogGamma/zeta'(-1) asymptotic continuation is not implemented).  N at an
integer order routes through the exact value and numericalize.  Used by
```

Product to recognise prod_{k=1}^{n-1} Gamma[k] = BarnesG[n].

Memory: honours the builtin ownership contract.

## Implementation notes

**Algorithm.** `builtin_barnesg` uses `G(1) = G(2) = 1`, `G(z+1) = Gamma(z) G(z)`. For an exact integer `z`: non-positive integers are zeros (`G(0) = G(-1) = … = 0`), and for `z >= 3` it builds the exact **superfactorial** `G(z) = prod_{k=1}^{z-2} k!` with GMP, guarded by a runaway cap (`BARNESG_MAX_N = 2000`). Barnes G has no elementary closed form off the integers, so a non-integer *inexact* argument (under `N`) routes to `barnesg_numeric`, which combines the recurrence `G(w) = exp(A(Z)) / prod_{j=0}^{m-1} Gamma(w+j)`, `Z = w + m - 1`, with the Barnes asymptotic expansion of `log G(Z+1)` (carrying twelve fixed Bernoulli-number terms `c_k = B_{2k+2}/((2k)(2k+2))` plus the Glaisher–Kinkelin constant). The whole continuation is assembled as a Mathilda expression and handed to `N[]`, so it inherits `Gamma`/`Log`/`Exp`'s complex and MPFR kernels and `Glaisher` at full precision; the shift `m` is chosen so the twelve terms reach the requested digits (~40 digits accuracy). Degenerate results that overflow a `double` (e.g. `BarnesG[100.]` ~ `3·10^6626`) are promoted to an extended-exponent MPFR real via `numeric_promote_result_if_degenerate`. Everything else stays symbolic.

**Data structures.** `Expr`; GMP `mpz_t` for the superfactorial product; the numeric continuation is an assembled `Expr` tree evaluated through `N`, so its backend is whatever `Gamma`/`Log`/`Exp` use (libm at machine precision, MPFR at arbitrary precision). No ND kernel and no `Compile[]` lowering (`CompileDiagnostics` reports `Compiled -> False`): `BarnesG` is not an element-wise machine primitive — it is exact-integer GMP work or an assembled high-precision continuation.

**Complexity / limits.** Exact only on integers up to the cap of 2000 (the superfactorial is `O(z)` big-integer multiplications). Non-integer values exist only under `N` (asymptotic continuation, twelve fixed Bernoulli terms, ~40 digits). Used by `Product` to recognise `prod_{k=1}^{n-1} Gamma[k] = BarnesG[n]`. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Product](../../calculus/Product/), [N](../../arithmetic/N/), [Gamma](../../special-functions/Gamma/), [Log](../../elementary-functions/Log/), [Exp](../../elementary-functions/Exp/)

- DLMF §5.17 — Barnes G-function: G(z+1) = Gamma(z) G(z), the Glaisher–Kinkelin constant, and the asymptotic expansion of log G.
- Source: [`src/special_functions/barnesg.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/barnesg.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndsolve_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_compile.c)

## Notes & additional examples

### Notes

The Barnes G-function satisfies `G(1) = G(2) = 1` and the functional equation
`G(z+1) = Gamma[z] G(z)`. On the positive integers it is the **superfactorial**
`G(n) = prod_{k=1}^{n-2} k!`, computed exactly with GMP; the non-positive
integers are its zeros.

Off the integers there is no elementary closed form, so a non-integer argument
has a value only under `N`, where an asymptotic continuation (twelve fixed
Bernoulli terms plus the Glaisher–Kinkelin constant) carries it to roughly forty
digits. `BarnesG` has no `NDArray` kernel and does not lower under `Compile[]`:
it is exact-integer GMP work or an assembled high-precision continuation, not an
element-wise machine primitive.

`Product` uses `BarnesG` to recognise `prod_{k=1}^{n-1} Gamma[k] = BarnesG[n]`.
