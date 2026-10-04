# CoshIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CoshIntegral[z]`**

gives the hyperbolic cosine integral Chi(z) = EulerGamma + Log\[z\] + Integral\_0^z (Cosh\[t\] - 1)/t dt.

**`CoshIntegral[0] = -Infinity, CoshIntegral[Infinity] = Infinity,`**

**`CoshIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0\]. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[CoshIntegral\[z\], z\] = Cosh\[z\]/z. Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

Logarithmic singularity at the origin

```mathematica
In[1]:= CoshIntegral[0]
Out[1]= -Infinity
```

Grows without bound on the positive real axis

```mathematica
In[2]:= CoshIntegral[Infinity]
Out[2]= Infinity
```

Arbitrary-precision value through the MPFR kernel

```mathematica
In[3]:= N[CoshIntegral[1], 20]
Out[3]= 0.837866940980208240895
```

The derivative is Cosh[x]/x

```mathematica
In[4]:= D[CoshIntegral[x], x]
Out[4]= Cosh[x]/x
```

Threads element-wise over the packed real list

```mathematica
In[5]:= CoshIntegral[{1.0, 2.0, 3.0}]
Out[5]= {0.837867, 2.45267, 4.96039}
```

## Algorithm

Mathilda -- the hyperbolic cosine integral

```text
  Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.

  CoshIntegral[z]   Chi(z)
```

Chi has a logarithmic singularity at z = 0 and a branch cut running along the

```text
negative real axis (-Infinity, 0].  It is the imaginary-axis sibling of the
cosine integral: Chi(z) = Ci(i z) - i Pi/2.  Evaluation is layered so each
```

argument takes the cheapest route:

```text
  exact special values     ->  -Infinity, Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Chi(z) = EulerGamma + Log(z) + Sum_{k>=1} z^(2k) / (2k (2k)!),
```

valid on the principal branch for all z != 0 -- the principal Log(z) supplies the correct +-i Pi jump across the cut, so no folding is needed for the

```text
convergent path.  Every series term is positive, so on the real axis there is
NO catastrophic cancellation and only a small fixed guard is needed.  For large
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Chi(z) ~ sinh(z) F(z) + cosh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

(the same F, G SinhIntegral computes; only the combination differs -- Shi uses

```text
cosh F + sinh G).  The bare part sinh F + cosh G is even and real on the
positive real axis (no constant).  Off the real axis, matching Mathematica's
```

Series[CoshIntegral[z],{z,Infinity,k}], a piecewise Stokes term restores the principal branch:

```text
  Chi(z) = B2(z) + i K,   B2(z) = sinh(z) F(z) + cosh(z) G(z),
    K = -Pi/2                       for Im(z) < 0,
    K = Pi sgn+(Re z) - Pi/2        for Im(z) > 0   (sgn+(0) = +1, from above).
```

(On the imaginary axis Im(z) > 0, Re(z) = 0 gives K = +Pi/2, so

```text
Chi(+- i Infinity) = +- i Pi/2.)  A negative real x is handled by the real
```

path as the from-above branch value Complex[Chi(|x|), Pi], matching the log jump; and Chi(-Infinity) -> Infinity (the real part dominates).

Machine-real results that overflow a C double (e.g. CoshIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_coshintegral` handles `CoshIntegral[z] = Chi(z)`, the imaginary-axis sibling of `Ci` (`Chi(z) = Ci(iz) - iπ/2`), with a logarithmic singularity at `0` and a branch cut along `(-Infinity, 0]`. Exact special values first: `Chi[0] = -Infinity`, `Chi[±Infinity] = Infinity`, `Chi[±I Infinity] = ±I π/2`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real or complex** argument routes to the MPFR kernel: for moderate `|z|` the **convergent Maclaurin series** (`Chi(z) = γ + Log(z) + Sum_{k>=1} z^{2k}/(2k(2k)!)` — the trig series with the alternating sign removed, so every term is positive and the real axis needs no cancellation guard, only a fixed 64-bit guard); the principal `Log(z)` supplies the `±iπ` cut jump. For large `|z|` the **asymptotic expansion** `sinh(z) F(z) + cosh(z) G(z)` plus a piecewise Stokes constant `K` (`-π/2` for `Im z < 0`, `π sgn⁺(Re z) - π/2` for `Im z > 0`) restores the principal branch. A negative real `x` gives the from-above value `Complex[Chi(|x|), Pi]`; machine-real results that overflow a `double` (Chi grows like `e^{|x|}`) are emitted as a 53-bit MPFR real. The complex path uses the shared `ncpx` toolkit; a `USE_MPFR=0` build uses a machine-double series/asymptotic.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`). The ND kernel is a real `REG_U` registration (`NDKU_CoshIntegral`): real buffers via `ndk_CoshIntegral_r` → `sf_machine_chi`, complex buffers via `coshintegral_machine_complex` (with the shared `sf_series_usable` cancellation gate). `Compile[]` lowers `CoshIntegral` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|` and precision. Logarithmic singularity at `0`, branch cut on the negative real axis. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [CosIntegral](../../special-functions/CosIntegral/), [Log](../../elementary-functions/Log/)

- DLMF §6.2.16 — the hyperbolic cosine integral Chi(z) = γ + Log(z) + Int_0^z (cosh t - 1)/t dt.
- Source: [`src/special_functions/coshintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/coshintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_coshintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coshintegral.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_series.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series.c)

## Notes & additional examples

### Notes

`CoshIntegral[z] = Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt`
is the hyperbolic sibling of `CosIntegral` (`Chi(z) = Ci(i z) - i Pi/2`). It has
a logarithmic singularity at `0` and a branch cut along the negative real axis;
a negative real argument returns the from-above branch value `Chi(|x|) + i Pi`.

The numeric kernel sums a convergent Maclaurin series for moderate `|z|` and an
asymptotic expansion for large `|z|` (where `Chi` grows like `e^{|z|}`;
machine-real results that overflow a C `double` are kept as extended-exponent
reals). `CoshIntegral` carries a real `NDArray` kernel and lowers under
`Compile[]` at both scalar and rank-1 array shapes.
