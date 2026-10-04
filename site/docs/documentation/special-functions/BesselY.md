# BesselY

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselY[n, z]`**

gives the Bessel function of the second kind Y\_n(z), the solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.

<details>
<summary>Notes</summary>

Y\_0(0) = -Infinity, Y\_n(0) = ComplexInfinity for integer n != 0; Y\_n has a logarithmic branch point at 0 and a branch cut along the negative real z axis, with Y\_{-n} = (-1)^n Y\_n for integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision; D\[BesselY\[n, z\], z\] = (BesselY\[n-1, z\] - BesselY\[n+1, z\])/2. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])
```

### Applications (5)

```mathematica
In[3]:= N[BesselY[1, 3.0]]
Out[3]= 0.324674

In[4]:= BesselY[1/2, z]
Out[4]= -Cos[z] Sqrt[2/(Pi z)]

In[5]:= N[BesselY[0, 1], 40]
Out[5]= 0.088256964215676957982926766023515162827815

In[6]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5], 30]
Out[6]= 0.127323954473516268615107010698

In[7]:= N[2/(5 Pi), 30]
Out[7]= 0.127323954473516268615107010698
```

## Implementation notes

**Algorithm.** `builtin_bessely` (via `bessely_two_arg`) handles `BesselY[n, z]`, the second-kind solution, singular at the origin. Exact `z == 0` with classified order gives `Y_0(0) = -Infinity`, `0` (a negative half-odd-integer order, where `cos(nu π) = 0` cancels the divergent term so `Y_{-1/2} = J_{1/2}`), `Indeterminate` (`Re nu = 0`, `nu != 0`), or `ComplexInfinity` otherwise. For a numeric call: integer order and real `z > 0` take the MPFR-native **`mpfr_yn`** fast path (correctly rounded); otherwise the unified core `by_core` routes among (i) the **connection** `Y_nu(z) = (J_nu(z) cos(nu π) - J_{-nu}(z)) / sin(nu π)` for small `|z|`, non-integer order, (ii) the **logarithmic series** DLMF 10.8.1 for small `|z|`, integer order, and (iii) the **asymptotic series** DLMF 10.17.4 (`Y_nu(z) ~ sqrt(2/(πz))[sin(w) A + cos(w) B]`) for large `|z|`, summed to optimal truncation. Half-integer→elementary rewrites live in `src/internal/bessel.m`; `Series`/`D` in `calculus/`. Everything else stays symbolic.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` at explicit working precision. The ND kernel is a binary `REG_B` registration (`NDKB_BesselY`): over real arrays with **integer** order it calls libc `yn` (declines non-integer order or any complex operand to the `List` path). `Compile[]` lowers `BesselY[n, z]` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes via that kernel.

**Complexity / limits.** `mpfr_yn` is `O(1)` at machine precision for integer order; the core's term counts and guard bits scale with `|z|` and precision. Branch cut along the negative real `z` axis. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [BesselK](../../special-functions/BesselK/), [BesselI](../../special-functions/BesselI/)

- DLMF §10.8.1 — the logarithmic series for integer-order Y_n(z) near the origin.
- DLMF §10.17.4 — the large-argument asymptotic expansion of Y_nu(z).
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

`BesselY[n, z]` is the Bessel function of the second kind, singular at the origin: `Y_0(0) = -Infinity`, `Y_n(0) = ComplexInfinity` for integer `n != 0`, with a logarithmic branch point at 0 and a branch cut along the negative real axis (`Y_{-n} = (-1)^n Y_n` for integer `n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`. Listable.
