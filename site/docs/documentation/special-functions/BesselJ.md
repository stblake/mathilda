# BesselJ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselJ[n, z]`**

gives the Bessel function of the first kind J\_n(z), a solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 regular at the origin.

**`D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

J\_0(0) = 1, J\_n(0) = 0 for integer n != 0. Has a branch cut along the negative real z axis for non-integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])
```

### Applications (6)

```mathematica
In[3]:= BesselJ[0, 0]
Out[3]= 1

In[4]:= BesselJ[1, 0]
Out[4]= 0

In[5]:= BesselJ[1/2, z]
Out[5]= Sin[z] Sqrt[2/(Pi z)]

In[6]:= Series[BesselJ[0, x], {x, 0, 6}]
Out[6]= 1 - 1/4 x^2 + 1/64 x^4 - 1/2304 x^6 + O[x]^7

In[7]:= N[BesselJ[0, 1], 40]
Out[7]= 0.7651976865579665514497175261026632209093

In[8]:= N[BesselJ[0, 10 + 5 I], 30]
Out[8]= -17.78959112945037151834426180967 + 0.2007116167212048509818027697064*I
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| BesselJ[0, .] over 10^6 | 3.46e+03 s | 1.74e+03 s | 53 s |
| Zeta over 10^6 | 8.17 s | 4.31e+03 s | 5.7 s |
| AiryAi over 10^6 | 3.14 s | 107 s | 58.9 s |
| PolyGamma[0, .] over 10^6 | 1.55 s | 151 s | 8.19 s |
| Gamma over 10^6 | 1.16 s | 1.35 s | 7.34 s |
| Erf over 10^6 | 0.968 s | 1.23 s | 7.43 s |

## Implementation notes

**Algorithm.** `builtin_besselj` (via `besselj_two_arg`) handles `BesselJ[n, z]`, the solution regular at the origin. Exact `z == 0` with classified order gives `J_0(0) = 1`, `J_n(0) = 0` (`Re nu > 0` or `n` a negative integer), `ComplexInfinity` (non-integer `Re nu < 0`), or `Indeterminate` (`Re nu = 0`, `nu != 0`). For a numeric call with an inexact argument: integer order and real `z` take the MPFR-native **`mpfr_jn`** fast path (correctly rounded); any other numeric order/argument goes to the unified complex-MPFR core `bj_core`, which picks the **power series** DLMF 10.2.2 for small/moderate `|z|` (one-Gamma recurrence `t_0 = (z/2)^nu/Gamma(nu+1)`, `t_k = t_{k-1}·(-(z/2)^2)/(k(nu+k))`, with `~2|z|/ln2` guard bits to absorb the `~e^{|z|}` partial-sum cancellation) or the **asymptotic series** DLMF 10.17.3 for large `|z|` (`J_nu(z) ~ sqrt(2/(πz))[cos(w) A - sin(w) B]`, summed to optimal truncation). Integer order is reduced to non-negative order via `J_{-n}(z) = (-1)^n J_n(z)` (the power series would otherwise hit Gamma poles). Half-integer→elementary rewrites live in `src/internal/bessel.m`; `Series` at 0/Infinity in `calculus/series.c`; `D[BesselJ[n,z],z] = (BesselJ[n-1,z] - BesselJ[n+1,z])/2` in `calculus/deriv.c`. Everything else stays symbolic, letting those DownValues and Series/D fire.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`) at explicit working precision. The ND kernel is a binary `REG_B` registration (`NDKB_BesselJ`): over real arrays with **integer** order it calls libc `jn` (declines non-integer order or any complex operand to the `List` path). `Compile[]` lowers `BesselJ[n, z]` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes via that kernel.

**Complexity / limits.** `mpfr_jn` is `O(1)` at machine precision for integer order; the core's series/asymptotic term counts and guard bits scale with `|z|` and precision. Non-integer order has a branch cut along the negative real `z` axis (from the `(z/2)^nu` factor). Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [BesselY](../../special-functions/BesselY/), [BesselI](../../special-functions/BesselI/), [BesselK](../../special-functions/BesselK/), [N](../../arithmetic/N/), [Gamma](../../special-functions/Gamma/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

- DLMF §10.2.2 — power series J_nu(z) = Sum_k (-1)^k (z/2)^{nu+2k}/(k! Gamma(nu+k+1)).
- DLMF §10.17.3 — the large-argument asymptotic expansion of J_nu(z).
- Source: [`src/special_functions/bessel.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/bessel.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m61_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m61_stress.c)

## Notes & additional examples

### Notes

`BesselJ[n, z]` is the Bessel function of the first kind, regular at the origin, with `J_0(0) = 1` and `J_n(0) = 0` for integer `n != 0`. Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2`. There is a branch cut along the negative real axis for non-integer order. Listable.
