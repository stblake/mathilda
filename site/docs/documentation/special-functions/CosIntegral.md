# CosIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CosIntegral[z]`**

gives the cosine integral Ci(z) = -Integral\_z^Infinity Cos\[t\]/t dt.

**`CosIntegral[0] = -Infinity, CosIntegral[Infinity] = 0,`**

**`CosIntegral[-Infinity] = I Pi, CosIntegral[+-I Infinity] = Infinity.`**

<details>
<summary>Notes</summary>

Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0\]. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[CosIntegral\[z\], z\] = Cos\[z\]/z. Listable.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= CosIntegral[2.8]
Out[1]= 0.186488

In[2]:= N[CosIntegral[2], 50]
Out[2]= 0.422980828774864995698565153198255894135737756306189

In[3]:= CosIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[3]= {I Pi, 0, Infinity, Infinity}

In[4]:= CosIntegral[-2.]
Out[4]= 0.422981 + 3.14159*I

In[5]:= CosIntegral[3. I]
Out[5]= 4.96039 + 1.5708*I

In[6]:= D[CosIntegral[x], x]
Out[6]= Cos[x]/x

In[7]:= Series[CosIntegral[x], {x, 0, 6}]
Out[7]= EulerGamma + Log[x] - 1/4 x^2 + 1/96 x^4 - 1/4320 x^6 + O[x]^7

In[8]:= Normal[Series[CosIntegral[x], {x, Infinity, 3}]]
Out[8]= -Cos[x]/x^2 + Sin[x] (1/x - 2/x^3)
```

## Algorithm

```text
 Mathilda -- the cosine integral  Ci(z) = -Int_z^Inf Cos[t]/t dt.

  CosIntegral[z]   Ci(z)
```

Ci has a logarithmic singularity at z = 0 and a branch cut running along the negative real axis (-Infinity, 0]. It is NOT entire and NOT odd, so -- unlike SinIntegral -- there is no odd-symmetry fold; instead negative / left-half-plane inputs are handled through the principal Log branch and an explicit reflection. Evaluation is layered so each argument takes the cheapest route:

```text
  exact special values     ->  -Infinity, 0, I Pi, Infinity, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (DLMF 6.6.5) is

```text
  Ci(z) = EulerGamma + Log(z) + Sum_{k>=1} (-1)^k z^(2k) / (2k (2k)!),
```

valid on the principal branch for all z != 0 -- the principal Log(z) already supplies the correct +-i Pi jump across the cut, so no folding is needed for the convergent path. Its partial sums can reach magnitude ~e^|z| before the O(1)-sized answer emerges, so the MPFR paths add ~|z|/ln2 guard bits to absorb that cancellation exactly. For large |z| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 6.12.4)

```text
  Ci(z) = sin(z) f(z) - cos(z) g(z)   (+ Stokes constant, see below),
    f(z) ~ Sum (-1)^k (2k)!   / z^(2k+1) = 1/z - 2!/z^3 + ...,
    g(z) ~ Sum (-1)^k (2k+1)! / z^(2k+2) = 1/z^2 - 3!/z^4 + ...,
```

(the same f, g SinIntegral computes; only the combination differs -- Si uses Pi/2 - cos f - sin g). The bare sin f - cos g asymptotic is valid in the open right half plane; on/left of the imaginary axis a piecewise-constant Stokes term restores the principal branch:

```text
  Ci(z) = sin(z) f(z) - cos(z) g(z) + i C,
    C =  0                for Re(z) > 0,
    C = +-Pi/2            for Re(z) = 0 (sign of Im z),
    C = +-Pi              for Re(z) < 0 (sign of Im z; +Pi from above on the cut).
```

(Equivalently, for Re(z) < 0 we fold w = -z into the right half plane, use the reflection Ci(z) = Ci(-z) + i Pi sign(Im z), and evaluate Ci(-z) by the clean right-half asymptotic.) For a negative real x the result is the from-above limit Complex[Ci(|x|), Pi], matching Ci(-Infinity) = I Pi.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_cosintegral` handles `CosIntegral[z] = Ci(z)`, which has a logarithmic singularity at `0` and a branch cut along `(-Infinity, 0]` (not entire, not odd). Exact special values first: `Ci[0] = -Infinity`, `Ci[Infinity] = 0`, `Ci[-Infinity] = I Pi`, `Ci[±I Infinity] = Infinity`, `ComplexInfinity`/`Indeterminate -> Indeterminate`. A **numeric real or complex** argument routes to the MPFR kernel: for moderate `|z|` the **convergent Maclaurin series** DLMF 6.6.5 (`Ci(z) = γ + Log(z) + Sum_{k>=1} (-1)^k z^{2k}/(2k(2k)!)`, with `|z|/ln2` guard bits), whose principal `Log(z)` supplies the correct `±iπ` jump across the cut with no folding; for large `|z|` the **asymptotic expansion** DLMF 6.12.4 (`sin(z) f(z) - cos(z) g(z)`) plus a piecewise Stokes/reflection constant (`0`, `±π/2`, or `±π`) restoring the principal branch (a negative real `x` gives the from-above value `Complex[Ci(|x|), Pi]`). The complex path uses the shared `ncpx` toolkit. A `USE_MPFR=0` build uses a machine-double series/asymptotic.

**Data structures.** `Expr`; the shared complex-MPFR toolkit `ncpx` (`numeric_complex.h`). The ND kernel is a real `REG_U` registration (`NDKU_CosIntegral`): real buffers via `ndk_CosIntegral_r` → `sf_machine_ci`, complex buffers via `cosintegral_machine_complex` (the file's own machine series, with a shared `sf_series_usable` cancellation gate that declines rather than returns garbage). `Compile[]` lowers `CosIntegral` at scalar (`Compiled -> True`, `ResultType -> Real`) and rank-1 array shapes.

**Complexity / limits.** `O(1)` per element at machine precision; MPFR term count and guard bits scale with `|z|` and precision. Logarithmic singularity at `0`, branch cut on the negative real axis. Symbolic arguments stay symbolic. Attributes: `Listable`, `NumericFunction`, `Protected`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Log](../../elementary-functions/Log/)

- DLMF §6.2.11, §6.6.5 — the cosine integral Ci(z) and its convergent Maclaurin series.
- DLMF §6.12.2 — the asymptotic expansion Ci(z) ~ sin(z) f(z) - cos(z) g(z).
- Source: [`src/special_functions/cosintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/cosintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_cosintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cosintegral.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_series.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series.c)

## Notes & additional examples

### Notes

`CosIntegral[z]` is the cosine integral `Ci(z) = -Integral_z^Infinity Cos[t]/t dt`.
Unlike its sibling [`SinIntegral`](SinIntegral.md) — which is entire and odd — `Ci`
has a logarithmic singularity at the origin (`CosIntegral[0] = -Infinity`) and a
branch cut running along the negative real axis `(-Infinity, 0]`. On the cut it
takes the from-above value, so for a negative real `x` the result is complex:
`Ci(x) = Ci(|x|) + I Pi`, matching `CosIntegral[-Infinity] = I Pi`. On the imaginary
axis `Ci(I y) = Chi(y) + I Pi/2` in terms of the hyperbolic cosine integral. Its
derivative is `Cos[z]/z`. Numeric evaluation uses a convergent series near the
origin and a trig-prefactored asymptotic expansion for large `|z|`, at machine or
arbitrary (MPFR) precision. Listable.
