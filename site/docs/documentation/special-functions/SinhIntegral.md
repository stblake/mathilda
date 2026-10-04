# SinhIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SinhIntegral[z]`**

gives the hyperbolic sine integral Shi(z) = Integral\_0^z Sinh\[t\]/t dt.

**`SinhIntegral[+-Infinity] = +-Infinity, SinhIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. SinhIntegral\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[SinhIntegral\[z\], z\] = Sinh\[z\]/z. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= SinhIntegral[0]
Out[1]= 0
```

Odd in its argument

```mathematica
In[2]:= SinhIntegral[-2 x]
Out[2]= -SinhIntegral[2 x]
```

```mathematica
In[3]:= SinhIntegral[Infinity]
Out[3]= Infinity
```

The integrand of the defining integral

```mathematica
In[4]:= D[SinhIntegral[x], x]
Out[4]= Sinh[x]/x
```

```mathematica
In[5]:= Series[SinhIntegral[x], {x, 0, 7}]
Out[5]= x + 1/18 x^3 + 1/600 x^5 + 1/35280 x^7 + O[x]^8

In[6]:= N[SinhIntegral[1], 20]
Out[6]= 1.05725087537572851458
```

A complex argument

```mathematica
In[7]:= N[SinhIntegral[2.0 + 1.0 I]]
Out[7]= 2.03968 + 1.67824*I
```

## Algorithm

```text
 Mathilda -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.

  SinhIntegral[z]   Shi(z)
```

Shi is entire and odd, with no branch cuts. It is the imaginary-axis sibling of the sine integral: Shi(z) = -i Si(i z). Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Shi(z) = Sum_{k>=0} z^(2k+1) / ((2k+1) (2k+1)!),
```

valid for all z. Unlike Si, every term is positive, so on the real axis there is NO catastrophic cancellation -- the partial sums climb monotonically to the O(e^|z|)-sized answer and only a small fixed guard is needed. (Complex inputs still cancel like Si, so the complex paths keep the ~|z|/ln2 guard.) For large

```text
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Shi(z) ~ cosh(z) F(z) + sinh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

summed to the smallest term (optimal truncation). F, G are Si's f, g without the (-1)^k. The odd-symmetry reduction Shi(-z) = -Shi(z) folds negative real / left-half-plane inputs onto Re >= 0, keeping the asymptotic within its valid sector. On the real axis the bare form above is exact (real, no constant); off it, matching Mathematica's Series[SinhIntegral[z],{z,Infinity,k}], a Stokes constant restores the analytic value:

```text
  Shi(z) = B(z) + i (Pi/2) sign(Im z),   B(z) = cosh(z) F(z) + sinh(z) G(z).
```

(B is odd; on the imaginary axis B ~ O(1/z) vanishes and the i Pi/2 sign(Im z) is all that survives, giving Shi(+- i Infinity) = +- i Pi/2.)

Machine-real results that overflow a C double (e.g. SinhIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_sinhintegral` evaluates `Shi(z) = Int_0^z Sinh[t]/t dt`,
entire and odd, the imaginary-axis sibling of `Si` (`Shi(z) = -i Si(i z)`). Exact
special values: `0 -> 0`, `±Infinity`, `±I Pi/2` (directed imaginary infinity),
`ComplexInfinity`/`Indeterminate -> Indeterminate`; a negative-leading `Times`
folds by odd symmetry. Numeric (MPFR): the convergent Maclaurin series — whose
real-axis terms are all positive, so there is **no** catastrophic cancellation
and only a small fixed guard is needed (complex inputs keep the `~ |z|/ln2`
guard); and for large `|z|` the asymptotic `cosh(z) F(z) + sinh(z) G(z)` summed
to its smallest term, with the Stokes constant `i (Pi/2) sign(Im z)` restored
off the real axis. A machine-real result that overflows a `double`
(`Shi` grows like `e^|x|`) is emitted as a 53-bit MPFR real. The machine kernel
`sinhintegral_machine_complex` carries a `sf_series_usable` cancellation gate.

**Data structures.** `Expr`; `mpfr_t` (real) and `ncpx` (`mpfr_t` re/im,
complex); a `double complex` fallback for `USE_MPFR=0` builds. ND: unary kernel
`NDKU_SinhIntegral = { sinhintegral_machine_complex, ndk_SinhIntegral_r, ... }`
(the real kernel is `sf_machine_shi`), registered `REG_U`, so `packed_aware`.
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Convergent regime needs no cancellation guard on the
real axis; off-axis and asymptotic regimes behave like `Si`.
`D[SinhIntegral[z], z] = Sinh[z]/z`. `Compile[]` lowers at both scalar and
rank-1 array shapes (`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [SinIntegral](../../special-functions/SinIntegral/)

- DLMF §6.2 — the hyperbolic sine integral Shi.
- Source: [`src/special_functions/sinhintegral.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/sinhintegral.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sinhintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sinhintegral.c)

## Notes & additional examples

### Notes

`SinhIntegral[z]` is the hyperbolic sine integral `Shi(z) = Int_0^z Sinh[t]/t
dt`, an entire, odd function and the imaginary-axis sibling of `SinIntegral`
(`Shi(z) = -I SinIntegral[I z]`). Its derivative is `Sinh[z]/z`, the integrand
of its defining integral.

Numeric evaluation sums the Maclaurin series for moderate `|z|` and switches to
an asymptotic expansion for large `|z|`; because every real-axis term is
positive there is no cancellation, and a real result that overflows a machine
double (`Shi` grows like `e^|z|`) is returned at extended exponent range. The
special values `SinhIntegral[±Infinity] = ±Infinity` and
`SinhIntegral[±I Infinity] = ±I Pi/2` are returned directly.
