# NSeries

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NSeries[f, {x, x0, n}]`**

gives a numerical approximation to the series expansion of f about x = x0, including the terms (x - x0)^-n through (x - x0)^n, as a SeriesData object.

<details>
<summary>Notes</summary>

f is sampled on a circle in the complex plane centred at x0 and a discrete Fourier transform of the samples recovers the Taylor or Laurent coefficients (Cauchy's integral formula). The region of convergence is the annulus, containing the sampled circle, where f is analytic. Works for essential singularities (e.g. Sin\[x + 1/x\]) where the symbolic Series cannot. Returns an incorrect result if the disk centred at x0 contains a branch cut of f; for a Laurent series the SeriesData neglects higher-order poles. No effort is made to justify the precision of the coefficients, and small spurious residuals are not recognised as zero. The sample count is grown adaptively (doubling the DFT size) until the coefficients settle within AccuracyGoal, and a NSeries::accgl warning is issued if the goal is not reached; Chop the result when needed. Options: Radius (radius of the sampled circle, default 1), WorkingPrecision (default MachinePrecision), AccuracyGoal (default MachinePrecision), PrecisionGoal.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= NSeries[Exp[x], {x, 0, 5}] // Chop
Out[1]= 1.0 + 1.0 x + 0.5 x^2 + 0.166667 x^3 + 0.0416667 x^4 + 0.00833333 x^5 + O[x]^6

In[2]:= NSeries[Exp[x], {x, I, 5}] // Chop
Out[2]= 0.540302 + 0.841471*I + (0.540302 + 0.841471*I) (x - I) + (0.270151 + 0.420735*I) (x - I)^2 + (0.0900504 + 0.140245*I) (x - I)^3 + (0.0225126 + 0.0350613*I) (x - I)^4 + (0.00450252 + 0.00701226*I) (x - I)^5 + O[x - I]^6

In[3]:= NSeries[Sin[x + 1/x], {x, 0, 10}] // Chop
Out[3]= 2.49234e-06/x^9 - 0.000174944/x^7 + 0.00703963/x^5 - 0.128943/x^3 + 0.576725/x + 0.576725 x - 0.128943 x^3 + 0.00703963 x^5 - 0.000174944 x^7 + 2.49234e-06 x^9 + O[x]^11
```

### Options (2)

```mathematica
In[4]:= NSeries[1/((1 + x) (3 + x)), {x, 0, 10}, Radius -> 5] // Chop
Out[4]= (9841.0 + 1.4949e-09*I)/x^10 + (-3280.0 - 6.31889e-10*I)/x^9 + (1093.0 + 1.09245e-10*I)/x^8 - 364.0/x^7 + 121.0/x^6 - 40.0/x^5 + 13.0/x^4 - 4.0/x^3 + 1.0/x^2 + O[x]^11

In[5]:= NSeries[Exp[x], {x, 0, 5}, WorkingPrecision -> 30] // Chop
Out[5]= 0.9999999999999999999999999999992 + 0.9999999999999999999999999999984 x + 0.4999999999999999999999999999992 x^2 + 0.1666666666666666666666666666657 x^3 + 0.0416666666666666666666666666655 x^4 + 0.008333333333333333333333333331009 x^5 + O[x]^6
```

### Applications (2)

```mathematica
In[1]:= Chop[NSeries[Cos[x], {x, 0, 6}]]
Out[1]= 1.0 - 0.5 x^2 + 0.0416667 x^4 - 0.00138889 x^6 + O[x]^7
```

```mathematica
In[1]:= Chop[NSeries[Sin[x + 1/x], {x, 0, 3}]]
Out[1]= -0.128943/x^3 + 0.576725/x + 0.576725 x - 0.128943 x^3 + O[x]^4
```

## Algorithm

nseries.c — NSeries[f, {x, x0, n}, opts]

Numerical Taylor/Laurent series expansion of `f` about x = x0, including the terms (x - x0)^-n through (x - x0)^n, returned as a SeriesData object.

```text
  NSeries[f, {x, x0, n}]
```

Method (Lyness & Sande, 1971; Bornemann, FoCM 2011). f is sampled at N equispaced points on a circle of radius r centred at x0,

```text
    z_j = x0 + r e^(2 pi i j / N),   j = 0 .. N-1,
```

and a discrete Fourier transform of the samples recovers the Laurent coefficients via Cauchy's integral formula:

```text
    c_k = (1/N) sum_j f(z_j) e^(-2 pi i j k / N)        (DFT of the samples)
    a_e = c_(e mod N) * r^(-e)         for e = -n .. n  (coeff of (x-x0)^e)
```

The upper-half DFT bins (k = N-m) supply the NEGATIVE-power coefficients, so one transform yields both the principal part and the analytic part. This is exact (no truncation) when f is analytic on an annulus containing the circle; it fails if the disk centred at x0 contains a branch cut of f.

N is chosen as a power of two with an oversampling margin,

```text
    N = 2^(ceil(log2 n) + 2),
```

so the leading aliased term a_(k +/- N) r^(+/- N) is pushed below the round-off floor. Because round-off breaks the conjugate symmetry of real-coefficient functions, the result carries tiny spurious residuals — Chop the result when needed.

A direct O(N^2) DFT is used rather than an FFT: N is small (<= a few hundred) and each sample requires a full symbolic evaluation of f, which dominates the runtime by orders of magnitude. The same code path serves both the machine (double _Complex) and arbitrary-precision (MPFR) computations; no double-precision FFT library could serve the MPFR path anyway.

Options (trailing Rule[...] in any order):

```text
  Radius           -> r                contour radius (default 1)
  WorkingPrecision -> MachinePrecision | digits
```

Memory: receives `res` owned by the evaluator. Returns a fresh Expr* (a SeriesData) on success or NULL (unevaluated). Never frees `res`. All temporary OwnValues are removed before returning, on every path.

## Implementation notes

**Attributes:** `Protected`.

## See also

[SeriesData](../../power-series/SeriesData/), [Series](../../power-series/Series/), [Chop](../../elementary-functions/Chop/), [N](../../arithmetic/N/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_nseries.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nseries.c)

## Notes & additional examples

### Notes

`NSeries[f, {x, x0, n}]` recovers Taylor or Laurent coefficients by sampling `f`
on a circle in the complex plane around `x0` and taking a discrete Fourier
transform (Cauchy's integral formula). The `Cos` example reproduces the familiar
machine-precision coefficients `1, -1/2, 1/24, -1/720`. The second case is the
function's headline capability: `Sin[x + 1/x]` has an essential singularity at
the origin, so the ordinary symbolic `Series` cannot expand it, yet `NSeries`
returns its full Laurent expansion. The coefficients are Bessel values
(`0.576725 = BesselJ[1, 2]`, the `±1` terms, and `-0.128943` for the `±3`
terms). `Chop` clears the small spurious residuals from the numerical transform.
