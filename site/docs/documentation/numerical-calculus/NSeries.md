# NSeries

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NSeries[f, {x, x0, n}]
    gives a numerical approximation to the series expansion of f about x = x0, including the terms (x - x0)^-n through (x - x0)^n, as a SeriesData object.

f is sampled on a circle in the complex plane centred at x0 and a discrete Fourier transform of the samples recovers the Taylor or Laurent coefficients (Cauchy's integral formula). The region of convergence is the annulus, containing the sampled circle, where f is analytic. Works for essential singularities (e.g. Sin[x + 1/x]) where the symbolic Series cannot. Returns an incorrect result if the disk centred at x0 contains a branch cut of f; for a Laurent series the SeriesData neglects higher-order poles. No effort is made to justify the precision of the coefficients, and small spurious residuals are not recognised as zero -- Chop the result when needed.

Options: Radius (radius of the sampled circle, default 1), WorkingPrecision (default MachinePrecision).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NSeries[Exp[x], {x, 0, 5}] // Chop
Out[1]= 1.0 + 1.0 x + 0.5 x^2 + 0.166667 x^3 + 0.0416667 x^4 + 0.00833333 x^5 + O[x]^6

In[2]:= NSeries[Exp[x], {x, I, 5}] // Chop
Out[2]= 0.540302 + 0.841471*I + (0.540302 + 0.841471*I) (x - I) + (0.270151 + 0.420735*I) (x - I)^2 + (0.0900504 + 0.140245*I) (x - I)^3 + (0.0225126 + 0.0350613*I) (x - I)^4 + (0.00450252 + 0.00701226*I) (x - I)^5 + O[x - I]^6

In[3]:= NSeries[Sin[x + 1/x], {x, 0, 10}] // Chop
Out[3]= 2.49234e-06/x^9 - 0.000174944/x^7 + 0.00703963/x^5 - 0.128943/x^3 + 0.576725/x + 0.576725 x - 0.128943 x^3 + 0.00703963 x^5 - 0.000174944 x^7 + 2.49234e-06 x^9 + O[x]^11

In[4]:= NSeries[1/((1 + x) (3 + x)), {x, 0, 10}, Radius -> 5] // Chop
Out[4]= (9841.0 - 1.0932e-10*I)/x^10 + (-3280.0 - 1.39072e-10*I)/x^9 + 1093.0/x^8 - 364.0/x^7 + 121.0/x^6 - 40.0/x^5 + 13.0/x^4 - 4.0/x^3 + 1.0/x^2 + O[x]^11

In[5]:= NSeries[Exp[x], {x, 0, 5}, WorkingPrecision -> 30] // Chop
Out[5]= 0.9999999999999999999999999999984 + 0.9999999999999999999999999999968 x + 0.5 x^2 + 0.1666666666666666666666666666661 x^3 + 0.04166666666666666666666666666515 x^4 + 0.008333333333333333333333333330171 x^5 + O[x]^6
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
