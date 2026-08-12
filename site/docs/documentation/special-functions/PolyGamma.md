# PolyGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolyGamma[z]`**

gives the digamma function psi(z) (rewritten as PolyGamma\[0, z\]).

**`PolyGamma[n, z]`**

gives the n-th derivative of the digamma function, psi^(n)(z).

**`LogGamma[z]. Listable.`**

<details>
<summary>Notes</summary>

Positive-integer arguments reduce to exact values: psi(m) to a rational minus EulerGamma, and psi^(n)(m) for odd n to a rational plus a rational multiple of Pi^(n+1); even orders stay symbolic. Non-positive integer arguments give ComplexInfinity. Inexact real and complex arguments evaluate numerically at machine or arbitrary (MPFR) precision. PolyGamma\[-1, z\] gives

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

The digamma function at positive integers reduces to an exact rational minus
Euler's constant:

```mathematica
In[1]:= PolyGamma[1]
Out[1]= -EulerGamma

In[2]:= PolyGamma[5]
Out[2]= 25/12 - EulerGamma
```

Higher derivatives `psi^(n)` at integer points give exact closed forms. The
trigamma at `1` is the Basel constant, and the odd-order values are rational
multiples of even powers of `Pi`:

```mathematica
In[1]:= PolyGamma[1, 1]
Out[1]= 1/6 Pi^2

In[2]:= PolyGamma[3, 1]
Out[2]= 1/15 Pi^4
```

The numeric paths cover arbitrary precision and complex arguments — here the
tetragamma value `psi^(2)(1+i)` to 30 digits:

```mathematica
In[1]:= N[PolyGamma[0, 3/2], 40]
Out[1]= 0.036489973978576520559023667001244432806843

In[2]:= N[PolyGamma[2, 1 + I], 30]
Out[2]= 0.3685529315879351717366345429807 + 0.7666528503450662124026953776316*I
```

The order `-1` is the special "integral of psi" case, the log-gamma function:

```mathematica
In[1]:= PolyGamma[-1, z]
Out[1]= LogGamma[z]
```

## Algorithm

Mathilda -- the PolyGamma function family.

```text
  PolyGamma[z]      digamma psi(z) = Gamma'(z)/Gamma(z)
                      (always rewrites to the two-argument form PolyGamma[0, z])
  PolyGamma[n, z]   n-th polygamma psi^(n)(z) = d^n/dz^n psi(z)
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route available:

```text
  order n = -1                -> LogGamma[z]   (inert; psi^(-1) = log-gamma)
  z a non-positive integer    -> ComplexInfinity (pole of every psi^(n))
  z a positive integer        -> exact closed form:
                                   n = 0          : H_{z-1} - EulerGamma
                                   n >= 1 odd     : rational + rational*Pi^(n+1)
                                                    (via zeta(n+1), even -> Bernoulli)
                                   n >= 1 even    : stays symbolic (zeta(odd), no
                                                    closed form)
  z inexact real (Real/MPFR)  -> numeric: mpfr_digamma for n = 0, otherwise the
                                   recurrence-shift + Bernoulli asymptotic series
  z inexact complex           -> the same asymptotic, in complex arithmetic
  everything else             -> stays symbolic (return NULL)
```

MPFR provides a digamma (n = 0) but no higher polygamma and no Hurwitz zeta, so the n >= 1 numeric kernel and the complex digamma are implemented here from the classical asymptotic expansion (Abramowitz & Stegun 6.4.11):

```text
  psi^(n)(w) ~ (-1)^(n-1) [ (n-1)!/w^n + n!/(2 w^(n+1))
                            + Sum_{k>=1} B_{2k} (2k+n-1)!/(2k)! w^(-(2k+n)) ]   (n>=1)
  psi(w)     ~ ln w - 1/(2w) - Sum_{k>=1} B_{2k}/(2k) w^(-2k)                   (n=0)
```

valid for large |w|; the recurrence psi^(n)(z) = psi^(n)(z+1) - (-1)^n n! z^(-(n+1)) shifts a small or negative argument up into that regime.

Attributes: Listable, NumericFunction, Protected.

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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[EulerGamma](../../mathematical-constants/EulerGamma/), [Series](../../power-series/Series/), [D](../../calculus/D/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)

## Notes & additional examples

### Notes

`PolyGamma[z]` is the digamma function ψ(z), stored internally as
`PolyGamma[0, z]`; `PolyGamma[n, z]` is its n-th derivative ψ⁽ⁿ⁾(z). Positive
integer arguments reduce exactly: ψ(m) to a rational minus `EulerGamma`, and
ψ⁽ⁿ⁾(m) for odd `n` to a rational plus a rational multiple of `Pi^(n+1)` (even
orders are left symbolic). Non-positive integer arguments give `ComplexInfinity`
(the poles). Inexact real and complex arguments evaluate numerically at machine
or MPFR precision, and `PolyGamma[-1, z]` returns `LogGamma[z]`. Listable.
