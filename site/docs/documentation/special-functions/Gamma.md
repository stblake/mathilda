# Gamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Gamma[z]
    is the Euler gamma function Gamma(z).
Gamma[a, z]
    is the upper incomplete gamma function Gamma(a, z).
Gamma[a, z0, z1]
    is the generalized incomplete gamma Gamma(a, z0) - Gamma(a, z1).
Integer and half-integer arguments reduce to exact values ((z-1)!,
and rational multiples of Sqrt[Pi]); non-positive integers give
ComplexInfinity. Machine and arbitrary-precision (MPFR) real inputs
evaluate numerically, as do machine-precision complex inputs. Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact reductions for `Gamma[z]`:
  - Positive integers: `Gamma[n] = (n-1)!` (exact, with GMP BigInt for large `n`).
  - Non-positive integers are poles: `Gamma[0]`, `Gamma[-n]` → `ComplexInfinity`.
  - Half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Gamma[5]
Out[1]= 24
```

```mathematica
In[1]:= Gamma[7/2]
Out[1]= 15/8 Sqrt[Pi]
```

```mathematica
In[1]:= Gamma[-1/2]
Out[1]= -2 Sqrt[Pi]
```

```mathematica
In[1]:= N[Gamma[1/3], 40]
Out[1]= 2.6789385347077476336556929409746776441289
```

```mathematica
In[1]:= N[Gamma[3 + 4 I], 20]
Out[1]= 0.00522553847136921419473 - 0.172547079294300187719*I
```

### Notes

`Gamma` is the Euler gamma function, the analytic continuation of the
factorial: `Gamma[n] = (n-1)!`, so `Gamma[5] = 24`. Half-integer arguments
collapse to exact rational multiples of `Sqrt[Pi]` — `Gamma[7/2] = 15/8 Sqrt[Pi]`
— and this continues through the poles at non-positive integers into negative
half-integers, where `Gamma[-1/2] = -2 Sqrt[Pi]`. For arguments with no
closed form it evaluates to arbitrary precision via MPFR (`Gamma[1/3]` to 40
digits) and across the complex plane (`Gamma[3 + 4 I]`). The two- and
three-argument forms give the upper incomplete and generalized incomplete gamma
functions.
