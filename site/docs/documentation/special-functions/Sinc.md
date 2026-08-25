# Sinc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sinc[z]`**

gives the cardinal sine Sin\[z\]/z, with Sinc\[0\] = 1.

**`D[Sinc[z], z] = Cos[z]/z - Sin[z]/z^2. Listable.`**

<details>
<summary>Notes</summary>

An entire, even function. Sinc\[+-Infinity\] = 0. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Sinc[0]
Out[1]= 1

In[2]:= Sinc[2.]
Out[2]= 0.454649

In[3]:= N[Sinc[2], 45]
Out[3]= 0.454648713412840847698009932955872421351127485

In[4]:= Sinc[1. + I]
Out[4]= 0.966711 - 0.331747 I

In[5]:= D[Sinc[x], x]
Out[5]= Cos[x]/x - Sin[x]/x^2

In[6]:= Series[Sinc[x], {x, 0, 6}]
Out[6]= 1 - x^2/6 + x^4/120 - x^6/5040 + O[x]^7
```

## Algorithm

```text
 Mathilda -- the cardinal sine  Sinc[z] = Sin[z]/z  (Sinc[0] = 1).
```

Sinc is entire and even, with a removable singularity at the origin. Each kind of argument takes the cheapest route:

```text
  exact special values   ->  1 (at 0), 0 (at +-Infinity), Indeterminate
  machine real           ->  libm sin(x)/x
  arbitrary real (MPFR)  ->  mpfr_sin(x)/x at the input precision
  complex (any prec)     ->  sin(z)/z via the shared ncpx toolkit
  everything else        ->  stays symbolic (return NULL)
```

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [SinIntegral](../../special-functions/SinIntegral/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sinc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sinc.c)

## Notes & additional examples

### Notes

`Sinc[z]` is the cardinal sine `Sin[z]/z`, with the removable singularity at the
origin filled in as `Sinc[0] = 1`. It is entire and even, and `Sinc[±Infinity] = 0`.
It appears as the derivative of the sine integral: `D[SinIntegral[z], z] = Sinc[z]`.
Numeric evaluation is at machine or arbitrary (MPFR) precision for both real and
complex arguments. Listable. See also [`SinIntegral`](SinIntegral.md).
