# Beta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Beta[a, b]`**

is the Euler beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a+b).

**`Beta[z, a, b]`**

is the incomplete beta function Integral\_0^z t^(a-1) (1-t)^(b-1) dt.

**`Beta[z0, z1, a, b]`**

is the generalized incomplete beta Beta\[z1, a, b\] - Beta\[z0, a, b\].

<details>
<summary>Notes</summary>

Exact for rational arguments (a positive integer gives a rational via Pochhammer); non-positive integer poles give ComplexInfinity. Machine and arbitrary-precision (MPFR) real and complex inputs evaluate numerically. The incomplete form reduces through Hypergeometric2F1. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Beta[3, 5]
Out[1]= 1/105
```

The central value is exactly `Pi`, and rational orders fold into Gamma quotients:

```mathematica
In[1]:= Beta[1/2, 1/2]
Out[1]= Pi

In[2]:= Beta[1/3, 1/3]
Out[2]= Gamma[1/3]^2/Gamma[2/3]
```

Positive-integer orders give the reciprocal binomial relation `1/B(7, 3) = 9 C(8, 2)`:

```mathematica
In[1]:= Beta[7, 3]
Out[1]= 1/252

In[2]:= 1/Beta[7, 3] - 9 Binomial[8, 2]
Out[2]= 0
```

The three-argument incomplete beta, and arbitrary-precision numerics:

```mathematica
In[1]:= Beta[5, 2, 3]
Out[1]= 1025/12

In[2]:= N[Beta[2.5, 3.5], 30]
Out[2]= 0.036815538909255388078101134397
```

## Algorithm

beta.c -- the Euler beta function and its incomplete / generalized forms.

```text
  Beta[a, b]         = Gamma(a) Gamma(b) / Gamma(a+b)
  Beta[z, a, b]      = Int_0^z t^(a-1) (1-t)^(b-1) dt          (incomplete)
  Beta[z0, z1, a, b] = Beta[z1, a, b] - Beta[z0, a, b]         (generalized)
```

Rather than re-deriving the transcendental machinery, Beta is assembled from functions that already exist and are well tested:

```text
  - Beta[a, b] is reduced to a ratio of Gamma calls and handed back to the
    evaluator. The existing Gamma builtin already closes exact integers,
    half-integers, rationals (-> Sqrt[Pi] forms), machine / arbitrary-
    precision reals, and complex arguments, so Beta inherits all of that for
    free. Non-positive-integer poles (where a, b, or a+b hits a gamma pole)
    are detected up front: a surviving pole gives ComplexInfinity, while a
    cancelling pair of poles reduces by the finite limit of the gamma ratio.

  - Beta[z, a, b] is reduced through Hypergeometric2F1 via
      B_z(a, b) = z^a / a * 2F1(a, 1-b; a+1; z),
    which terminates to an exact closed form when b is a positive integer
    (1-b a non-positive integer) and evaluates numerically (real or complex)
    otherwise.
```

Memory: builtin_beta takes ownership of res. It returns a freshly built tree (the evaluator frees res) or NULL to leave the call unevaluated (the evaluator keeps res). It never frees res itself. Every intermediate tree is consumed by eval_and_free.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## See also

[Gamma](../../special-functions/Gamma/), [Pi](../../mathematical-constants/Pi/), [Hypergeometric2F1](../../special-functions/Hypergeometric2F1/), [PolyGamma](../../special-functions/PolyGamma/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_integrate_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_beta.c)
- Tests: [`tests/test_integrate_residue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_residue.c)

## Notes & additional examples

### Notes

`Beta[a, b] = Gamma[a] Gamma[b]/Gamma[a+b]` is the Euler beta function. `Beta[z, a, b]` is the incomplete beta integral, and `Beta[z0, z1, a, b]` the generalized incomplete form. Exact for rational arguments via Pochhammer; non-positive integer poles give `ComplexInfinity`; the incomplete form reduces through `Hypergeometric2F1`. Real and complex inputs evaluate at machine or MPFR precision. Listable.
