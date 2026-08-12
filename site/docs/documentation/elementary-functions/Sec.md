# Sec

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sec[z]`**

gives the secant of z (= 1 / Cos\[z\]).

<details>
<summary>Notes</summary>

Sec is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= Sec[Pi/3]
Out[1]= 2

In[2]:= Sec[0]
Out[2]= 1

In[3]:= N[Sec[1]]
Out[3]= 1.85082
```

Exact special values come out in closed form, including the golden-ratio-related
`Sec[Pi/5]`:

```mathematica
In[1]:= Sec[Pi/4]
Out[1]= Sqrt[2]

In[2]:= Sec[Pi/5]
Out[2]= -1 + Sqrt[5]
```

An imaginary argument folds onto the hyperbolic secant via `Sec[I z] = Sech[z]`:

```mathematica
In[1]:= Sec[I]
Out[1]= Sech[1]
```

The Maclaurin series of `Sec` exposes the secant (Euler) numbers `1, 5, 61, ...`
in its coefficients:

```mathematica
In[1]:= Series[Sec[x], {x, 0, 6}]
Out[1]= 1 + 1/2 x^2 + 5/24 x^4 + 61/720 x^6 + O[x]^7
```

High-precision evaluation is available through `N`:

```mathematica
In[1]:= N[Sec[1], 40]
Out[1]= 1.8508157176809256179117532413986501934704
```

## Implementation notes

**Algorithm.** `builtin_sec` follows the `src/trig.c` cascade but uses even symmetry: `strip_inverse_call(arg, "ArcSec")` for `Sec[ArcSec[x]] -> x`; `even_fold` for `Sec[-x] -> Sec[x]` when the argument is superficially negative; `trig_i_fold(arg, "Sech", 0)` for `Sec[I y] -> Sech[y]`; and `Sec[0] = 1`. Exact values at rational multiples of `Pi` are recognised by `extract_pi_multiplier` and produced by `exact_sec`.

**Numeric.** MPFR arguments use `numeric_mpfr_apply_unary(..., mpfr_sec)` (complex fallback `mpfr_complex_sec`); otherwise `get_approx` computes `1.0 / ccos(c)` for inexact inputs, yielding `EXPR_REAL` or `Complex`. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_integrate_derivdivides.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_derivdivides.c)
- Tests: [`tests/test_integrate_dispatch.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_dispatch.c)

## Notes & additional examples

### Notes

`Sec[z]` is `1/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Sec` is Listable.
