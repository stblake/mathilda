# Sin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sin[z]`**

gives the sine of z (argument in radians).

<details>
<summary>Notes</summary>

Sin is Listable. Numeric inputs are evaluated via libm (Real) or MPFR (arbitrary precision); rational multiples of Pi reduce to exact values.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (9)

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2

In[2]:= N[Sin[1]]
Out[2]= 0.841471

In[3]:= Sin[{0, Pi/6, Pi/2}]
Out[3]= {0, 1/2, 1}

In[4]:= Sin[ArcSin[x]]
Out[4]= x

In[5]:= Sin[Pi/10]
Out[5]= 1/4 (-1 + Sqrt[5])

In[6]:= Sin[Pi/12]
Out[6]= 1/4 (Sqrt[6] - Sqrt[2])

In[7]:= Sin[I]
Out[7]= I Sinh[1]

In[8]:= TrigExpand[Sin[3 x]]
Out[8]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[9]:= N[Sin[1], 40]
Out[9]= 0.84147098480789650665250232163029899962254
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Interpolation evaluate, 20000 points | 18.9 s | 18.1 s | 0.073 s |
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| v^2.5 over 4x10^6 | 6.66 s | 3.79 s | 16.3 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| Sin[Exp[Log[v]]] fused? | 6.52 s | 5.61 s | 42.5 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |

## Implementation notes

**Algorithm.** `builtin_sin` is a single-argument cascade tried in order: (1) `strip_inverse_call(arg, "ArcSin")` collapses `Sin[ArcSin[x]] -> x`; (2) `try_simp_forward_of_inverse` rewrites `Sin` of the *other* inverse trig functions to radical forms (`Sin[ArcCos[x]] -> Sqrt[1-x^2]`, `Sin[ArcTan[x]] -> x/Sqrt[1+x^2]`); (3) `odd_fold` uses the odd symmetry `Sin[-x] -> -Sin[x]` whenever `expr_is_superficially_negative(arg)`; (4) `trig_i_fold` extracts an imaginary unit, `Sin[I y] -> I Sinh[y]`; (5) `Sin[0] = 0`. Exact special values come from `extract_pi_multiplier`, which detects `Pi` or `Times[Rational[n,d], Pi]`, handing `(n,d)` to `exact_sin`.

`exact_sin` reduces the angle into `[0, Pi/2]` (tracking a sign through the `[0,2Pi)` and `[0,Pi]` foldings), reduces the fraction by `gcd`, and switches on the denominator: closed radical forms are tabulated for `d` in {1,2,3,4,5,6,10,12}, e.g. `d==3 -> Sqrt[3]/2`, `d==5`/`d==10` give the golden-ratio nested-radical values, `d==12` gives `(Sqrt[6]±Sqrt[2])/4`. Anything outside the table returns `NULL` and stays symbolic.

**Numeric.** With MPFR built, an MPFR-valued argument is evaluated via `numeric_mpfr_apply_unary(..., mpfr_sin)`, falling back to `mpfr_complex_sin` for complex MPFR values, preserving the input precision. Otherwise `get_approx` extracts a `double complex` and, only if the value is genuinely inexact, applies `csin`; a real input yields `EXPR_REAL`, a complex one a `Complex[Real, Real]`.

`Sin` is registered with `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`; element-wise threading over lists is handled generically by the evaluator.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)

## Notes & additional examples

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm (Real) or MPFR. `Sin` is Listable, so it threads over lists element-wise.
