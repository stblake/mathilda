# Tan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Tan[z]`**

gives the tangent of z. Equivalent to Sin\[z\] / Cos\[z\].

<details>
<summary>Notes</summary>

Tan is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= Tan[Pi/4]
Out[1]= 1

In[2]:= N[Tan[1]]
Out[2]= 1.55741

In[3]:= Tan[Pi/2]
Out[3]= ComplexInfinity
```

```mathematica
In[1]:= Tan[Pi/12]
Out[1]= 2 - Sqrt[3]
```

```mathematica
In[1]:= N[Tan[1], 40]
Out[1]= 1.5574077246549022305069748074583601730872
```

```mathematica
In[1]:= Simplify[Tan[ArcSin[x]]]
Out[1]= x/Sqrt[1 - x^2]
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NSolve all roots, degree 50 | 15.9 s | 0.907 s | 0.184 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| integrate Tan[x]^3 | 4.07 s | 0.196 s | 3.62 s |
| integrate 1/(1+Exp[x]) | 2.43 s | 0.109 s | 3.31 s |
| integrate Log[x]^3 | 1.59 s | 1.85 s | 10.6 s |

## Implementation notes

**Algorithm.** `builtin_tan` mirrors the `Sin` cascade in `src/trig.c`: `strip_inverse_call(arg, "ArcTan")` for `Tan[ArcTan[x]] -> x`; `try_simp_forward_of_inverse` for `Tan` of the other inverse trig functions (`Tan[ArcSin[x]] -> x/Sqrt[1-x^2]`, `Tan[ArcCos[x]] -> Sqrt[1-x^2]/x`, `Tan[ArcCot[x]] -> 1/x`); `odd_fold` for the odd symmetry `Tan[-x] -> -Tan[x]`; `trig_i_fold` for `Tan[I y] -> I Tanh[y]`; and `Tan[0] = 0`. Exact rational-multiple-of-`Pi` values are detected by `extract_pi_multiplier` and computed by `exact_tan` (a denominator-switch table analogous to `exact_sin`).

**Numeric.** MPFR-valued arguments go through `numeric_mpfr_apply_unary(..., mpfr_tan)` with an `mpfr_complex_tan` complex fallback; otherwise `get_approx` plus `ctan` produces a machine-precision real or `Complex` result, only when the argument is inexact. Symbolic arguments return `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_condition_downvalue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_condition_downvalue.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`Tan[z]` is equivalent to `Sin[z]/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Tan` is Listable.
