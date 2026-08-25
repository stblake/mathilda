# Cos

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Cos[z]`**

gives the cosine of z (argument in radians).

<details>
<summary>Notes</summary>

Cos is Listable. Numeric inputs route to libm / MPFR; rational multiples of Pi reduce to exact values.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= Cos[Pi/3]
Out[1]= 1/2

In[2]:= Cos[{0, Pi/3, Pi}]
Out[2]= {1, 1/2, -1}

In[3]:= Cos[ArcCos[x]]
Out[3]= x

In[4]:= Cos[Pi/12]
Out[4]= 1/4 (Sqrt[2] + Sqrt[6])

In[5]:= Cos[Pi/5]
Out[5]= 1/4 (1 + Sqrt[5])

In[6]:= TrigExpand[Cos[a + b]]
Out[6]= Cos[a] Cos[b] - Sin[a] Sin[b]

In[7]:= Series[Cos[x], {x, 0, 8}]
Out[7]= 1 - 1/2 x^2 + 1/24 x^4 - 1/720 x^6 + 1/40320 x^8 + O[x]^9

In[8]:= N[Cos[1], 40]
Out[8]= 0.54030230586813971740093660744297660373228
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NSolve all roots, degree 50 | 15.9 s | 0.907 s | 0.184 s |
| Simplify quartic-to-Cos[4x] | 2.35 s | 0.003 s | 7.45 s |
| TrigReduce product of 4 sines | 0.894 s | 0.14 s | 16.5 s |
| FullSimplify nested radical | 0.185 s | 0.003 s | 1.28 s |
| Simplify log-exp collapse | 0.059 s | 0.004 s | 1.35 s |
| Simplify trig Pythagorean | 0.048 s | 0.004 s | 3.69 s |

## Implementation notes

**Algorithm.** `builtin_cos` (`src/trig.c`) applies a fixed cascade of reductions, returning the first that fires. (1) `strip_inverse_call` folds `Cos[ArcCos[x]] -> x`. (2) `try_simp_forward_of_inverse` handles `Cos[ArcSin[x]] -> Sqrt[1-x^2]` and `Cos[ArcTan[x]] -> 1/Sqrt[1+x^2]`, building the unevaluated tree and letting the evaluator canonicalise. (3) `even_fold` uses evenness `Cos[-x] -> Cos[x]`. (4) `trig_i_fold` rewrites `Cos[I y] -> Cosh[y]`. (5) `Cos[0] -> 1`. (6) For a rational-multiple-of-Pi argument (recognised by `extract_pi_multiplier` matching `Pi` or `Times[Rational[n,d], Pi]`), `exact_cos` reduces n/d mod 2π using even symmetry, maps into `[0, π/2]` tracking a sign, and returns the closed surd form from a hardcoded table for denominators 1,2,3,4,5,6,10,12. (7) Numeric fallback: MPFR via `mpfr_cos`/`mpfr_complex_cos` for arbitrary-precision args, else `get_approx` + C99 `ccos` for already-inexact (real or complex) inputs. Anything else returns `NULL` (stays symbolic).

**Data structures.** Plain `Expr*` trees throughout; exact values are assembled with `make_times`/`make_plus`/`make_sqrt`/`make_rational` helpers. Pi multiples are carried as `int64_t n, d`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)

## Notes & additional examples

### Notes

The argument is in radians; rational multiples of `Pi` reduce to exact values while numeric inputs route to libm / MPFR. `Cos` is Listable.
