# Abs

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Abs[z] gives the absolute value (modulus) of numeric z, Sqrt[Re[z]^2 + Im[z]^2] for complex z.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Abs[-5]
Out[1]= 5

In[2]:= Abs[3 + 4 I]
Out[2]= 5

In[3]:= Abs[-3/4]
Out[3]= 3/4

In[4]:= Abs[{-1, 2, -3}]
Out[4]= {1, 2, 3}

In[5]:= Abs[(1 + I)^10]
Out[5]= 32

In[6]:= Abs[Sqrt[2] + Sqrt[3] I]
Out[6]= Sqrt[5]

In[7]:= N[Abs[Gamma[1/3 + 2 I]], 30]
Out[7]= 0.09665959425732664141022797859867
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| NI 2-D ridge | 2.54 s | 43.5 s | 0.685 s |
| NI oscillatory k=40 | 0.515 s | 4.05 s | 0.002 s |
| NI oscillatory k=200 | 0.477 s | 21.5 s | 0.002 s |
| NI oscillatory k=1000 | 0.389 s | 140 s | 0.002 s |
| NI oscillatory k=1001 nonzero | 0.382 s | 1.03 s | 0.626 s |

## Implementation notes

`builtin_abs` handles each numeric kind directly: `mpz_abs` for `EXPR_BIGINT`, sign flips for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR` (via `mpfr_abs`), and `|n|/d` for rationals. For a `Complex[re, im]` literal (or an expression `complex_decompose` splits into numeric real/imag parts) it builds the symbolic modulus `Power[Plus[re^2, im^2], 1/2]`; when MPFR components are present it folds directly through `mpfr_hypot` at the combined working precision instead, which is also numerically stable across disparate magnitudes. Symbolic arguments return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [ReIm](../../arithmetic/ReIm/), [Sign](../../arithmetic/Sign/), [Conjugate](../../arithmetic/Conjugate/), [Arg](../../arithmetic/Arg/), [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)

## Notes & additional examples

### Notes

For complex `z`, `Abs[z]` returns the modulus `Sqrt[Re[z]^2 + Im[z]^2]`. Abs is Listable, so it threads element-wise over lists. Exact arguments give exact results: `Abs[(1 + I)^10]` is `32` (since `|1 + I| = Sqrt[2]` and `(Sqrt[2])^10 = 32`), and `Abs[Sqrt[2] + Sqrt[3] I]` collapses to `Sqrt[5]`. The modulus of a complex value of a special function such as `Gamma` is evaluated to the requested precision under MPFR.
