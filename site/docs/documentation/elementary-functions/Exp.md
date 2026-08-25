# Exp

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Exp[z]`**

gives the exponential E^z.

<details>
<summary>Notes</summary>

Exp is Listable. Exp\[0\] = 1, Exp\[Log\[x\]\] -\> x, Exp\[I Pi\] = -1. Numeric inputs route to libm / MPFR.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Exp[0]
Out[1]= 1

In[2]:= Exp[I Pi]
Out[2]= -1

In[3]:= Exp[Log[x]]
Out[3]= x

In[4]:= Exp[2 Log[x]]
Out[4]= x^2

In[5]:= D[Exp[Sin[x]], x]
Out[5]= Cos[x] E^Sin[x]

In[6]:= Series[Exp[x], {x, 0, 6}]
Out[6]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + 1/720 x^6 + O[x]^7

In[7]:= N[Exp[1], 50]
Out[7]= 2.71828182845904523536028747135266249775724709369996
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NI 50-digit Gaussian | 6.87 s | 2.11 s | 1.1 s |
| v^2.5 over 4x10^6 | 6.66 s | 3.79 s | 16.3 s |
| integrate Sin[x] Exp[x] | 6.62 s | 0.726 s | 5.47 s |
| Sin[Exp[Log[v]]] fused? | 6.52 s | 5.61 s | 42.5 s |
| integrate Exp[x]/x (non-elementary) | 5.38 s | 0.318 s | 33.5 s |
| Sqrt over 4x10^6 | 5.35 s | 0.444 s | 0.925 s |

## Implementation notes

**Algorithm.** `builtin_exp` (1-arg). Exact special values: `Exp[0] = 1`, `Exp[-Infinity] = 0`, `Exp[Infinity] = Infinity`. The notable closed-form path is **Euler's formula** for `Exp[I q Pi]`: when the argument is `Times[Complex[0,q], Pi]` with `q` an integer or rational, it rewrites to `Cos[q Pi] + I Sin[q Pi]` (which the trig kernels then evaluate to exact algebraic values). Numeric arguments go through MPFR when `USE_MPFR` — `mpfr_exp` for pure reals, the `exp(a)(cos b + i sin b)` complex helper otherwise — or `cexp` on a `double complex` for machine inputs, collapsing to a real result when the imaginary part is zero. Anything that matches none of these is returned as `Power[E, z]` (the canonical internal form for an unevaluated exponential), so `Exp` is effectively a thin front-end that produces `E^z` and lets the `Power` machinery carry the symbolic case.

**Data structures.** `Expr*` trees throughout; numeric evaluation uses `double complex` (machine) or `mpfr_t`/complex-MPFR helpers (arbitrary precision). The Euler path scans the `Times` argument list for a `Pi` factor and a single pure-imaginary `Complex[0, q]` coefficient.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/logexp.c`](https://github.com/stblake/mathilda/blob/main/src/logexp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)

## Notes & additional examples

### Notes

`Exp[z]` is `E^z`; it inverts `Log` and evaluates Euler's identity exactly.
Logarithmic arguments collapse (`Exp[2 Log[x]] -> x^2`), it differentiates and
series-expands symbolically, and numeric inputs route to libm / MPFR for
arbitrary precision. Exp is Listable.
