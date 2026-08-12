# Sqrt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sqrt[z]`**

represents the principal square root of z.

<details>
<summary>Notes</summary>

Sqrt is Listable. Sqrt\[z\] is canonicalised to Power\[z, 1/2\]; perfect integer / rational squares reduce to exact form, negative real inputs yield I \* Sqrt\[-x\], and numeric inputs (Real / MPFR / Complex) are evaluated directly. Branch cut along the negative real axis.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (8)

```mathematica
In[1]:= Sqrt[50]
Out[1]= 5 Sqrt[2]
```

```mathematica
In[1]:= Sqrt[-9]
Out[1]= 3*I
```

Like terms combine and products of surds simplify back to integers:

```mathematica
In[1]:= Sqrt[2] + Sqrt[8]
Out[1]= 3 Sqrt[2]

In[2]:= Sqrt[12]*Sqrt[27]
Out[2]= 18
```

Nested radicals denest under `FullSimplify`:

```mathematica
In[1]:= Sqrt[3 + 2 Sqrt[2]] // FullSimplify
Out[1]= 1 + Sqrt[2]
```

Principal value of a complex root, both symbolically and to 40 digits:

```mathematica
In[1]:= Sqrt[I]
Out[1]= (1 + I)/Sqrt[2]

In[2]:= N[Sqrt[I], 40]
Out[2]= 0.70710678118654752440084436210484903928487 + 0.70710678118654752440084436210484903928487*I
```

The Puiseux series of `Sqrt[1 + x]` is delivered exactly:

```mathematica
In[1]:= Series[Sqrt[1 + x], {x, 0, 5}]
Out[1]= 1 + 1/2 x - 1/8 x^2 + 1/16 x^3 - 5/128 x^4 + 7/256 x^5 + O[x]^6
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

`builtin_sqrt` is a thin wrapper: it rewrites `Sqrt[x]` to `Power[x, Rational[1, 2]]` (via `make_rational(1, 2)`) and returns that, letting the full `Power` machinery handle all simplification (exact perfect squares, `Sqrt[8] -> 2 Sqrt[2]` radical extraction, numeric/MPFR evaluation, infinity algebra). `Sqrt` carries `LISTABLE | NUMERICFUNCTION | PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[FactorInteger](../../number-theory/FactorInteger/)

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on square-free factorization and radical simplification.
- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on integer square roots.
- Source: [`src/power.c`](https://github.com/stblake/mathilda/blob/main/src/power.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_arc_exact.c`](https://github.com/stblake/mathilda/blob/main/tests/test_arc_exact.c)
- Tests: [`tests/test_assuming.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assuming.c)

## Notes & additional examples

### Notes

`Sqrt[n]` is `Power[n, 1/2]`, so it inherits the perfect-square extraction logic:
the largest square factor is pulled out of the radical, reducing `Sqrt[50]` to
`5 Sqrt[2]` and `Sqrt[18]` to `3 Sqrt[2]`. Perfect squares collapse fully, and a
perfect-square rational like `1/4` yields the exact rational `1/2`. Negative
arguments produce a pure imaginary result, with `Sqrt[-9]` giving `3*I`. The
surd is kept in symbolic form when no square factor remains, e.g. `Sqrt[2]`.
Denesting of compound radicals like `Sqrt[3 + 2 Sqrt[2]]` is not automatic but is
recovered by `FullSimplify`. The branch cut runs along the negative real axis,
so `N[Sqrt[I], 40]` returns the upper-half-plane principal value.
