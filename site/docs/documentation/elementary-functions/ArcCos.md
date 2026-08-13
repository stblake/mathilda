# ArcCos

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArcCos[z]`**

gives the principal inverse cosine of z, in \[0, Pi\] for real z in \[-1, 1\].

<details>
<summary>Notes</summary>

ArcCos is Listable. Branch cuts run along the real axis with |z| \> 1.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= ArcCos[0]
Out[1]= 1/2 Pi

In[2]:= ArcCos[1/2]
Out[2]= 1/3 Pi

In[3]:= ArcCos[1]
Out[3]= 0

In[4]:= N[ArcCos[0.5]]
Out[4]= 1.0472

In[5]:= ArcCos[-1/2]
Out[5]= 2/3 Pi

In[6]:= ArcCos[Sqrt[2]/2]
Out[6]= 1/4 Pi

In[7]:= N[ArcCos[2], 20]
Out[7]= 0.0 + 1.31695789692481670862*I
```

## Implementation notes

**Algorithm.** `builtin_arccos` (`src/trig.c`): (1) `arc_pi_minus_fold` applies the reflection `ArcCos[-x] -> Pi - ArcCos[x]` for superficially-negative arguments. (2) `arccos_i_fold` applies the imaginary-axis identity `ArcCos[I y] -> Pi/2 - I ArcSinh[y]`. (3) Exact inversion via `exact_arccos`, which scans n in `[0,d]` for d in {1,2,3,4,5,6,10,12}, computes `exact_cos(n,d)`, and returns `n/d * Pi` on an `expr_eq` match. (4) Numeric fallback: MPFR via `mpfr_acos`/`mpfr_complex_acos`, else `get_approx` + C99 `cacos`; for real x>1 the imaginary part is negated to match the `+i*acosh(x)` branch convention. Otherwise `NULL`.

**Data structures.** `Expr*` trees; exact inversion reuses the forward `exact_cos` table.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_arc_exact.c`](https://github.com/stblake/mathilda/blob/main/tests/test_arc_exact.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)

## Notes & additional examples

### Notes

`ArcCos[z]` gives the principal inverse cosine, in `[0, Pi]` for real `z` in `[-1, 1]`. Special exact angles are recognised: `ArcCos[-1/2]` is `2/3 Pi` and `ArcCos[Sqrt[2]/2]` is `1/4 Pi`. For `|z| > 1` the result moves off the real axis onto the branch cut, so `ArcCos[2]` is purely imaginary, `1.3169... I`. `ArcCos` is Listable.
