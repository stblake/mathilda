# ArcSin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ArcSin[z]
    gives the principal inverse sine of z, in [-Pi/2, Pi/2] for real z
    in [-1, 1].
ArcSin is Listable. Branch cuts run along the real axis with |z| > 1.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_arcsin` (`src/trig.c`): (1) `odd_fold` uses oddness `ArcSin[-x] -> -ArcSin[x]`. (2) `trig_i_fold` applies the principal-branch identity `ArcSin[I y] -> I ArcSinh[y]`. (3) Exact inversion via `exact_arcsin`, which brute-forces the forward table: for each denominator d in {1,2,3,4,5,6,10,12} and each n it computes `exact_sin(n,d)` and, if it `expr_eq`-matches the argument, returns `n/d * Pi`. (4) Numeric fallback: MPFR via `mpfr_asin`/`mpfr_complex_asin`, else `get_approx` + C99 `casin` for inexact inputs. On the real branch cut x>1 the imaginary part is negated to match Mathematica's lower-side convention (C99 lands on the upper side). Otherwise `NULL`.

**Data structures.** `Expr*` trees; the exact-inversion search reuses the forward `exact_sin` table rather than a separate inverse table.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ArcSin[1/2]
Out[1]= 1/6 Pi

In[2]:= ArcSin[1]
Out[2]= 1/2 Pi

In[3]:= N[ArcSin[0.5]]
Out[3]= 0.523599

In[4]:= ArcSin[Sin[x]]
Out[4]= ArcSin[Sin[x]]
```

### Notes

`ArcSin[z]` gives the principal inverse sine, in `[-Pi/2, Pi/2]` for real `z` in `[-1, 1]`. The inverse-of-forward composition `ArcSin[Sin[x]]` is deliberately *not* folded to `x`, since that holds only on the principal branch. `ArcSin` is Listable.
