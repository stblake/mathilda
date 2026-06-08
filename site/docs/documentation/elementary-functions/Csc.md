# Csc

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Csc[z]
    gives the cosecant of z (= 1 / Sin[z]).
Csc is Listable. Singularities at z = k Pi yield ComplexInfinity.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_csc` (`src/trig.c`) applies: (1) `strip_inverse_call` folds `Csc[ArcCsc[x]] -> x`. (2) `odd_fold` for oddness `Csc[-x] -> -Csc[x]`. (3) `trig_i_fold` rewrites `Csc[I y] -> -I Csch[y]`. (4) `Csc[0] -> ComplexInfinity`. (5) For a rational multiple of Pi (`extract_pi_multiplier`), `exact_csc` returns the closed-form surd from the table for denominators 1,2,3,4,5,6,10,12. (6) Numeric fallback: MPFR via `mpfr_csc`/`mpfr_complex_csc`, else `get_approx` + `1/csin(c)` for inexact inputs. Otherwise `NULL`. (Unlike Cos/Tan, Csc has no forward-of-inverse fold step.)

**Data structures.** `Expr*` trees via the `make_*` helpers; Pi multiples as `int64_t n, d`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Csc[Pi/6]
Out[1]= 2

In[2]:= N[Csc[1]]
Out[2]= 1.1884

In[3]:= Csc[Pi]
Out[3]= ComplexInfinity
```

### Notes

`Csc[z]` is `1/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Csc` is Listable.
