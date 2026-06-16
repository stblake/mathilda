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
```

```mathematica
In[1]:= Csc[Pi/12]
Out[1]= Sqrt[2] (1 + Sqrt[3])
```

```mathematica
In[1]:= N[Csc[1], 40]
Out[1]= 1.1883951057781212162615994523745510035279
```

```mathematica
In[1]:= Series[Csc[x], {x, 0, 5}]
Out[1]= 1/x + 1/6 x + 7/360 x^3 + 31/15120 x^5 + O[x]^6
```

```mathematica
In[1]:= Csc[I]
Out[1]= -I Csch[1]
```

### Notes

`Csc[z]` is `1/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Csc` is Listable. Exact special angles are returned in closed radical form, the Laurent expansion gives the cosecant's pole at the origin, and imaginary arguments map onto the hyperbolic cosecant.
