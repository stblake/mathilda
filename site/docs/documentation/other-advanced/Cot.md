# Cot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cot[z]
    gives the cotangent of z. Equivalent to Cos[z] / Sin[z].
Cot is Listable. Singularities at z = k Pi yield ComplexInfinity.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_cot` (`src/trig.c`) runs the same cascade as the other trig heads. (1) `strip_inverse_call` folds `Cot[ArcCot[x]] -> x`. (2) `try_simp_forward_of_inverse` handles `Cot[ArcTan[x]] -> 1/x`. (3) `odd_fold` uses oddness `Cot[-x] -> -Cot[x]`. (4) `trig_i_fold` rewrites `Cot[I y] -> -I Coth[y]`. (5) `Cot[0] -> ComplexInfinity`. (6) For a rational multiple of Pi (via `extract_pi_multiplier`), `exact_cot` reduces n/d mod π into `[0, π/2]` with sign tracking and returns the table value (denominators 1,2,3,4,5,6,10,12; `ComplexInfinity` at multiples of π). (7) Numeric fallback: MPFR via `mpfr_cot`/`mpfr_complex_cot`, else `get_approx` + `1/ctan(c)`. Otherwise `NULL`.

**Data structures.** `Expr*` trees built with the `make_*` helpers; Pi multiples carried as `int64_t n, d`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cot[Pi/4]
Out[1]= 1

In[2]:= N[Cot[1]]
Out[2]= 0.642093

In[3]:= Cot[Pi]
Out[3]= ComplexInfinity
```

### Notes

`Cot[z]` is equivalent to `Cos[z]/Sin[z]`. Singularities at `z = k Pi` yield `ComplexInfinity`. `Cot` is Listable.
