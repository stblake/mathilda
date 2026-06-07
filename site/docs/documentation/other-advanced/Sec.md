# Sec

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sec[z]
    gives the secant of z (= 1 / Cos[z]).
Sec is Listable. Singularities at z = Pi/2 + k Pi yield ComplexInfinity.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_sec` follows the `src/trig.c` cascade but uses even symmetry: `strip_inverse_call(arg, "ArcSec")` for `Sec[ArcSec[x]] -> x`; `even_fold` for `Sec[-x] -> Sec[x]` when the argument is superficially negative; `trig_i_fold(arg, "Sech", 0)` for `Sec[I y] -> Sech[y]`; and `Sec[0] = 1`. Exact values at rational multiples of `Pi` are recognised by `extract_pi_multiplier` and produced by `exact_sec`.

**Numeric.** MPFR arguments use `numeric_mpfr_apply_unary(..., mpfr_sec)` (complex fallback `mpfr_complex_sec`); otherwise `get_approx` computes `1.0 / ccos(c)` for inexact inputs, yielding `EXPR_REAL` or `Complex`. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sec[Pi/3]
Out[1]= 2

In[2]:= Sec[0]
Out[2]= 1

In[3]:= N[Sec[1]]
Out[3]= 1.85082
```

### Notes

`Sec[z]` is `1/Cos[z]`. Singularities at `z = Pi/2 + k Pi` yield `ComplexInfinity`. `Sec` is Listable.
