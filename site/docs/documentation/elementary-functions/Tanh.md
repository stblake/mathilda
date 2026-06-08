# Tanh

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Tanh[z]
    gives the hyperbolic tangent of z, Sinh[z] / Cosh[z].
Tanh is Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_tanh` follows the same hyperbolic cascade in `src/hyperbolic.c`: `strip_inverse_call(arg, "ArcTanh")` for `Tanh[ArcTanh[x]] -> x`; `try_simp_forward_of_inverse_hyp` for `Tanh` of the other inverse hyperbolics (`Tanh[ArcSinh[x]] -> x/Sqrt[1+x^2]`, `Tanh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]/x`, `Tanh[ArcCoth[x]] -> 1/x`); `odd_fold` for `Tanh[-x] -> -Tanh[x]`; `hyp_i_fold(arg, "Tan", +1)` for `Tanh[I y] -> I Tan[y]`. Special points: `Tanh[0] = 0`, `Tanh[Infinity] = 1`, `Tanh[-Infinity] = -1`.

**Numeric.** MPFR values use `numeric_mpfr_apply_unary(..., mpfr_tanh)` (complex fallback `mpfr_complex_tanh`); otherwise `get_approx` + `ctanh` covers inexact real/complex inputs. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/hyperbolic.c`](https://github.com/stblake/mathilda/blob/main/src/hyperbolic.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Tanh[0]
Out[1]= 0

In[2]:= N[Tanh[1]]
Out[2]= 0.761594

In[3]:= Tanh[ArcTanh[z]]
Out[3]= z
```

### Notes

`Tanh[z]` is the hyperbolic tangent, `Sinh[z]/Cosh[z]`. `Tanh` is Listable.
