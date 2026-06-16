# Sinh

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sinh[z]
    gives the hyperbolic sine of z, (Exp[z] - Exp[-z]) / 2.
Sinh is Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_sinh` is the hyperbolic analogue of `builtin_sin`: `strip_inverse_call(arg, "ArcSinh")` for `Sinh[ArcSinh[x]] -> x`; `try_simp_forward_of_inverse_hyp` for `Sinh` of the other inverse hyperbolics (`Sinh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]`, `Sinh[ArcTanh[x]] -> x/Sqrt[1-x^2]`); `odd_fold` for `Sinh[-x] -> -Sinh[x]`; `hyp_i_fold(arg, "Sin", +1)` for `Sinh[I y] -> I Sin[y]`. Special points: `Sinh[0] = 0`, `Sinh[Infinity] = Infinity`, `Sinh[-Infinity] = -Infinity`.

**Numeric.** MPFR values evaluate via `numeric_mpfr_apply_unary(..., mpfr_sinh)` with an `mpfr_complex_sinh` complex fallback; otherwise `get_approx` + `csinh` yields a real or `Complex` result for inexact arguments. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/hyperbolic.c`](https://github.com/stblake/mathilda/blob/main/src/hyperbolic.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sinh[0]
Out[1]= 0

In[2]:= N[Sinh[1]]
Out[2]= 1.1752

In[3]:= Sinh[-x]
Out[3]= -Sinh[x]

In[4]:= Sinh[ArcSinh[x]]
Out[4]= x
```

```mathematica
In[1]:= Sinh[I x]
Out[1]= I Sin[x]

In[2]:= TrigExpand[Sinh[x + y]]
Out[2]= Cosh[x] Sinh[y] + Sinh[x] Cosh[y]

In[3]:= N[Sinh[1], 40]
Out[3]= 1.1752011936438014568823818505956008151557
```

### Notes

`Sinh[z]` is the hyperbolic sine, `(Exp[z] - Exp[-z])/2`. It is odd, so negative arguments pull the sign out front. `Sinh` is Listable.
