# Sinh

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sinh[z]`**

gives the hyperbolic sine of z, (Exp\[z\] - Exp\[-z\]) / 2.

<details>
<summary>Notes</summary>

Sinh is Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Sinh[0]
Out[1]= 0

In[2]:= N[Sinh[1]]
Out[2]= 1.1752

In[3]:= Sinh[-x]
Out[3]= -Sinh[x]

In[4]:= Sinh[ArcSinh[x]]
Out[4]= x

In[5]:= Sinh[I x]
Out[5]= I Sin[x]

In[6]:= TrigExpand[Sinh[x + y]]
Out[6]= Cosh[x] Sinh[y] + Sinh[x] Cosh[y]

In[7]:= N[Sinh[1], 40]
Out[7]= 1.1752011936438014568823818505956008151557
```

## Implementation notes

**Algorithm.** `builtin_sinh` is the hyperbolic analogue of `builtin_sin`: `strip_inverse_call(arg, "ArcSinh")` for `Sinh[ArcSinh[x]] -> x`; `try_simp_forward_of_inverse_hyp` for `Sinh` of the other inverse hyperbolics (`Sinh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]`, `Sinh[ArcTanh[x]] -> x/Sqrt[1-x^2]`); `odd_fold` for `Sinh[-x] -> -Sinh[x]`; `hyp_i_fold(arg, "Sin", +1)` for `Sinh[I y] -> I Sin[y]`. Special points: `Sinh[0] = 0`, `Sinh[Infinity] = Infinity`, `Sinh[-Infinity] = -Infinity`.

**Numeric.** MPFR values evaluate via `numeric_mpfr_apply_unary(..., mpfr_sinh)` with an `mpfr_complex_sinh` complex fallback; otherwise `get_approx` + `csinh` yields a real or `Complex` result for inexact arguments. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/hyperbolic.c`](https://github.com/stblake/mathilda/blob/main/src/hyperbolic.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_besseli.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besseli.c)
- Tests: [`tests/test_complexexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complexexpand.c)
- Tests: [`tests/test_coshintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coshintegral.c)

## Notes & additional examples

### Notes

`Sinh[z]` is the hyperbolic sine, `(Exp[z] - Exp[-z])/2`. It is odd, so negative arguments pull the sign out front. `Sinh` is Listable.
