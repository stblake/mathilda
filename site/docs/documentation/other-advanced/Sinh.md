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
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
