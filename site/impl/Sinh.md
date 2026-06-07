---
source: src/hyperbolic.c
---
**Algorithm.** `builtin_sinh` is the hyperbolic analogue of `builtin_sin`: `strip_inverse_call(arg, "ArcSinh")` for `Sinh[ArcSinh[x]] -> x`; `try_simp_forward_of_inverse_hyp` for `Sinh` of the other inverse hyperbolics (`Sinh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]`, `Sinh[ArcTanh[x]] -> x/Sqrt[1-x^2]`); `odd_fold` for `Sinh[-x] -> -Sinh[x]`; `hyp_i_fold(arg, "Sin", +1)` for `Sinh[I y] -> I Sin[y]`. Special points: `Sinh[0] = 0`, `Sinh[Infinity] = Infinity`, `Sinh[-Infinity] = -Infinity`.

**Numeric.** MPFR values evaluate via `numeric_mpfr_apply_unary(..., mpfr_sinh)` with an `mpfr_complex_sinh` complex fallback; otherwise `get_approx` + `csinh` yields a real or `Complex` result for inexact arguments. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.
