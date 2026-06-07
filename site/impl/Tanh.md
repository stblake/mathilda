---
source: src/hyperbolic.c
---
**Algorithm.** `builtin_tanh` follows the same hyperbolic cascade in `src/hyperbolic.c`: `strip_inverse_call(arg, "ArcTanh")` for `Tanh[ArcTanh[x]] -> x`; `try_simp_forward_of_inverse_hyp` for `Tanh` of the other inverse hyperbolics (`Tanh[ArcSinh[x]] -> x/Sqrt[1+x^2]`, `Tanh[ArcCosh[x]] -> Sqrt[x-1] Sqrt[x+1]/x`, `Tanh[ArcCoth[x]] -> 1/x`); `odd_fold` for `Tanh[-x] -> -Tanh[x]`; `hyp_i_fold(arg, "Tan", +1)` for `Tanh[I y] -> I Tan[y]`. Special points: `Tanh[0] = 0`, `Tanh[Infinity] = 1`, `Tanh[-Infinity] = -1`.

**Numeric.** MPFR values use `numeric_mpfr_apply_unary(..., mpfr_tanh)` (complex fallback `mpfr_complex_tanh`); otherwise `get_approx` + `ctanh` covers inexact real/complex inputs. Symbolic input returns `NULL`. Attributes: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.
