---
source: src/contfrac.c
---
**Algorithm.** `builtin_from_continued_fraction` reconstructs the value from a `List` of terms. A non-periodic list uses the standard convergent recurrence `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}` (`fcf_simple`), evaluating each step so numeric terms collapse and symbolic terms stay in convergent form; the result is `h_{n-1}/k_{n-1}`. A trailing sub-list marks the cyclic (period) block, requiring all integer terms: `fcf_periodic` builds the period's convergents in GMP, forms the quadratic `A x^2 + B x + C = 0` for the purely-periodic tail (with `A=k_{k-1}`, `B=k_{k-2}-h_{k-1}`, `C=-h_{k-2}`), solves via the discriminant with largest-square extraction (`fcf_extract_square`, trial division up to `10^6`), then applies the leading terms as a Möbius transform `(Hx+H')/(Kx+K')` and rationalises to `(P + Q Sqrt[R])/S` (`fcf_qirr_to_expr`).

**Data structures.** Convergents are GMP `mpz_t` registers; symbolic non-periodic convergents are `Expr` trees via `eval_and_free`. Output is an Integer/Rational or a rationalised quadratic-irrational expression.

**Complexity / limits.** Linear in the number of terms. The periodic path requires exact integer terms; a residual `p^2 q` with both `p, q > 10^6` is left un-reduced (astronomically unlikely for reconstructed CF data).
