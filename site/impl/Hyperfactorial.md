---
source: src/special_functions/hyperfactorial.c
references:
  - "DLMF §5.17 — Barnes G-function (the analytic continuation Hyperfactorial[z] = Gamma[z+1]^z / BarnesG[z+1])."
---
**Algorithm.** `builtin_hyperfactorial` computes `H(n) = prod_{k=1}^{n} k^k`
(with `H(0) = H(1) = 1`). A non-negative exact integer order is done exactly in
GMP — a running `mpz_t` accumulating `mpz_ui_pow_ui(k, k)` — capped at
`HYPERFACTORIAL_MAX_N = 20000` to protect memory. Negative-integer, exact
non-integer, and symbolic orders are left unevaluated (`NULL`). For an inexact
argument (under `N[...]`) it falls back to the Barnes-G analytic continuation
`H(z) = Gamma[z+1]^z / BarnesG[z+1]`, built as an `Expr` tree and evaluated so it
inherits `Gamma`/`BarnesG`'s numeric kernels (machine and MPFR); a result that
overflows a `double` is promoted to an extended-exponent MPFR real or made
`Overflow[]` by `numeric_promote_result_if_degenerate`. `Product` relies on this
head to recognise `prod k^k`.

**Data structures.** `Expr`; GMP `mpz_t` for the exact integer product; the
continuation is an `Expr` over `Gamma`/`BarnesG` evaluated through the ordinary
numeric stack (libm / MPFR). No ND kernel and not on `pack.c`'s `AWARE` list;
`Compile[]` does not lower it (`CompileDiagnostics` reports `Compiled -> False`).
Attributes: `Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** The exact path is `O(n)` growing-bignum multiplies,
hard-capped at `n = 20000`; the continuation costs one `BarnesG`/`Gamma`
evaluation. Exact non-integer and negative-integer orders stay symbolic (the
`K`-function continuation is reached only numerically, via Barnes G). No
closed-form derivative rule.
