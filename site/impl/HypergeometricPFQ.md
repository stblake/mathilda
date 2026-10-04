---
source: src/special_functions/hypergeopfq.c
references:
  - "DLMF §16 — generalized hypergeometric functions pFq."
---
**Algorithm.** `builtin_hypergeometric_pfq` handles `pFq(a;b;z)` with the series
`Sum_k (prod_i (a_i)_k / prod_j (b_j)_k) z^k / k!`. The pipeline is: validate
(arity 3, first two args `List`); a packed/`NDArray` `z` buffer maps the real
`double` series `sf_machine_pfq` over the elements **only** when every parameter
is an exact `Real` (an exact-integer or complex buffer, or a parameter that is
not a bare `Real`, degrades via `ndarray_delist_and_reeval`, since exact
parameters hit the exact/closed-form steps that follow); thread over a `List`
`z`; **generic cancellation** of a parameter common to the upper and lower lists
(kept if it is a non-positive integer); `z == 0 -> 1`; **termination** to the
explicit degree-`n` polynomial when an upper parameter is a non-positive integer
`-n`; **reduction** to elementary closed forms (`0F0 = Exp`, `1F0 = (1-z)^-a`,
`0F1` cosh/sinh, `1F1(1;2;z)`, several `2F1`, and a table-backed very-well-poised
`4F3` Ramanujan `1/Pi` value); then **numeric** direct series summation. Numeric
uses a machine `double complex` sum, or an MPFR `(re,im)` sum when the governing
(minimum-contagion) precision exceeds 53 bits; a badly-cancelling machine sum
whose peak term dwarfs the answer is re-summed in MPFR at `53 + lost_bits` and
demoted back, so argument precision — not series conditioning — fixes the output
type. Convergence gate: `p <= q` (entire), `p == q+1` (`|z| < 1`); otherwise
(and for `|z| >= 1` at `p = q+1`) the call stays unevaluated — analytic
continuation is deliberately deferred.

**Data structures.** `Expr`; machine `double complex` path; an MPFR `cpx_t`
(`mpfr_t` re/im pair, using `mpfr_complex_div`) path; `Pochhammer`/`Factorial`
`Expr` builders for the terminating polynomial. The shared machine kernel
`sf_machine_pfq` (`src/special_functions/sf_machine.c`) backs both the ND
element map and the `Compile[]` VM. ND: variadic N-ary kernel
`NDKN_HypergeometricPFQ`; `pack.c` explicitly lists all four hypergeometric heads
on `AWARE` (the convenience heads rewrite to `PFQ`, so the gate must not
materialise at the rewrite). Attributes: `NumericFunction`, `Protected` — **not**
`Listable` (threading over a `List` `z` is done by hand in the builtin).

**Complexity / limits.** Term caps `HGPFQ_MACHINE_MAX_TERMS = 200000`,
`HGPFQ_MPFR_MAX_TERMS = 1000000`. Numeric only in the convergent regime.
`Compile[]` lowers the head at scalar shape (`Compiled -> True`) but **not** at
rank-1 array shape (`Compiled -> False`), where the registered ND kernel still
serves the packed / visible-`NDArray` fast path at the REPL.
