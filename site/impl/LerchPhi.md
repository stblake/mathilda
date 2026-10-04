---
source: src/special_functions/lerchphi.c
references:
  - "DLMF §25.14 — the Lerch transcendent."
  - "A. Erdelyi et al., Higher Transcendental Functions, Vol. I, §1.11(8) — the large-z continuation."
---
**Algorithm.** `builtin_lerchphi` evaluates `LerchPhi[z, s, a] = Sum_{k>=0}
z^k/(k+a)^s`, the common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`.
Exact reductions: `z = 0 -> a^-s`; `s = 0 -> 1/(1-z)`; `z = 1 -> Zeta[s, a]`;
`z = -1 -> 2^-s (Zeta[s, a/2] - Zeta[s, (a+1)/2])` (with special handling of
`s = 1` via a digamma difference, and of half-integer `a` via a Dirichlet-beta
closed form); `a` a positive integer `m -> z^-m (PolyLog[s, z] - Sum_{j<m}
z^j j^-s)`; `a` a non-positive integer reduces onto `PolyLog`; `s` a negative
integer `-n -> (z d/dz + a)^n [1/(1-z)]`, a rational function built with the
Euler operator and `Together`. Options `IncludeSingularTerm -> True` (at a
non-positive-integer `a` gives `ComplexInfinity`) and `DoublyInfinite -> True`
(`Phi(z,s,a) + z^-1 Phi(1/z,s,1-a)`). Numeric (at least one inexact operand):
`|z| < 1` (or `|z| = 1` with `Re s > 1`) sums the complex-MPFR power series
(symmetric power `((k+a)^2)^(-s/2)` for `Re a < 0`); `|z| > 1` off the cut with
`|Log z| < 2 pi` uses the Erdelyi large-`z` continuation through `HurwitzZeta`
and `Gamma`, else stays symbolic.

**Data structures.** `Expr`; `lcx` (`mpfr_t` re/im) MPFR-complex toolkit; the
exact and continuation paths reuse the `PolyLog`/`Zeta`/`HurwitzZeta`/`Gamma`
builtins. ND: N-ary kernel `NDKN_LerchPhi` (arity 3), element-wise over the `z`
buffer; `packed_aware` from kernel registration. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Series term cap `wp*64 + 100000`; the large-`z`
continuation truncates at `NMAX = 2000` with an optimal-truncation safeguard.
`|z| > 1` is generally left symbolic outside the Erdelyi domain (integer `s` on
the cut is not implemented). `Compile[]` lowers at scalar shape
(`Compiled -> True`) but **not** at rank-1 array shape (`Compiled -> False`); the
ND kernel still serves the packed / visible-`NDArray` fast path at the REPL.
