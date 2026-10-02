---
references:
  - "B. C. Carlson, *Computing elliptic integrals by duplication*, Numer. Math. **33** (1979) 1-16."
  - "B. C. Carlson, *Numerical computation of real or complex elliptic integrals*, Numer. Algorithms **10** (1995) 13-26 — `R_J` and the `p < 0` transformation."
  - "W. H. Press et al., *Numerical Recipes in C*, 2nd ed. (Cambridge, 1992), §6.11 — the `rj`/`rc` loops."
  - "DLMF §19.25.14 — the third kind in terms of `R_F` and `R_J`."
source: src/special_functions/elliptic.c
---
**Algorithm.** `builtin_ellipticpi` dispatches on argument count (2 complete, 3 incomplete;
anything else emits `EllipticPi::argt`). Both forms run the numeric path first, then the
reductions: `Π[0, m] = K(m)`, `Π[n, 0] = π/(2√(1−n))`, `Π[1, m] = ComplexInfinity`,
`Π[n, ∞] = Π[∞, m] = 0`, and for the incomplete form `Π[n, 0, m] = 0`,
`Π[0, φ, m] = F(φ, m)`, `Π[n, π/2, m] = Π[n, m]`,
`Π[n, φ, 0] = ArcTanh[√(n−1) tan φ]/√(n−1)`, plus the oddness fold in the amplitude.

The machine kernels are Carlson compositions:

- complete — `Π(n|m) = R_F(0, 1−m, 1) + (n/3) R_J(0, 1−m, 1, 1−n)`, checked against mpmath
  at ≤ 2 ulp over `n, m ∈ (0,1)` and ≤ 3.3 ulp with `n ∈ (−6, 1)`;
- incomplete — `Π(n; φ|m) = s R_F(c², 1−m s², 1) + (n/3) s³ R_J(c², 1−m s², 1, 1−n s²)`
  with `s = sin r`, `c = cos r` on the principal strip, ≤ 3.8 ulp over 400 random
  `(n, φ, m)`, extended by the quasi-period `Π(n; φ+kπ|m) = Π(n; φ|m) + 2k Π(n|m)`
  (verified exactly against mpmath at `k = −1, 1, 2`). `c²` is computed as `c*c` for the
  same reason as in `EllipticF`.

`R_J` needs `R_C` (the degenerate `R_F`), which is why both landed together.

**`n > 1` is a decline, not a principal value.** The standing reason for having no kernel
here was that `n > 1` needs a Cauchy principal value, and that reason was simply wrong:
the value past the pole at `sin²t = 1/n` is genuinely complex —
`Π(3/2 | 1/2) = −0.456720313453 − 2.72069904635 i`, mpmath and Arb agreeing. So what a
`double` kernel owes there is a decline (`carlson_rj` refuses `p ≤ 0`), exactly as `K`, `E`
and `F` decline outside their real domains, and `flint_num_elliptic_pi` answers through
Arb's `acb_elliptic_pi`.

**The pole belongs to the complete form only.** `Π[1, m]` diverges because the integrand
carries a `1/cos²t` and the upper limit is `π/2` — verified divergent as `1/√ε`: 21.5, 221,
2221, 22214 at `ε = 10⁻², 10⁻⁴, 10⁻⁶, 10⁻⁸`. The incomplete form has no such pole, and
`Π[1, 1, 1/2]` is `1.73199154202`. The test sits ahead of the numeric path, or Arb's
non-finite ball would come back as an unevaluated `EllipticPi[1., 1/2]`.

**Data structures.** `double` scalars through `carlson_rj`/`carlson_rc`; `acb_t` above
machine precision. Registration is split by arity: a binary (`REG_B`) kernel for the
complete form and an n-ary (`REG_N`) one for the incomplete. The element-wise NDArray layer
tops out at arity 2, so only the complete form rides that buffer, while `Compile[]` lowers
both — the three-argument shape through the generic n-ary `OP_KERNN` path, which reserves a
consecutive register block and reads the kernel from the symbol's `ndarray_nary_kernel`.

**A sequencing trap, recorded.** `packed_aware` is a property of the *symbol*, not of one
arity, so registering the two-argument kernel also stopped the transparency gate
materialising packed `List`s for the three-argument form. The `ndarray_delist_and_reeval`
fallback in all three argument positions is what keeps that correct, and it had to land
first — a visible `NDArray` left unevaluated is a wrong answer, not a slow one.

**Complexity / limits.** `O(1)` per element: 41 ns/element for the complete form at 2 ulp
(a 970× speedup over the Arb-per-element path it replaced, and 7.7× faster than the
equivalent SciPy composition), 1 ulp for the incomplete one. Real principal domain only
(`m < 1`, `n < 1`, and `1 − n sin²φ > 0`). The `n`- and `m`-derivatives of the *incomplete*
form stay inert `Derivative[…]`: a least-squares fit over the natural candidate basis does
not recover them (residual 1.5 relative, no simple rational coefficients), so the closed
form involves terms that basis does not span — and an inert derivative is honest where a
guessed one corrupts every caller silently.
