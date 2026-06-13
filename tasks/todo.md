# Fix: ExpIntegralEi MPFR precision overflow on large |z|

## Bug
`N[ExpIntegralEi[-10^60 I] + I Pi, 20]` aborts:
  init2.c:53: MPFR assertion failed (precision out of range)
Root cause: ei_mpfr_complex (expintegralei.c:410) and ei_mpfr_real (x<0, :377)
compute guard bits as `(long)(|z|/ln2)`; |z|~1e60 overflows long -> bogus wp ->
mpfr_init2 assertion. Convergent series is also infeasible at that magnitude.

## Plan
- [ ] Add ecx_exp (complex exp) and ecx_div (complex divide) helpers.
- [ ] Add ei_ecx_asymptotic: Ei(z) ~ (e^z/z) Sum k!/z^k + i*pi*sign(Im z),
      truncated at the smallest term.
- [ ] Add ei_real_asymptotic for the real x<0 large path.
- [ ] Route to asymptotic when |z| large enough that its min term (~e^-|z|)
      is below target precision; else keep the convergent series (unchanged).
      This bounds convergent wp -> no overflow.
- [ ] Validate asymptotic vs convergent in the overlap region (|z|~40-60).
- [ ] Verify crash cases now return finite values; run special-function tests.
- [ ] Docs: changelog + calculus.md note.

## Review (done 2026-06-13)
All plan items complete.
- Root cause: ei_mpfr_complex:410 / ei_mpfr_real:377 computed guard bits as
  `(long)(|z|/ln2)`; |z|~1e60 overflowed long -> out-of-range mpfr_init2 prec.
- Fix: asymptotic expansion Ei(z) ~ (e^z/z)Sum k!/z^k + i*pi*sign(Im z) for
  large |z|; routing threshold also bounds the convergent path's precision so
  the overflow cannot recur. Convergent path otherwise unchanged.
- New helpers: ecx_exp, ecx_div, ei_ecx_asymptotic, ei_real_asymptotic.
- Verified: no crash on the reported case + real-negative analogue; asymptotic
  vs convergent agree on the overlap; correct branch signs all quadrants;
  machine + arbitrary precision both correct; 0 leaks; ei/li suites pass.
- Docs: special-functions.md + changelog/2026-06-08.md updated.
