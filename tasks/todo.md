# Task: Implement LerchPhi[z, s, a]

Lerch transcendent Phi(z,s,a) = Sum_{k>=0} z^k / (k+a)^s. Generalizes Zeta,
HurwitzZeta and PolyLog. Template: src/special_functions/{zeta,hurwitzzeta,polylog}.c.

## Plan

- [ ] `src/special_functions/lerchphi.{c,h}` — builtin_lerchphi + lerchphi_init
  - Option parsing: 3 positional args + DoublyInfinite / IncludeSingularTerm;
    argrx for <3 args, nonopt for bad options beyond position 3.
  - Exact reductions (built symbolically, eval_and_free): z==0 -> a^-s; s==0 ->
    1/(1-z); z==1 -> Zeta[s,a]; z==-1 -> 2^-s(Zeta[s,a/2]-Zeta[s,(a+1)/2]);
    a positive integer m -> z^-m(PolyLog[s,z]-Sum_{j=1}^{m-1} z^j j^-s);
    s negative integer -n -> (z d/dz + a)^n[1/(1-z)] Together'd;
    IncludeSingularTerm->True & a non-positive integer -> ComplexInfinity.
  - DoublyInfinite->True -> LerchPhi[z,s,a] + z^-1 LerchPhi[1/z,s,1-a].
  - Numeric single-sum kernel (complex MPFR `lcx`), |z|<1, machine + MPFR, complex.
- [ ] sym_names.{c,h}: SYM_LerchPhi, SYM_DoublyInfinite, SYM_IncludeSingularTerm
- [ ] core.c: include + lerchphi_init()
- [ ] info.c: docstring
- [ ] calculus/deriv.c: d/dz, d/da elementary rules; d/ds generic Derivative
- [ ] tests/CMakeLists.txt + tests/test_lerchphi.c
- [ ] docs/spec/builtins/special-functions.md + changelog/2026-06-15.md
- [ ] Build, run tests, valgrind clean.

## Scope note
Numeric continuation for |z|>1 (general a) is NOT implemented (stays symbolic).

## Review (done)
- All planned pieces implemented and wired. `lerchphi_tests` (16 cases) pass;
  zeta/hurwitzzeta/polylog/deriv neighbours still green. Main binary builds clean
  under -std=c99 -Wall -Wextra. Valgrind: no LerchPhi/src frames in any leak
  stack (only the documented ~12,800 B macOS dyld/Accelerate baseline).
- Bug found & fixed: negative-integer-s operator collapsed when z was a concrete
  number (D of a constant); now differentiates against a fresh placeholder symbol
  then substitutes z (LerchPhi[2,-1,a] -> 2-a, was -a).
- Numerics cross-checked independently: z·LerchPhi[z,s,1] == PolyLog[s,z] and
  LerchPhi[1,s,a] == Zeta[s,a] agree across the distinct kernels.
