# DiffUnderInt — round of improvements (Rounds 1 + 2)

## Round 1 — half-line sinc / complex-Log cluster
- [ ] R1.1 Mute spurious arithmetic warnings (In[5] flood): wrap stage_quadrature
      in arith_warnings_mute_push/pop.
- [ ] R1.2 Non-even real-rational half-line family `rational_halfline_general`
      (linear any m, irreducible quadratic m=1; Log-at-inf coeffs must cancel).
      Wire as fallback in laplace_sinc_halfline and inner_definite. Unlocks In[9].
- [ ] R1.3 `diui_finalize` output cleanup via 1-arg PowerExpand + re-verify.
      Cleans In[7].

## Round 2 — Gaussian moment family (B1)
- [ ] R2.1 `gaussian_halfline`: c x^n e^{-p x^2} {1,Cos qx,Sin qx} closed forms.
      Relax contains_gaussian_exp decline so Gaussians reach the family.
- [ ] R2.2 Gaussian parameter back-integration -> Erf (engine can't). Unlocks In[10].

## Verify
- [ ] In[5] no message flood; In[7] clean; In[9] closes; In[10] closes.
- [ ] Regression: In[3],In[6],In[8] still correct; doc examples #1,2,3,5,6,7,13,14,17,20,21,22.
- [ ] valgrind clean on a representative case.
- [ ] docs/spec + changelog updated.

## Review (2026-07-08)
DONE — both rounds landed in src/calculus/integrate_diffunderint.c:
- R1.1 mute: In[5] Sin[a x]Sin[b x]/x^2 now declines SILENTLY (0 stderr lines).
- R1.2 rational_halfline_general: In[9] Exp[-c x](1-Cos[a x])/x^2 -> a ArcTan[a/c]
  - c/2 Log[1+a^2/c^2] (0.55s). Non-even sinc via real ArcTan, no complex logs.
- R1.3 diui_finalize (PowerExpand+re-verify): In[7] now clean
  -1/2 Log[a^2+c^2] + 1/2 Log[b^2+c^2].
- R2.1/R2.2 gaussian_halfline + Erf back-integration: In[10] Exp[-x^2]Sin[a x]/x
  -> Pi/2 Erf[a/2] (0.14s). Top-level Gaussian decline relaxed (families still
  never hand Gaussians to the general engine).
Verify: all 10 batch examples correct/clean; test_integrate_diffunderint (8 incl
2 new) + integrals/dispatch/derivdivides/goursat/residue/newton_leibniz/ramanujan
/limit suites all pass. valgrind: no new leaks (residual = pre-existing
Goursat/Factor evaluate-class). Docs + changelog updated.
Still open (future rounds): In[1] trig-period (B2b), In[2] log-rational base (B6),
In[4] finite-radical (B2b/B6), In[5] piecewise VALUE Pi/2 Min[a,b] (B5).
