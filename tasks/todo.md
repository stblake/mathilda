# P2 — Real rational-trigonometric integration in the GENUINE Risch engine

Roadmap P2 (Bronstein §5.10 hypertangent + the complex-exponential tower).

## DIRECTIVE (2026-07-14, user)
**Do NOT patch Risch deficiencies with the Weierstrass (t = Tan[x/2]) substitution.**
Rational trigonometric integrands must be integrated by the genuine Risch machinery —
the complex-exponential tower (TrigToExp → integrate over the E^(ix) monomial →
ExpToTrig, with REAL reconstruction of the log part) and the real §5.10 hypertangent
case — producing real, clean closed forms. The half-angle substitution and the
Jeffrey–Rich continuous-Weierstrass integrator are NOT acceptable as the Risch answer.

An earlier attempt (rt_weierstrass_case dispatch + a top-level HypertangentTrig
promotion) was REVERTED in full (branch reset to main) as a Weierstrass patch.
See tasks/lessons.md + memory `project_risch_no_weierstrass_patch`.

## Baseline on clean main (genuine routes, verified 2026-07-14)
GOOD (genuine, clean): Sin[x]^2 -> (2x-Sin[2x])/4, Sin[x]^3 -> (Cos[3x]-9Cos[x])/12,
  Tan[x] -> -Log[Cos[x]], Cot[x] -> Log[Sin[x]], Tanh/Coth clean.
DEFICIENT (the real work):
  - I-laden: Csc[x] -> Log[..+I Sin..]; Tan[x]/(3+Tan[x]^2) -> (2I)x - Log[..I..].
  - Declined: Sec[x], Sec[x]^3, 1/(2+Cos[x]).

## STATUS (2026-07-14): Fix 1 ✅, Fix 2 ✅, Fix 3 ✅ (all committed on feature/risch-p2-weierstrass)
Fix 1 (3d74197): §5.10 quadratic residues -> real ArcTan.
Fix 2+3: rt_exp_ratreduce_case wired into rt_trig_frontend (closes Sec/Csc/1-over-(2+Cos));
  cx_reim + rt_realify real reconstruction of the exp-tower log part (numeric diff-back gate);
  2-arg ArcTan[u,v] derivative in deriv.c. All rational trig now REAL (verified sweep).
REMAINING completeness item: 1/(a+b Sin[x]) still declines (rt_exp_ratreduce_case sin-denominator).

## Fix 1 — §5.10 gate: admit irreducible-quadratic residues (genuine, low-risk)
The direct hypertangent case (rt_hypertangent_case, Tan/Cot/…) routes through
rt_hypertan_family, whose normal-pole pre-gate requires TOTAL denominator degree
<= 1. The §5.10 residue driver ALREADY realises irreducible-quadratic residues as
real ArcTan (verified: driver on 1/(3+t^2) -> ArcTan[t/Sqrt3]/Sqrt3). Relax the gate
to admit a normal part whose every IRREDUCIBLE factor over Q has degree <= 2
(new helper rt_max_irr_degree via Factor). This is the genuine §5.10 driver, NOT
Weierstrass. Closes Tan[x]/(3+Tan[x]^2) -> real ArcTan.
- [ ] rt_max_irr_degree + gate change in rt_hypertan_family.
- [ ] Test: assert_tan_real on irreducible-quadratic-residue Tan integrands.

## Fix 2 — Real reconstruction of the complex-exponential log part
When the E^(ix) rational integral's log part comes back as conjugate complex logs
Log[E^(ix) - alpha] (alpha complex), combine each conjugate pair into a real
Log/ArcTan (the LogToReal step, as intrat's rt_frac_lrt already does for the
rational base field). Applies to Csc and any rational-in-E^(ix) with a real
integrand. Investigate: rt_exp_ratreduce_case / rt_frac_case / the exp frac log
part — where the complex logs are emitted and why real reconstruction is skipped
on the exponential tower.
- [ ] Real-reconstruct the exp-tower log part; Csc[x] -> real (Log[Tan[x/2]] or
      Log[1-Cos[x]]-Log[Sin[x]] class), no I.

## Fix 3 — Completeness of the rational-function-of-E^(ix) route
Sec, Sec^3, 1/(2+Cos[x]) are rational functions of E^(ix) that the route currently
DECLINES. Diagnose why rt_exp_ratreduce_case (kernelize -> rational integral in t=E^(ix)
-> back-substitute) declines them (F free of x, u'=i free of x should qualify) and
extend so they close, with real reconstruction (Fix 2). Expected real forms:
  Sec[x] -> Log-of-real / ArcTanh[Sin[x]] class; 1/(2+Cos[x]) -> real ArcTan.
- [ ] Close Sec/Sec^3/1/(2+Cos) via the genuine exp route, real output.

## NON-NEGOTIABLE constraints
- No t = Tan[x/2] Weierstrass substitution anywhere in the Risch path.
- No routing rational trig through Jeffrey-Rich as the Risch answer.
- Every closed form correct by construction / diff-back gated; real (Complex-free)
  for real integrands. Clean multiple-angle / ArcTan / Log forms, not raw kernels.

## Verification
Clean -std=c99 -Wall -Wextra; scoped *_tests green (grep FAIL:); REPL pins via JSON
pipe (python parser, sed mangles Floor/`/`); valgrind vs Sin[1.0] baseline; diff-back
Complex-free. Commit per genuine fix.
