# Integrate`RischMacsyma — Phase B, fourth increment: genuine one-extension recursion

Goal: replace the flat "SolveAlways over all tower vars at once" ansatz with the
genuine **one-extension-at-a-time recursion** (Bronstein/Maxima risch.lisp): peel
the top kernel, integrate the polynomial/Laurent part by recursing into the lower
field for each coefficient, verify by construction. Closes what the flat towers
cannot: **mixed exp/log towers** and **rational lower-field coefficients**.

Target closures (currently decline, verified elementary):
- T1  ∫ (E^x/x + E^x Log[x]) dx = E^x Log[x]         (independent mixed)
- T2' ∫ (Log[1+E^x] + x E^x/(1+E^x)) dx = x Log[1+E^x] (nested log-over-exp)
- T3  ∫ (1/(x^2 Log[x]) - Log[Log[x]]/x^2) dx = Log[Log[x]]/x (rational coeff)

## Design (validated in REPL)
- `RmTower`: ordered kernels (kind LOG/EXP, kernel, arg, t-symbol, Dcoef), + subrules.
- `rm_tower_build`: collect logs+exps, order innermost-first (containment; tiebreak
  EXP-deeper so primitive recursion sits on top / RDEs bottom out in C(x)),
  structure-theorem check (each Dcoef in K_{i-1}: triangular + no foreign kernel).
- `rm_field_integrate(F, T, L, x)`: recursive.
  - L<0 → BronsteinRational (base field C(x)).
  - LOG top → `rm_int_primitive_poly`: q_i' + (i+1)q_{i+1}Dt = p_i, each solve is
    `rm_limited_field_integrate` at level L-1 (recursion) with new-log fold-back.
  - EXP top → `rm_int_hyperexp_poly`: Laurent; i=0 recurse L-1, i≠0 `rm_field_rde`.
  - proper rational part (nonzero remainder) → DECLINE this increment (tower Hermite/
    Rothstein-Trager later); all three targets have zero proper part at every level.
- `rm_field_rde`: base (Dcoef,p in C(x)) → rm_solve_rde; else NULL (general field
  RDE later).
- Wrapper `rm_recursive_tower_case`: build, substitute, whole-tower rational gate,
  integrate, back-substitute t→kernels, **diff-back verify** (search, not decision).
- Wire after rm_exp_tower_case, before rm_trig_frontend.

## Steps
- [ ] RmTower + rm_tower_build + rm_tower_free (structure theorem)
- [ ] rm_field_integrate + rm_int_primitive_poly + rm_limited_field_integrate
- [ ] rm_int_hyperexp_poly + rm_field_rde
- [ ] rm_recursive_tower_case wrapper + wire in
- [ ] tests (T1/T2'/T3 diff-back=0 + decline regressions) + valgrind
- [ ] docs: RISCH_STATUS §3.12/§6, calculus.md, changelog; memory
- [ ] commit + push to main

## Non-goals (this increment)
- Nonzero proper rational part at any level (tower Hermite / Rothstein-Trager).
- General field RDE (mixed exp-top towers with lower-field structure).
- Algebraically dependent kernels the evaluator merged (E^x E^(E^x)→E^(x+E^x)).
