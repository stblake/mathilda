# DiffUnderInt Improvement Campaign — todo

Plan: `~/.claude/plans/similar-to-our-contour-nested-patterson.md`
Target file: `src/calculus/integrate_diffunderint.c`. Start version: 0.284.

## Phase 1 — Correctness: hang + wrong answer  [DONE v0.285]
- [x] 1a. Reproduced hang: `integrate_over_param` → engine on `Log[1/(a²+4b²)^(1/4)]`
- [x] 1a. `integrate_over_param` now PowerExpands + gates (radical/trig/gaussian-of-p → decline)
- [x] 1a. Regression test added; squared-sinc closes for b>0, declines clean for Reals
- [x] 1b. Audited: `Sin·Sin` was non-real for b<c, Indeterminate for b=c (verify can't catch base)
- [x] 1b. `has_trig_product_of_x`+TrigReduce in inner_definite; zero-base-first pass in stage_quadrature
- [x] 1. tests (test_phase1_robustness) + docs + changelog + version 0.285

## Phase 2 — Stage B: first-order linear ODE in parameter (λ≠0)
- [ ] `collect_atoms`, `boundary_value`, `detect_linear_ode` (SolveAlways)
- [ ] `solve_linear_ode` (homogeneous branch first)
- [ ] `stage_linear_ode` + hook in `stage_quadrature`
- [ ] Closes `Exp[-a^2 x^2]Cos[bx]` = (√π/2a)Exp[-b^2/4a^2]; verify; tests/docs/version

## Phase 3 — Self-similar substitution recognizer
- [ ] `detect_selfsimilar_ode` for `Exp[-a^2 x^2 - b^2/x^2]` = (√π/2a)Exp[-2ab]; tests/docs/version

## Phase 4 — Finite-period rational-trig inner family (residue reuse)
- [ ] `is_rational_trig_in_x`, `rational_trig_inner` (symmetry folds + gates) via `integrate_residue_try`
- [ ] Hook into `inner_definite`; closes the two Log-trig; tests/docs/version

## Phase 5 — Reverse recognition (Γ′)
- [ ] `strip_log_x_power`, `find_log_derivative_param`, `known_integral_closed_form`, `stage_reverse_feynman`
- [ ] Closes `Exp[-x]x^(a-1)Log[x]` = Γ[a]PolyGamma[0,a]; tests/docs/version

## Phase 6 — Repeated differentiation / power reduction
- [ ] `inv_x_power`, `multi_base`, `stage_repeated_quadrature`
- [ ] Closes `Sin[ax]^3/x^3` = 3πa^2/8 via DiffUnderInt; tests/docs/version

## Phase 7 — Output-quality finalize (assumption-aware)
- [ ] Extend `diui_finalize`: Sqrt collapse, ArcTan reciprocal fold, Log contraction, ParamBound→assumptions
- [ ] Cleans In[2]/In[3]; tests/docs/version

## Review
- (to be filled in as phases land)
