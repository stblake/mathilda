# Task: Clear rational denominators in multivariate Reals QE (unblock Resolve class A)

Root cause (verified): chained QE emits rational atoms (`m <= 1/del`) at intermediate
levels; every multivariate real engine declines any atom with a non-constant
denominator (`nonconst_denom`), and the Reals preprocessing has no denominator-clearing
pass (it clears radicals + selectors only). CAD handles the parameter×variable PRODUCT
fine — the RESOLVE_IMPROVEMENTS.md "CAD blow-up" hypothesis was wrong.

Fix: a fully algorithmic denominator-clearing preprocessing pass
(`p/q REL 0 -> p·q REL 0 [&& q!=0]`, NNF-walked) in the multivariate Reals path.

## Implementation
- [ ] `src/solve/reduce_realfn.c`: Preprocessing 6 — rational-denominator clearing
  - [ ] `is_var_denom_power` + `has_var_denom` (cheap structural gate)
  - [ ] `clearfrac_has_branch_cut` (skip Together on Log / inverse-trig)
  - [ ] `clearfrac_relation` (Together/Numerator/Denominator, Expand[num*den], table)
  - [ ] `clearfrac_tree` (NNF walk, parallel to rationalize_tree)
  - [ ] `reduce_stmt_has_fraction` (public gate)
  - [ ] wire into `reduce_piecewise_preprocess` loop: selector -> fraction -> radical
- [ ] `src/solve/reduce_realfn.h`: declare `reduce_stmt_has_fraction`
- [ ] `src/solve/reduce.c:543-544`: add `|| reduce_stmt_has_fraction(...)` to gate
- [ ] `src/solve/reduce_qe.c:237-238`: add `|| reduce_stmt_has_fraction(...)` to gate

## Verification
- [ ] Build clean (`make -j`)
- [ ] §2 [2],[4] -> True; real limit statements (1/x->0 at inf; 1/x^2->oo)
- [ ] base engine rational atoms decide (m<=1/del; 1/x+1/(x+1)<=m; m==1/x)
- [ ] soundness: wrong limit (1/x=5 at inf) NOT True; (3x-1)=6 stays False
- [ ] no regression: §1 passing rows; univariate rational suite; test_reduce*,
      reduce_corpus; add new corpus/unit rows
- [ ] `make check-c99`, `make check-messages`

## Bookkeeping
- [ ] version.h bump v0.284 -> v0.285, tag v0.285
- [ ] changelog + docs/spec builtins note
- [ ] correct RESOLVE_IMPROVEMENTS.md hypothesis, move fixed rows
- [ ] update memory project_resolve_param_times_var_decline

## Review
(to be filled on completion)
