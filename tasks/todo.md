# Reduce improvement campaign (per tasks/reduce_deficiencies.md)

Plan: `/Users/user/.claude/plans/let-s-improve-our-implementation-squishy-tide.md`
Full campaign: M1 + M2 + M3 + D5a + D5b. Land in order; each a tagged commit.

## M1 — fix 1-arg/2-arg `Reduce[expr, domain]` misparse  (D1+D2)  → v0.319 ✅ DONE
- [x] `reduce_is_domain_symbol` helper (Reals/Complexes/Integers/Rationals/Booleans/Primes)
- [x] `reduce_infer_var_list` free-var collector (exclude domains/constants; alphabetical)
- [x] arg-bind rewrite + relax 1-arg guards; `owned_vars` freed at every `owned_list` site
- [x] REPL proof + regression pins; extended tests/test_reduce.c (409 pass, corpus 174/174)
- [x] docs + changelog + bump v0.319; leak-clean (420 blocks == startup baseline)

## M2 — TimeConstrained-preemptible CAD, leak-free  (D6)  → v0.320 ✅ DONE
- [x] `tc_deadline_passed()` non-jumping query (core.h/core.c)
- [x] `CAD_POLL()` at CAD loop heads → existing clean-return; reduce_cad.c #include core.h
- [x] REPL proof (D1 grind now returns $Aborted at 2.5s; normal CAD unchanged); valgrind clean (420 blocks)
- [x] fixed latent ifun-suppress leak across TC abort (tc_run_guarded saves/restores ifun depth)
- [x] test_timeconstrained_preempt added (410 pass); corpus 174/174; groebner/core tests pass; bump v0.320

## M3 — QE Phase-6e augment-retry  (partial D3)  → v0.321
- [ ] extract `cad_build_augmented`; rewire reduce_cad_nvar + reduce_cad_qe
- [ ] #28/#36 resolve; QE corpus unchanged; valgrind; tests; bump/tag

## D5a — Handelman exact-LP emptiness certificate  → v0.322
- [ ] reduce_sos.{c,h}; refutation reduction; GBPoly + ideal reduction
- [ ] exact rational `{Ax=b,x>=0}` LP oracle (from reduce_fm core)
- [ ] hook at reduce.c:737; tests; valgrind; bump/tag

## D5b — SOS/Putinar numeric-guided + exact-verified  → v0.323
- [ ] affine-subspace build; alternating-projection SDP guide (LAPACK)
- [ ] rational rounding + exact correction + rational LDL^T PSD verify
- [ ] #34/#40 → False; non-empty region not wrongly False; USE_LAPACK=0 degrades; bump/tag

## Close-out
- [ ] update tasks/reduce_deficiencies.md (what closed / what remains)
- [ ] rebuild code-review graph; tasks/lessons.md if any correction

## Review
_(to be filled as phases complete)_
