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

## M3 — QE Phase-6e augment-retry  (partial D3)  → v0.321 ✅ DONE
- [x] extracted `cad_build_augmented`; rewired reduce_cad_nvar + reduce_cad_qe
- [x] no regression (410 tests, corpus 174/174, parametric QE works); leak-neutral (A/B valgrind)
- [NOTE] #28/#36 NOT unblocked — walled by qqbar degree cap / projection, not nullification.
        That wall is the D5 / further-D3 frontier (SOS certificate).

## D5 — SOS / Positivstellensatz emptiness certificate  → v0.322 ✅ DONE
(Built as ONE unified SOS/Putinar engine — the SDP tier covers the strictly-positive
polytope class too, so a separate Handelman-LP tier was unnecessary for the targets.)
- [x] reduce_sos.{c,h}; refutation reduction (q<0 target, K={g>=0,h==0}); linear-eq elimination
- [x] Putinar SOS assembly over exponent/mpq bookkeeping (GBPoly)
- [x] numeric SDP (alternating projection, LAPACK dsyev/dgesv) + eigenvalue-floor interior
- [x] numeric zero-find (penalty descent) → exact rational verify → facial reduction
- [x] exact rational null-space rounding + exact LDL^T PSD verify (sound-or-decline)
- [x] hook at reduce.c BEFORE cad; #34/#40 → False; strictly-pos → False
- [x] SOUNDNESS: non-empty regions never False (sweep); Minimize #40 picks it up -> {0,centroid}
- [x] leak-clean (success+decline == baseline); USE_LAPACK=0 degrades (clean stub); c99/messages gates
- [x] tests (test_sos_emptiness, updated M2 preempt); docs + changelog; bump v0.322

## Close-out
- [ ] update tasks/reduce_deficiencies.md (what closed / what remains)
- [ ] rebuild code-review graph; tasks/lessons.md if any correction

## Review
_(to be filled as phases complete)_
