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

## Review (campaign complete — v0.319–v0.322)

**Reframe:** D1+D2 were not what the seed doc assumed (fast-give-up vs route-through-QE)
— they were ONE correctness bug: the 2-arg/1-arg `Reduce[expr, domain]` misparsed the
domain symbol as a variable. Fixing the root cause (M1) subsumed both and was far simpler.

**Delivered (each its own tagged commit):**
- **v0.319 (M1)** — `Reduce[expr]` / `Reduce[expr, dom]` infer vars + detect domain → run
  the correct 3-arg engine. Broad correctness win; closes D1+D2.
- **v0.320 (M2)** — CAD preemptible by `TimeConstrained` (non-jumping `tc_deadline_passed`
  + clean-return polls); closes D6. Also fixed a latent ifun-suppress leak across TC abort.
- **v0.321 (M3)** — QE shares the Phase-6e augment-retry with plain Reduce; partial D3.
- **v0.322 (D5)** — new `reduce_sos.c`: exact SOS/Positivstellensatz emptiness certificate
  (numeric SDP + exact rational verify + facial reduction for interior zeros). #34/#40 →
  False; Minimize picks them up.

**Verification:** reduce_tests 414 pass, corpus 174/174, refine/minimize/solve/groebner/core
pass. Leak-clean on all new paths (fixed-init 420-block baseline). Soundness swept (SOS
never False on a non-empty region). check-messages / check-c99 clean. USE_LAPACK=0 degrades.

**Open frontier** (documented in reduce_deficiencies.md): full-D3 universal QE over a
continuum (#28/#36, qqbar-degree-cap wall), D4 Max/Abs under a quantifier, non-polytope SOS.

**Not pushed** — commits + tags v0.319..v0.322 are LOCAL (push deferred per the
outward-facing-action guard; `git push --follow-tags` when ready).
