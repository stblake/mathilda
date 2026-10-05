# Task: General algorithmic improvements to definite integration (stress-test campaign)

Plan: `/Users/user/.claude/plans/here-are-some-stress-indexed-nygaard.md`
(full phased roadmap: Phase 1 reliability/speed, Phase 2 new transform methods,
Phase 3 deep research). 20 stress integrals → 7 general gaps.

## Phase 1 — reliability & speed (v0.291) — COMPLETE, verified, awaiting commit

### 1a. Confidence-first cascade reorder (`src/calculus/integrate.c`)
- [x] Move `integrate_beta_try` / `integrate_trigpower_try` / `integrate_ramanujan_try`
      to run BEFORE `integrate_newton_leibniz_try_pv` (after residue/symmetry).
- [x] Verified: #5 54s→0.007s, #8 19.5s→0.0015s, #20 54s→0.005s; spurious
      `Integrate::nonelem` gone; integer-power/elementary cases still owned by FTC.

### 1b. Rational half-line Mellin accepts the true strip (`integrate_residue.c`)
- [x] `g_coupled` flag decouples "coupled FindInstance mode" from `g_nbounds`; record
      the sound `absorb_fact` bounds even in the coupled path (convergence gate on a
      directly-bounded exponent now discharges).
- [x] Integer-`s` guard reads literal `s` (`res_reim_direct`) so a symbolic exponent
      survives a representative point that lands on an integer.
- [x] Verified: #4 `x^a/(x²+2bx+c)`, `-1<a<1 && c>b² && b>0` → closes 0.42s (was 10s
      decline); numeric oracle match.

### 1c. Sound assumption-aware residue cleanup (`integrate_residue.c`, `refine.c`)
- [x] `res_refine_final` (top of `integrate_residue_try`, decoupled from `g_all_pos`):
      `Refine` → `Simplify` → guarded `FullSimplify` (radical/trig-free, node-capped),
      keep-if-strictly-smaller. NO ComplexExpand (explodes #11 → hangs Simplify).
- [x] `refine.c deep_positivity_walk`: `Arg[neg-real compound] → Pi` (sound, via
      `deep_sign`) — the gap that blocked #13.
- [x] Verified: #13 `Log[x]/(x²-a²)` PV → `π²/(4a)` (was ~200-node Arg soup);
      #11/#16 correct & unchanged (no hang); #17 clean Sec form.

### Verification (all green)
- [x] `make -j` clean (only pre-existing -Wmisleading-indentation warnings, not my code).
- [x] Suites: integrate_residue / ramanujan / symmetry / beta / newton_leibniz /
      principalvalue / dispatch / integrals all pass. (ramanujan's 2 soft-FAILs are the
      pre-existing Hypergeometric1F1/2F1 Mellin cases — `rec_pfq` gap, unrelated.)
- [x] `make check-c99`, `make check-messages` clean.
- [x] All closed forms cross-checked against an independent numeric oracle.
- [x] version.h → 0.291; changelog `docs/spec/changelog/2026-10-05.md`; docs
      `docs/spec/builtins/calculus.md`.

### Not committed — awaiting user go-ahead (commit + tag v0.291).
Suggested single commit "feat(integrate): Phase 1 — cascade reorder + rational-Mellin
strip + sound residue cleanup ; v0.291" (or split 1a/1b/1c with sequential bumps).

## Phase 2 — integral-representation recognizers + Euler→Beta·2F1 (v0.292) — COMPLETE, awaiting commit

### 2a. `src/calculus/integrate_intrep.c` (new) — half-line integral reps
- [x] Laplace-Bessel `E^(-p x)BesselJ[nu,q x]` → #2 `1/Sqrt[a^2+c^2]` (general nu).
- [x] Bessel-K cosh `E^(-A Cosh[x])Cosh[n x]` → #15 `BesselK[n,a]` (+ n=0 bare → K_0).
- [x] Bessel-K exp `x^(nu-1)E^(-A x-B/x)` → #19 `Sqrt[Pi/a]e^(-2Sqrt[ab])` (TrigToExp LAST, no
      trailing Simplify — it re-folds exp→cosh-sinh). General nu → BesselK.
- [x] Airy `Cos[p x^3+q x]` → #12 `Pi AiryAi[a]`.
- [x] Wired: METHOD_INTEGRAL_REP, cascade after Ramanujan/before NL, init, CMake mathilda_common.

### 2b. Euler→Beta·2F1 (`integrate_euler_2f1_try`, `integrate_beta.c`)
- [x] `x^(a-1)(1-x)^(b-1)(alpha+beta x)^e` on [0,1] → `alpha^e Beta[a,b] 2F1[-e,a,a+b,-beta/alpha]`.
      #18 → `Beta[a,b](1+1/c)^(-a)c^(-(a+b))`. Gate via `prove_ref` (Refine, not Simplify).
- [x] Stale beta-suite control updated (1/(1+x) Method->Beta now = Log[2], a valid Euler close).

### Verification
- [x] new `tests/test_integrate_intrep.c` (self-certifying Simplify[res-ref]===0) — PASS.
- [x] all definite suites green; check-c99 / check-messages clean.
- [x] version.h → 0.292; changelog + calculus.md docs updated.

### Awaiting user go-ahead to commit + tag v0.292.

## Phase 3 (next) — higher-order hyperbolic residue (#10), symmetry-aware Feynman (#9,#14), Malmsten (#3)

Still unevaluated (unbroken): #3 #6 #9 #10 #14. (#6 Hurwitz-zeta Mellin is a candidate too.)
#10 still ~55s to decline (unchanged baseline; Phase 3 turns it into a fast close).
