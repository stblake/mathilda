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

## Phase 2 (next) — integral-representation recognizers (#2,#12,#15,#19) + Euler→Beta·2F1 (#18)
## Phase 3 — higher-order hyperbolic residue (#10), symmetry-aware Feynman (#9,#14), Malmsten (#3)

Still unevaluated (unbroken, future phases): #2 #3 #6 #9 #10 #12 #14 #15 #18 #19.
#10/#15 still ~55s to decline (unchanged from baseline; narrow FTC pre-decline deferred —
Phase 2/3 turn them into fast closes).
