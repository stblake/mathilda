# Unified Risch-structure rational normalization (Together / Cancel / Simplify)

Plan: `/Users/user/.claude/plans/it-seems-like-we-fancy-bumblebee.md`

## Key finding from exploration
The FLINT canonicaliser engines already exist and are correct; the mess is the
**routing** — 8 terminal paths, `extension_autodetect` run 3×, duplicated
radical recombination, dependent-generator bail guards, Simplify re-implementing
the routing. The redesign is a **single field-structure classifier + one
dispatch** (`rat_canon_*`) reusing the existing engines, validated by a shadow
parity harness before switching the builtins over.

Reusable engines (`src/poly/flint_bridge.c`, public in `flint_bridge.h`):
- `flint_rational_together/cancel`  — plain Q + transcendental kernels (km)
- `flint_algebraic_field_canonical` — rationalised form over algebraic tower
- `flint_algebraic_field_together`  — WL-faithful, cube+ roots, symbolic radicand
- `flint_parametric_field_normalize`— Q(t..)(√k), symbolic radicand √k
- `flint_algebraic_field_normalize` — identically-zero test over tower (fix d≥5 SIGSEGV)
- `extension_autodetect` (`qafactor.h`) — algebraic constant field (one call)

## Phase 0 — Parity harness  [DONE 2026-07-22]
- [x] Classifier + dispatcher live in `src/rat.c` (NOT a new `poly/ratcanon.c` —
      the reused helpers `flint_cancel_fraction`/`flint_gaussian_*`/`extract_num_den`
      are static to rat.c, so a separate module can't reach them). Functions:
      `rat_canon_classify`, `rat_canon_dispatch`, `rat_canon_shadow`, `rc_math_equal`.
- [x] No CMake change needed (rat.c already in the build).
- [x] Shadow hook in `builtin_cancel`/`builtin_together` wrappers gated by
      `MATHILDA_RATCANON=shadow` (env cached; reentrancy-guarded; atexit summary).
      Oracle = numerator of `Together[a−b]==0`, RootReduce fallback for cyclotomic
      disguised-zeros. Only DIVERGENCE lines + summary printed.
- [x] Parity proven: DIVERGED=0 across risch_transcendental(7919 cov),
      goursat(859), simplify(144), rat/poly/trigrat/radical. All 46 initial
      goursat "divergences" were false positives (sign convention + cyclotomic
      disguised-zeros, verified 0 by RootReduce). Normal-mode suites unchanged.

## Phase 1 — Field-structure classifier + dispatch  [PARTIAL]
- [x] One-pass classify: PLAIN_Q | TRANSCENDENTAL | ALG_CONST | ALG_FUNC | DECLINE.
- [x] PLAIN_Q + TRANSCENDENTAL: real FLINT dispatch (flint_rational_*), proven parity.
- [ ] ALG_CONST + ALG_FUNC dispatch is INCOMPLETE — declines many Plus shapes the
      old cascade handles (e.g. Together[1/(b-a√k)+1/(b+a√k)], Cancel[…√2…,Ext->Auto]).
      Declining is safe (falls back to classical, no divergence) but must be COVERED
      before Phase 2 so the classical algebraic paths can be deleted. Next: trace how
      the old cascade routes these (flint_cancel_fraction gate on a Plus, the
      autodetect/together_recursive_ext number-field combine) and reproduce as single
      dispatch. Then run FULL suite in shadow to expand coverage.
- [ ] Fix `flint_algebraic_field_normalize` unbounded recursion (Q(√d) d≥5 SIGSEGV).
      NOTE: 1/(x^2-5)+1/(x-Sqrt[5]) and d=7 did NOT crash via the old path in testing;
      need the exact partial-fraction trigger from TOGETHER_ALGEBRAIC_OVERFLOW.md.
- [ ] Together vs Cancel read-back policy (combine vs per-additive-term) — dispatch
      currently uses flint_rational_together/cancel gates; verify per-term Cancel.

## Phase 2 — Switch builtins  [DONE 2026-07-22]
- [x] `builtin_cancel_compute`/`builtin_together_compute` call `rat_canon_dispatch`
      first (gated `rat_canon_enabled()`, default ON; `MATHILDA_RATCANON=off` reverts).
- [x] `rc_sign_normalize`: match classical denom-sign convention (flip only the
      FLINT all-negative-denom -(N)/(-D) artifact; term-wise `rc_negate`).
- [x] All direct suites pass; every integrate suite identical off-vs-on; simplify
      soft-fail 4=4 pre-existing; valgrind definitely-lost identical off/on (macOS
      dyld baseline only). User-facing shadow formdiff=0.
- Result: `Together[1/(x-I)+1/(x+I)]` -> `(2x)/(1+x^2)` (was ugly intermediate).

## Phase 3 — Delete dead cascade paths  [BLOCKED — needs broader dispatch coverage]
- Dispatch declines many regimes (alg-const rationalization, cube roots, poly
  radicands) to the classical cascade, which is STILL the needed fallback.
  Deleting the cascade requires those regimes covered by dispatch first (FLINT
  engines producing classical-matching forms). Do incrementally, shadow-verified.

## Phase 4 — Simplify: drop duplicated top-level autodetect+Together cascade. Pending.

## Also pending
- Fix `flint_algebraic_field_normalize` Q(√d) d≥5 SIGSEGV (exact trigger from
  TOGETHER_ALGEBRAIC_OVERFLOW.md; `1/(x^2-5)+1/(x-Sqrt[5])` did NOT reproduce).
- 509 internal risch formdiffs (Gaussian coeff-clearing (1/2 I)/(I+Log) ->
  I/(2I+2Log)) — harmless (integrate suites off==on); reconcile if Phase 3 needs it.

## Verification anchors
- `Cancel[(a+b Sqrt[k]-a^2 k)/(b-a Sqrt[k]+b^2 k-a^2 k^2), Extension->Automatic]` <1s
- `Together[1/(b-a Sqrt[k])+1/(b+a Sqrt[k])]` -> (2b)/(b^2-a^2 k)
- `Together`/`Cancel`/`Simplify` over Q(Sqrt d), d=5..8 — no SIGSEGV
- Regression: `Cancel[(1+Sqrt[2])/(1-Sqrt[2]),Extension->Automatic]` -> -3-2Sqrt[2];
  cyclotomic unchanged; `Together[(E^(2x)-1)/(E^x-1)]` -> 1+E^x; no Expand[Power[Plus,n]]
- Tests: rat/poly/simp/integrate `*_tests` only; valgrind on reproducers.
