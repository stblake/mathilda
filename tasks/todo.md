# FindInstance stress-test round: 10 hard cases

Goal: general, algorithmic, **sound** (verify-gated) methods so FindInstance
solves the 10 stress cases instead of declining.

## Diagnosis (empirically confirmed)
- **1** `Sqrt[z^2]!=z`, ℂ — Reduce/Solve decline; z=-1 already verifies True. Need a sampling fallback.
- **2** `Log[x y]!=Log[x]+Log[y]`, ℂ — need sampling AND `0 != 2 I Pi` must evaluate True.
- **3** `Sin[x+I y]==2`, ℝ — needs exact symbolic trig-system solving; numeric witness fails exact `==`. (defer)
- **4** `(x^3-1)/(x-1)==0 && x!=1`, ℂ — Reduce RETURNS the roots; only verification fails (`frac==0`, `root!=1` don't fold).
- **5** `ab+bc+cd+de==0 && a^2+..+e^2==5`, ℤ — integer box unbounded; the sum-of-squares equality is ignored.
- **6** returns `{}` — already correct.
- **7** two float eqns + `x!=y` + `a>0`, ℝ — needs robust numeric feasibility with tolerance verify (inexact input). (stretch)
- **8** `Tan[x]==x && x>10^6`, ℝ — no transcendental root search; a bracket root near (n+½)π verifies True.
- **9**, **10** already work.

## Plan

### Core (high confidence, sound, general)
- [ ] **A. `src/comparisons.c` — Equal/Unequal exact zero-test fallback.**
  When operands are not structurally equal and `compare_numeric` cannot decide,
  but BOTH are numeric constants (`is_numeric_quantity`), decide via
  `zero_test_decide(a-b)`: TRUE⇒equal, FALSE⇒unequal. Gate on numeric-constant
  so free-symbol `x==0` stays unevaluated. Expose `is_numeric_quantity` from
  core.c. Fixes `I==0`→False, `I!=0`/`2 I Pi!=0`→True, and cases **2**, **4** verify.
- [ ] **B. `src/solve/reduce_companions.c` — complex/real sampling witness layer.**
  Verify-gated candidate grid over ℂ/ℝ (values incl. -1,1,±2,±I,1±I,0,1/2,…),
  cross-product for small nv, capped otherwise. Late fallback only. Fixes **1**, **2**.
- [ ] **C. `src/solve/reduce_companions.c` — sum-of-even-powers integer bounding.**
  Derive per-variable bounds from an equality `Σ cᵢ vᵢ^(2k) (+ nonneg) == C`
  (cᵢ>0): `|vᵢ| ≤ (C/cᵢ)^(1/2k)`. Feed into the integer box. Fixes **5**.
- [ ] **D. `src/solve/reduce_companions.c` — 1-var real transcendental root scan.**
  Single real var + equality + a bound: scan `h=lhs-rhs` from the bound for a
  finite-endpoint sign change (skip pole spikes), bisect, verify. Fixes **8**.

### Stretch
- [ ] **E. case 7** — numeric feasibility (residual-minimising) + numeric-tolerance
  verify, gated to inexact-input ℝ systems.
- [ ] **F. case 3** — documented as future work (complexify-split is the entry point;
  needs a real trig-system solver).

### Verification & housekeeping
- [ ] Run the 10 cases; diff vs before.
- [ ] `tests/test_reduce.c` / `test_comparisons` / `test_find_instance` — add asserts; run full suite.
- [ ] valgrind the new paths.
- [ ] Docs: `docs/spec/builtins/*`, weekly changelog, version bump.
- [ ] `make check-c99`.

## Review (v0.107, 2026-08-27)

**Outcome:** 9/10 stress cases solved. Case 3 (`Sin[x+I y]==2`, exact transcendental
system) deferred — needs a real trig-system solver; complexify-split is the entry point.

- **A. `src/comparisons.c`** — `Equal`/`Unequal` fall back to `zero_test_decide(a-b)`
  when `compare_numeric` can't order and both sides are `NumericQ`
  (`expr_is_numeric_quantity`, newly exported from `src/core.c` via `src/numeric.h`).
  Fixes `I==0→False`, `I!=0`/`2 I Pi!=0→True`, and cases **2**, **4** verification.
  Free-symbol `x==0` stays symbolic (soundness gate).
- **B. sampling** (`fi_sample_search`, `fi_make_candidates`) — ℂ/ℝ candidate grid,
  verify-gated last resort. Cases **1**, **2**.
- **C. sum-of-even-powers bounding** (`fi_sos_upper_bound`, `fi_sos_solve_m`) — bounds
  the integer box from `Σ cᵢ vᵢ^(2k)==C`. Case **5**.
- **D. 1-var real root scan** (`fi_real_root_search` + `fi_smooth_residual` pole-clearing
  + high-precision bisection). Case **8**.
- **E. residual feasibility** (`fi_numeric_feasibility` + `fi_num_true` tolerance verify)
  for inexact-input systems; transc gate extended (`PolyGamma`/`LogGamma`/`Zeta`/`E^`).
  Case **7**.

**Verification:** full test suite green (comparisons, reduce, reduce_corpus 0/158,
solve*/root*/rootreduce/ratcanon/numeric*/simplify/…); new asserts in
`test_comparisons.c` + `test_reduce.c`; one stale test updated (`x!=0` now witnesses,
matching Mathematica); reduce-corpus verifier taught to skip out-of-domain sample points
(a pole / `Sqrt` of a negative). `make check-c99` clean; valgrind at macOS baseline (no
Mathilda frames in leak stacks). Timing: all 10 cases in ~1.1s.
