# Residue method: general algorithmic extensions

Plan: `~/.claude/plans/pasted-content-id-e7e1-the-following-humming-bachman.md`

Policy: tests stay correct-by-construction (closed-form pins; no committed
NIntegrate). I verify each case numerically in the REPL during dev. Bump
`src/version.h` + tag per substantive phase. Route warnings via `mth_message`.

## Phase 0 — PV plumbing (prereq)
- [ ] Thread parsed `principal_value` into `integrate_residue_try` (new param)

## Phase 1 — Tier A correctness (FOUNDATIONAL)
- [x] Baseline-repro: Case 2 wrong; `Residue[z^(1/3)/(1+z^2)^2,{z,I}]` → 0 wrong
- [x] `series.c`: symbolic-centre `pad` scales with pole multiplicity (`series_max_neg_power`)
- [x] `residue.c`: two-order agreement loop; dropped-pole retry via shifted form
- [x] Audit Case 1 (matches NIntegrate); symbolic double/triple + z^(1/3)/z^(2/3) tests
- [x] version 0.268; docs (Residue Method para) + changelog
- [~] FULL test suite running (dsolve corpus slow); then commit + tag v0.268

### Phase 2 entry findings (from REPL prototyping, current binary)
- Keyhole-log algorithm CONFIRMED: `I_m = Limit[D[G(q),{q,m}], q->a]`, G=∫x^q R.
  Case 4 via stripped-Beta derivative gives exactly -Sqrt[3]Pi^2/18.
- CANNOT reuse `Integrate[x^q R]` builtin: wraps in ConditionalExpression (breaks
  D), and non-power denominators (Case 5, x^2+x+1) don't produce a Beta. Must
  build G(q) from the residue sum inside mellin_core, symbolic q, then D + Limit.
- Symbolic-exponent mellin currently DECLINES because `mellin_core` calls
  `param_interval(s)` with s=a+1 (a Plus, not a bare param) -> unbounded ->
  convergence gate rejects. FIX: derive s-interval from `param_interval(p)+1`.
  (g_inst IS built; `Inequality[-1,Less,a,Less,2]` is absorbed fine.)

## Phase 2 — Tier B keyhole-log + PV axis pole
- [x] symbolic-exponent mellin fix (param_interval(p)+1): **Case 12 solved**
      (`x^a/(x+1)^3` -> `Pi a(1-a)/(2 Sin[Pi a])`), Case 3 still OK
- [x] ALGORITHM VALIDATED (REPL): a=0 keyhole via (log z)^(k+1) contour +
      triangular solve from residues only. Case 4 -> -Sqrt[3]Pi^2/18,
      Case 5 -> 3.5361. (The differentiate-G(q) route is fragile for m=2 — DO
      NOT use it; the (log)^(k+1) contour is the robust method.)
- [x] IMPLEMENT `residue_family_mellin_log` (a=0 / integer p): klog branch,
      S_{k+1}=Sum Res[(klog z)^(k+1) R], triangular solve; Cases 4, 5 DONE
      (simple-pole shortcut for algebraic locations; order-2 via shifted form)
- [x] FIXED self-inflicted perf regression: param_interval s-vs-p asymmetry broke
      rectangular family -> 120s cascade hang (A/B caught it). Suite back to 2.1s.
- [x] v0.269, docs + changelog + lesson; tests green; COMMIT + tag
- [x] non-integer-a + log (Case 14): `keyhole_log_general_a` with (1-e^{2πia})
      triangular system. `Sqrt[x]Log[x]/(x^2+1)^2 -> Pi(Pi-4)/(8 Sqrt[2])`. v0.270
- [x] ROOT-CAUSE FIX (v0.270): `Arg` now evaluates on root-of-unity constants
      (`Arg[(-1)^(1/6)]->Pi/6`, incl. imaginary-unit + negative-base factors) in
      `src/complex.c`. This un-declined Cases 4 & 4b (were returning unevaluated,
      so their committed `Chop[N[...]]` tests passed VACUOUSLY via the N[]
      NIntegrate fallback). a0_keyhole_log closing switched from
      Simplify[RootReduce[Re[.]]] (buried poles in Root objects Arg couldn't
      touch) to Simplify[ComplexExpand[.]]. C4=-Sqrt[3]Pi^2/18, C4b=-Sqrt[2]Pi^2/16
      now CLEAN. Tests rewritten to assert FreeQ[r,Integrate] (non-vacuous).
- [ ] PV indentation for axis pole (Case 18)  [Phase 2b]

## Phase 3 — Tier C parametrized contour + essential-sing
- [ ] `residue_family_contour_param`; `res_ess0` primitive
- [ ] Cases 9(param), 19, 20, 22; tests/docs/bump/tag

## Phase 4 — Tier D unit-circle order-n
- [ ] Verify/extend family C under assumptions; Case 16; tests/docs/bump/tag

## Phase 5 — Tier E special contours
- [ ] `residue_family_fresnel` (10)
- [ ] `residue_family_hyperbolic` (7)
- [ ] `residue_family_chebyshev_weight` (8); tests/docs/bump/tag

## Phase 6 — Case 11 Mellin–Barnes vertical line
- [ ] `I*Infinity` / vertical-line classification; `residue_family_mellin_barnes`
- [ ] Case 11 → 2πi e^(-x); tests/docs/bump/tag

## Phase 7 — Tier F hard tail
- [ ] Gaussian-Fourier recognizer (15)
- [ ] Best-effort 17 & 21, else documented honest declines; tests/docs/bump/tag

## Review
_(filled in as phases complete)_
