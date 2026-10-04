# Residue method: general algorithmic extensions

Plan: `~/.claude/plans/pasted-content-id-e7e1-the-following-humming-bachman.md`

Policy: tests stay correct-by-construction (closed-form pins; no committed
NIntegrate). I verify each case numerically in the REPL during dev. Bump
`src/version.h` + tag per substantive phase. Route warnings via `mth_message`.

## Phase 0 — PV plumbing (prereq)  [DONE v0.275]
- [x] Thread parsed `principal_value` into `integrate_residue_try` (new 7th param,
      carried via file-scope g_pv saved/restored per frame). integrate.c passes it;
      ContourResidue builtin passes false; even-half-line recursion inherits g_pv.

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
- [x] PV indentation for axis pole (Case 18)  [DONE v0.275]: a0_keyhole_log admits a
      SIMPLE axis pole as the average of the two keyhole-branch residues
      (1/2 Res(R,z0)[(Log z0)^k + (Log z0+2Pi i)^k]). Removable z0=1 (Log 1=0) needs no
      option; genuine z0!=1 needs PrincipalValue->True. C18 Log[x]/(x^3-1)=4Pi^2/27;
      Log[x]/(x^2-4) PV = Pi^2/8. Hand-derived + NIntegrate symmetric-exclusion checked.

## Phase 3 — Tier C parametrized contour + essential-sing  [DONE v0.271]
- [x] `residue_family_contour_param` (t->-I Log[u], G=F/(I u), 2 Pi i Sum Res
      inside |u|<1; finite-scalar gate since value is complex); `res_ess0`
      (w^1 coeff of Series[G/.u->1/w,{w,0,1}]). Gated on Exp presence; wired
      after residue_family_trig in the full-period branch.
- [x] Cases 19 (2 Pi I), 20 (2 Pi I(16E-128/3)), 22 (0); Case 9 delivered via
      the parametrized spelling (Case 19) -- literal Circle contour stays OOS.
      test_contour_param (FreeQ+Chop non-vacuous); docs + changelog; v0.271.
      Fixed a u=0 double-count (denominator root 0 vs the always-added 0 cand).

## Phase 4 — Tier D unit-circle order-n  [DONE v0.272]
- [x] build_instantiation FindInstance fallback for COUPLED assumptions (a>b>0 left
      `a` bounded only by `b` -> was declining). Fixes order-1 symbolic trig.
- [x] residue_family_trig symbolic branch: pole_order (counts vanishing derivatives
      at the instantiated root, since solve_roots dedups multiplicity) + analytic-part
      derivative residue (fast for radical poles; Series is >30s and form-sensitive).
      Case 16 1/(a+b Cos)^3 -> Pi(2a^2+b^2)/(a^2-b^2)^(5/2); order-1 -> 2Pi/Sqrt[a^2-b^2].
      test_trig_symbolic + order-2 numeric regression; docs + changelog; v0.272.

## Phase 5 — Tier E special contours
- [x] `residue_family_chebyshev_weight` (Case 8) v0.274: R(x)/Sqrt[1-x^2] on (-1,1),
      x=Cos[t] -> (1/2) trig-family. FIXES WRONG SIGN (NL gave -Pi/Sqrt2; correct +Pi/Sqrt2).
      Runs before NL so it pre-empts the bad antiderivative. test_chebyshev_weight.
- [~] Case 10 (sin(x^3)): ALREADY correct via integrate_fresnel.c (plain Integrate =
      Pi/(3 Sqrt[3] Gamma[2/3]) = Gamma[4/3]/2). Residue-path gap is cosmetic -> SKIP.
- [x] Case 7 x/Sinh = Pi^2/2  [DONE v0.275]: NO new family -- generalized
      residue_family_rectangular. After w=Exp[x], a Sinh/Cosh/... of Log[w] rationalises
      under TrigToExp+Together and a polynomial x^k becomes (Log w)^k, so x/Sinh[x] ->
      2 Log[w]/(w^2-1), routed to residue_family_mellin_log (removable w=1 pole). Pure-Exp
      strip path untouched (gated on a hyperbolic head). Works Automatic + pinned.

## Phase 6 — Case 11 Mellin–Barnes vertical line  [DONE v0.273]
- [x] CORE fix (src/plus.c classify_plus_term): a Times[Complex[0,b],Infinity]
      (non-real coeff) is a DIRECTED infinity, not real ±Infinity -> `1/2 - I Infinity`
      stays a Plus (was collapsing to Infinity, so the limits died before Integrate).
- [x] is_vertical_line + residue_family_mellin_barnes: const*Gamma[A+Bs]*X^(P+Qs),
      close left, 2 Pi i Sum_k Res at s=-(A+k)/B (Sum closes the series).
      Case 11 Gamma[s]x^-s -> 2 Pi I e^-x (was silent 0!); Gamma[2s]x^-s -> I Pi e^-Sqrt[x].
      test_mellin_barnes + test_directed_infinity; docs + changelog; v0.273.
      NOTE: declining vertical-line integrals (e.g. 1/s) emit a cosmetic
      "0 Infinity" message from a downstream cascade stage -- harmless, out of scope.

## Phase 7 — Tier F hard tail
- [x] Gaussian-Fourier recognizer (15)  [DONE v0.276]: residue_family_gaussian --
      shifted-rectangle (entire E^(quadratic) kernel). TrigToExp+Expand -> Sum of
      E^(linear); complete the square per term: Int e^(A x^2+B x+D) = Sqrt[-Pi/A]
      e^(D-B^2/(4A)), Re A<0 (sign verified over the guaranteed region). Runs FIRST
      on the whole-line branch (gated on contains_exp_of_var), so it also fixes the
      HANG (was spinning the nonelem antiderivative search). C15 E^-x^2 Cos[2ax] =
      Sqrt[Pi]e^-a^2 (symbolic + concrete + Automatic). test_gaussian.
- [x] 17 & 21 documented honest declines  [DONE]: both decline CLEANLY under strict
      Method->"Residue" (no hang, no wrong value). C17 Log[x]/Cosh[x] = Mellin-of-sech
      differentiated at s=1 -> Dirichlet-beta derivative (Gamma[1/4] constant), NOT a
      residue sum; w=Exp[x] fails (Log x -> Log[Log w]). C21 Hankel fragment ambiguous/
      divergent. Pinned as negative controls (test_honest_declines) + docs note. Tests-
      and-docs-only (no src change) -> NO version bump/tag per CLAUDE.md.

## Review
_(filled in as phases complete)_
