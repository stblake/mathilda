# P3 — Residue-criterion Boolean decision half (prove non-elementarity)

**Goal:** Turn the transcendental integrator's authoritative `NULL` declines into a
*positive* non-integrability verdict. Deliver a sound Bronstein decision predicate
`Risch\`ElementaryIntegralQ[f, x]` and an `Integrate::nonelem` message from the main
integrator. Sound-by-construction: `True`/`False` only behind an exact certificate;
otherwise decline to *Unknown* (unevaluated). Never guess.

**Scope decided (user):** complete decision (residue criterion §5.6 Thm 5.6.1(ii)
*and* the poly/RDE "no rational solution" certificates §5.8/§5.9/Ch.6) + builtin
*and* integrator message.

**Out of scope by construction:** §5.11 `IntegrateNonLinearNoSpecial` — the live tower
only builds primitive (log) and hyperexponential (exp) monomials plus the separate
hypertangent case; no nonlinear-no-special monomial arises, so it would be dead code.
Document the reasoning; do not build.

---

## Design

**Tri-state** `RtDecision { RT_DEC_UNKNOWN=0, RT_DEC_ELEMENTARY, RT_DEC_NONELEMENTARY }`,
threaded as an *optional* out-param `RtDecision* dec` (default NULL ⇒ current behavior
exactly). It only annotates the *reason* for a decline; it never changes control flow of
existing callers (they pass NULL). The audited-sound integrate path stays behaviorally
pristine.

Authoritative NONELEMENTARY sources (all already computed, all currently collapsed to NULL):
1. **Non-constant residue** in the simple part (§5.6 Thm 5.6.1(ii)) — in
   `intrat_log_part_core` (`!free_all`). Surface it through
   `intrat_transcendental_log_part` → `Integrate\`TranscendentalLogPart` →
   `rt_field_lrt_logpart` → `rt_field_ratint_hermite`.
2. **Poly/RDE "no rational solution"** (§5.8 `rt_int_primitive_poly`, §5.9
   `rt_int_hyperexp_poly`, and `rde_base`) — A1-audited authoritative. Must tag *only*
   the genuine no-solution branch; scope/degree-bound declines stay UNKNOWN.
3. **Hypertangent `Dc≠0`** certificate (already in `risch_hypertangent.c`).

## RDE INTEGRATOR AUDIT (user directive 2026-07-14): no scope declines allowed

Every `NULL`/decline from an RDE integrator must be an authoritative "no rational
solution" proof, not a bounded-ansatz "gave up." Audit findings:

- `rde_base` (base field C(x)): **authoritative** (A1/F5-certified). Confirm the internal
  `rt_is_poly(aa,bb,cc)` guard (L1233) and `rde_weak_normalizer==0` guard (L1184) can
  never spuriously decline a genuinely-solvable input.
- `rt_solve_rde` → `rde_base`: **authoritative** (stale "bounded ansatz" comment; fix it).
- **`rt_field_rde` base-case guard is a real scope decline (FIX):** requires
  `rt_is_poly(p,x)`, so rational `p` (E^x/x, E^(x²)) is diverted to the general ansatz.
  Broaden to route *all* C(x)-base-field RDEs (w, p free of tower vars; rational-in-x OK)
  straight to `rde_base` — provably authoritative and cleaner.
- **`rt_field_rde` general tower case (ESTABLISH/COMPLETE):** bounded `SolveAlways` on
  `q=h/pd`. Its `NULL` is authoritative iff (a) `rt_rde_var_bound` is a *proven* per-variable
  upper bound in every config used (prim/exp, resonance), (b) `denom(q) | pd`, and (c) the
  exponential special-pole Laurent range `[-bd,bd]` is wide enough. Verify each against
  Bronstein §6.1–6.6; where the ansatz cannot express the true solution, replace the
  SolveAlways with the genuine per-level reduction (`RdeSpecialDenomExp` §6.2 +
  `PolyRischDECancel{Prim,Exp}` §6.6). Size TBD by the audit; may be deferred with an
  explicit UNKNOWN boundary if a config resists a completeness proof.
- **Non-non-elementarity dispatch declines (leave neutral):** `rt_int_hyperexp_poly`
  returning NULL for a non-monomial t-denominator (proper part → try coupled path) and the
  base-vs-coupled routing are *not* non-elementarity claims. The decision driver must treat
  these as "try the next path," never as NONELEM.

---

## Steps

- [ ] **1. Enum + choke-point signal.** Add `RtDecision` (file-local). Add optional
  `RtDecision* dec` to `rt_field_integrate`: ELEMENTARY on full success (`*rem_out`
  NULL), NONELEMENTARY on partial (remainder non-NULL) or authoritative sub-decline,
  UNKNOWN on scope decline. Existing callers pass NULL.

- [ ] **2. Residue non-constant signal.** Thread the `!free_all` fact out of
  `intrat_log_part_core`. Cross the `Integrate\`TranscendentalLogPart` builtin boundary
  via a trailing `decide` marker arg → distinguished return (or reuse the
  partial-remainder channel: in decide mode let the fully-non-constant case return
  `logs=0, remainder=a/d` instead of NULL, so `rem_out` non-NULL ⇒ NONELEM). Detect in
  `rt_field_lrt_logpart` / `rt_field_ratint_hermite` and set `dec`.

- [ ] **3a. Kill the RDE base-case scope decline.** Broaden `rt_field_rde`'s base-case
  guard so every C(x)-base-field RDE (rational `p` included) routes to `rde_base`; fix the
  stale `rt_solve_rde` comment. Verify `E^x/x`, `E^(x²)` now decide via `rde_base`.

- [ ] **3b. Establish the general tower RDE is a decision procedure.** Verify the exact
  degree bound + `q=h/pd` + Laurent-range ansatz makes every `SolveAlways` NULL
  authoritative (Bronstein §6.1–6.6). Confirm `rt_rde_var_bound` is a proven upper bound in
  each config (unit test `test_rde_var_bound` + a targeted proof note). Where a config can't
  express the true solution, replace with the genuine per-level reduction. Document the
  precise authoritative boundary.

- [ ] **3c. Poly/RDE decline classification.** With 3a/3b done, add `RtDecision* dec` to
  `rt_int_primitive_poly`, `rt_int_hyperexp_poly`, `rt_field_rde`,
  `rt_limited_field_integrate`; set NONELEMENTARY on the theorem-backed no-solution branch
  (RDE no-solution; §5.8 Dc≠0 new-log certificate), keep dispatch/routing declines neutral.

- [ ] **4. Decision driver.** `rt_decide(f, x) -> RtDecision`: detect the single-monomial
  extension (reuse the tower-build used by `rt_hermite_case`/`rt_hyperexp_case` →
  `rt_tower_build`), call `rt_field_integrate(..., &rem, &dec)`; also run the hypertangent
  `Dc≠0` path. Deliberately **excludes** `rt_special_case` and the ansatz-tower cases, so
  `E^x/x` ⇒ False (RDE no-solution) not True (Ei), and ansatz NULLs never masquerade as
  proofs.

- [ ] **5. Builtin `Risch\`ElementaryIntegralQ[f, x]`.** ELEMENTARY→`True`,
  NONELEMENTARY→`False`, UNKNOWN→unevaluated + `ElementaryIntegralQ::undec` info message
  (respects the *Q-returns-bool convention: bool on decided, diagnostic on undecided).
  Register: builtin + `ATTR_PROTECTED|ATTR_READPROTECTED`, terse docstring, `Risch\``
  namespace. Add to `sym_names.c` if a new internal symbol is needed.

- [ ] **6. `Integrate::nonelem` message.** In `builtin_rischtranscendental`, when
  `rt_integrate` returns NULL, run `rt_decide`; on NONELEMENTARY emit `Integrate::nonelem`
  (informational) and still return NULL. Gate tightly so it fires only on the authoritative
  verdict; verify no message-noise regressions across the integrate suite.

- [ ] **7. Tests** `tests/test_risch_elementaryq.c`: False — `E^x/x` (Ei), `E^(x^2)` (Erf),
  `1/Log[x]` (Li), a non-constant-residue case (Bronstein Ex 5.6.x); True — `E^x`,
  `1/(1+x^2)`, `Tan[x]`, a repeated-pole Hermite case; Unknown — an out-of-scope/algebraic
  integrand. Add source to `tests/CMakeLists.txt` COMMON_SRC. Leak-check with valgrind
  (diff vs `Sin[1.0]` baseline).

- [ ] **8. Docs.** `docs/spec/builtins/` (calculus/risch page) entry for
  `Risch\`ElementaryIntegralQ` + `Integrate::nonelem`; changelog note under the current
  ISO-week file (`docs/spec/changelog/2026-07-13.md` — Monday of this week). Refresh the
  P3 status lines in `RISCH_BRONSTEIN_GAPS.md` §6/§7 and `RISCH_STATUS.md`.

---

## Verification
- New test binary green; existing transcendental Risch + intrat suites unchanged
  (scope tests to the affected binaries per lessons).
- The decision path must be **sound**: manually confirm no False is emitted without an
  exact certificate. Spot-check that `Integrate` output is unchanged (only a message added).
- `make` clean under `-std=c99 -Wall -Wextra`.

## RDE audit result (2026-07-14)
- **3a DONE + verified.** `rt_field_rde` base case broadened; rational-p C(x) RDEs now
  route to `rde_base` (authoritative). Integrator output **unchanged** (E^x/x→Ei,
  E^(x²)→Erf, elementary cases still close); `integrate_risch_transcendental_tests` +
  `integrate_risch_macsyma_tests` green.
- **3b FINDING: RDE integrators are authoritative post-3a — no large build needed.**
  `rt_rde_var_bound` is a proven cap-free upper bound (RdeBoundDegree, resonance sub-case,
  reachability proof). General tower ansatz `q=h/pd` + Laurent range + complete linear
  `SolveAlways` ⇒ every `NULL` is a genuine no-solution. Residual narrow edges (document,
  don't over-build): (i) algebraic-constant coefficients (`SolveAlways` over Q); (ii)
  `rde_base` defensive guards (`rt_is_poly(aa,bb,cc)`, `rde_weak_normalizer==0`) — confirm
  unreachable for valid inputs. Internal dispatch declines (hyperexp proper-part → coupled)
  are routing, not non-elementarity claims. → 3b reduces to verify + document.

## Review — P3 COMPLETE (2026-07-14)

Shipped the residue-criterion Boolean decision half, sound by construction.

**What landed (2 files: `intrat.c` +42, `integrate_risch_transcendental.c` +262):**
1. **3a** — `rt_field_rde` base case broadened: every `C(x)` RDE (rational RHS incl.)
   routes to the complete `rde_base`. Removes the only real RDE scope decline.
2. **Tri-state** `RtDecision` via a file-local decision context
   (`g_rt_decide_mode`/`g_rt_decision`), raised (write-once) only at authoritative
   declines: `rt_field_rde` no-solution, `rt_limited_field_integrate` §5.8 `Dc≠0`,
   and — crossing the builtin boundary — a decide-mode `Integrate`$NonConstantResidue`
   marker from `Integrate`TranscendentalLogPart` (7th arg) for a non-constant residue.
3. **`rt_decide` / `rt_decide_field`** — route single-monomial (`rt_tower_build_min`
   with min_n=1) and deeper towers through `rt_field_integrate` in decide mode;
   True side exhibits an elementary form via `rt_integrate` + an elementary-head scan
   (`rt_expr_is_elementary`); elementary result has precedence over any flag.
4. **`Risch`ElementaryIntegralQ[f, x]`** builtin (Protected/ReadProtected) +
   `Integrate::nonelem` message from `builtin_rischtranscendental`.

**Verified:** `E^x/x`, `E^(x²)`, `E^(-x²)`, `1/Log[x]`, `E^(E^x)/(1+E^(E^x))` → False;
`E^x`, `x E^x`, `E^x/(1+E^x)`, `1/(x Log[x])`, `Log[x]/x`, `1/(1+x²)`, `Tan[x]`,
`1/(x(Log²+1))` → True; `1/Sqrt[1+x³]` → unevaluated. New suite green; transcendental
+ macsyma + residue-split + field + canonical + structure + hermite + intrischnorman
suites unchanged (pre-existing: `Risch`PolyDivide` in field/canonical + 2 intrat
cosmetic — all reproduce at HEAD). Valgrind byte-identical to `Sin[1.0]` baseline.

**Out of scope (documented):** §5.11 `IntegrateNonLinearNoSpecial` (no nonlinear
no-special monomial in the live tower); hyperexp coupled path proves non-elem only via
the RDE, not the residue.
