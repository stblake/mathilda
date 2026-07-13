# Bronstein transcendental Risch — audit + forward implementation

Two tracks, both grounded in Bronstein *Symbolic Integration I* 2nd ed.
(PDF in repo root; book page N ≈ PDF page N+17).

## Track 1 — Faithfulness audit of existing integrate code (user directive)
"All existing integrate code should be thoroughly reviewed against Bronstein's
book for inconsistencies. When in doubt, re-implement using the book."

Read-only audits (subagents), each producing a severity-rated findings list in
`RISCH_AUDIT_FINDINGS.md`:
- [ ] A1 — Ch.6 Risch DE stack: `rde_base` + `rde_spde` + `rde_weak_normalizer`
      + `rde_polyrischde_nocancel1/2` + base degree bound (integrate_risch_transcendental.c
      L828–1272) vs §6.1–6.6 (WeakNormalizer, RdeNormalDenominator, RdeBoundDegree,
      SPDE, PolyRischDENoCancel1/2). Correctness-critical.
- [ ] A2 — Ch.5 reductions/residue: Hermite (`rt_hermite_try`), Rothstein–Trager/LRT
      (`rt_frac_lrt`, `rt_field_lrt_logpart`), primitive/hyperexp poly cases
      (`rt_int_primitive_poly`, `rt_int_hyperexp_poly`, `rt_log_poly_case`,
      `rt_exp_poly_case`) vs §5.3, §5.6, §5.8, §5.9.
- [ ] A3 — Rational function integration `intrat.c` vs Ch.2 (Hermite + Lazard–
      Rioboo–Trager / Rothstein–Trager), and the new foundation modules
      (`risch_field`, `risch_canonical`, `risch_structure`, `risch_hypertangent`)
      vs §3.4–3.5, §9.3, §5.4, §5.10.
- [ ] Consolidate findings; fix / re-implement inconsistencies (book as basis).

## Track 2 — Forward: coupled system + hypertangent reduced case (§5.10 + §8)
Verbatim algorithms extracted from the book (see agent report in session).
- [ ] Expose base-field Risch DE as `Risch`RischDE[f, g, x]` (wraps `rde_base`);
      verify it accepts Gaussian C(i)(x) coefficients. Standalone tests.
- [ ] `Risch`CoupledDESystem[f1, f2, g1, g2, x]` (a=-1): reduce to a single
      Risch DE over C(i)(x), split Re/Im. Book §8.1 "trivial route" — correct &
      real-valued for base field C(x). Tests vs Ex 8.4.1 sub-solves, Ex 5.10.2/3.
- [ ] `Risch`IntegrateHypertangentReduced[p, t, deriv]` (§5.10, p.169) on top of
      CoupledDESystem. Tests vs Ex 5.10.2 (sin x/x → no solution) and 5.10.3.
- [ ] (next) `CoupledDECancelTan` full Chapter-8 real recursion for deep towers;
      `IntegrateHypertangent` driver (Hermite+Residue+reduced+poly); route `tan`
      through it in the live integrator instead of TrigToExp.

## Review

### Track 1 — audit (DONE)
Three parallel Bronstein-faithfulness audits ran read-only over the integrate
code; consolidated in `RISCH_AUDIT_FINDINGS.md` (+ `A1/A2/A3`).
- **No CRITICAL / no soundness defects anywhere** — every closed form is
  certificate-gated, so wrong answers cannot be produced.
- Ch. 6 RDE stack, `intrat.c`, and the four foundation modules are faithful ports.
- Two MAJOR *completeness* gaps recorded → roadmap P1/P3: (A2-1) `rt_hermite_try`
  is a polynomial-coefficient ansatz, not literal `HermiteReduce` over the tower;
  (A2-2) the residue criterion discards the elementary `r_s` logs when the residue
  resultant mixes constant/non-constant roots.
- One follow-up: audit the FLINT fast path `flint_rde_base_solve_fg` (A1-F5).

### Track 2 — forward implementation (DONE)
- `Risch\`RischDE[f,g,x]` exposes `rde_base`; verified on Gaussian `C(i)(x)`
  (reproduces Bronstein Ex 5.10.3's solution exactly).
- `Risch\`CoupledDESystem[f1,f2,g1,g2,x]` (§8.1 reduction) — Ex 5.10.3 `(c,d)`,
  Ex 8.4.1 `(s1,s2)` exact; declines the `Ei`-obstruction cases.
- `Risch\`IntegrateHypertangentReduced[p,t,deriv]` (§5.10) — pole peeling with
  round-trip `p−D[q]∈k[t]` verified; Ex 5.10.2 `∫ sin x/x` → non-elementary.
- Bug caught + fixed during testing: `rc_mult_t2p1` double-freed on a zero
  polynomial (0 divisible forever → runaway guard freed `q`==`cur`); added a
  `p==0` short-circuit and Together-based normalization (Cancel alone doesn't
  combine the monomial-derivation output, which broke the pole valuation).
- `tests/test_risch_coupled.c`: extensive + stress (constructed-solution
  round-trips, high-multiplicity poles, scaled/half-angle monomials). All green;
  sibling risch + transcendental suites still green; leak-clean.

### Track 3 — literal HermiteReduce over the tower (audit A2-1) (DONE)
- `src/calculus/risch_hermite.{c,h}` — Bronstein's quadratic `HermiteReduce`
  (§5.3, Thm 5.3.1) verbatim over the monomial derivation, exposed as
  `Risch\`HermiteReduce[f,t,deriv] -> {g,h,r}`, invariant `f = D[g]+h+r`.
  Reuses the now-exposed `risch_squarefree_t` (Yun) + `Diophantine`. Reproduces
  Example 5.3.1 exactly; handles arbitrary rational C(x) coefficients (the ansatz
  gap). `tests/test_risch_hermite.c` extensive (invariant across log/exp/tan,
  rational coeffs, higher multiplicity, h-simple/r-reduced, robustness). All
  green, no regressions, leak-clean.
- Fixed a helper leak during authoring (`rh_mul(expr_new_integer(j), c)` leaked
  the temp integer since `rh_mul` copies its args).

### Track 4 — wire literal Hermite into rt_field_ratint, remove ansatz (DONE)
- `rt_field_ratint` (log/primitive proper part) now delegates to a new
  `rt_field_ratint_hermite`: build the tower `RischDeriv` from `RtTower`, run
  `risch_hermite_reduce` → `(g,h,r)`, integrate the simple `h` via the residue
  criterion `rt_field_lrt_logpart`, return `g + logpart` with a tower diff-back
  self-verify. The ~150-line SolveAlways ansatz is **deleted** (spliced out).
- Full transcendental integrator suite still green with the ansatz gone;
  rational-coefficient repeated-pole integrals (`∫D[x/((x+1)Log x)]`,
  `∫(Log x−1)/Log²x`) now close. Leak-clean.
- The 2 failing `intrat.c` (Ch.2 rational) tests are PRE-EXISTING — confirmed
  identical at HEAD by stashing my changes and rebuilding (cosmetic string diff
  `1/(-a+x)` vs `-1/(a-x)`; corpus DIFF-NONZERO on parametric/Sqrt[2] rationals).

### Track 5 — wire literal Hermite into the exponential top, remove its ansatz (DONE)
- `rt_field_hyperexp_coupled` now delegates to `rt_field_hyperexp_hermite`, the
  literal Bronstein hyperexponential pipeline: `HermiteReduce` (§5.3) → residue
  logs on the simple part (§5.6) → coupling reconciliation `P = h + r − D_tower[L]`
  (for an exp monomial `D[Log g]=D[g]/g` is improper, so the logs spill a Laurent
  polynomial that this subtraction reconciles) → `IntegrateHyperexponentialPolynomial`
  (§5.9, `rt_int_hyperexp_poly`, per-coefficient Risch DE) → `g + L + Q`, tower
  diff-back self-verified. The unified `SolveAlways` ansatz is **deleted**.
- Strictly more complete: repeated exponential poles with rational lower-field
  coefficients now close (`∫ D[1/((1+x)(1+E^x)²)]`, `∫ D[x²/(1+E^x)³]`,
  `∫ x E^x/(1+E^x)²`). Both proper-part paths (log/primitive + exponential) are
  ansatz-free. New value-add cases added to `test_hyperexponential_case`; full
  transcendental suite green, leak-clean (zero new-code frames), intrat baseline
  unchanged (pre-existing 5-above-2, all Ch.2 rational).

### Track 6 — audit A1-F5 (FLINT RDE fast path `flint_rde_base_solve_fg`) (DONE)
- Audited `flint_rde_base_solve_fg` (`src/poly/flint_bridge.c` L791–880)
  line-by-line against Bronstein §6.1–6.5. **Verdict: faithful `fmpq_poly` port;
  every "0 = no solution" is authoritative** (RdeNormalDenominator guard `eₙ∤h²`
  by Cor. 6.1.1(ii); RdeBoundDegreeBase = proven upper bound p.199; SPDE declines
  = Thm 6.4.1; nocancel1/integrate peels decline only outside the degree window).
  `nr=1` solutions are exact by construction; `nr=-1` correctly defers when `f`
  is rational or `g` non-univariate-over-ℚ. No defect, **no code change**.
- Empirical: constructed poly & rational-`g` repeated-pole RDEs solve and
  self-verify `D[y]+fy−g≡0`; declines exactly on Bronstein's own non-elementary
  Ex. 6.1.1 (`eᵗ/t`) and Ex. 6.3.2 (`e^{−t²}`). Corroborated by the Expr fallback
  being a line-for-line mirror of the same stack. Track-1 audit now fully closed.
- Findings recorded in `RISCH_AUDIT_FINDINGS.md` (A1-F5 RESOLVED, Disposition).

### Track 7 — CoupledDECancelTan (Bronstein §8.4, tangent cancellation) (DONE)
Real recursion for the tangent-cancellation case, keeping the hypertangent
monomial `t` REAL (only the base field `k` goes complex, via CoupledDESystem's
C(i)(x) route) — this is the book's whole point vs the trivial C(i)(x)(t) route.
- `rc_cancel_tan(b0,b2,c1,c2,t,d,x,η,n)` (`src/calculus/risch_coupled.c`) = the
  p.265 box verbatim: base `n=0` → `CoupledDESystem(b0,b2,c1,c2)`; else eval
  `c1,c2` at `t=√-1`, Gaussian-split `z1+z2√-1`; `(s1,s2)←CoupledDESystem(b0,
  b2-nη,z1,z2)`; `c=((c1-z1+nη(s1 t+s2))+(c2-z2+nη(s2 t-s1))√-1)/(t-√-1)` split
  `d1+d2√-1`; recurse `(b0,b2+η,d1,d2,n-1)`; return `(h1 t+h2+s1, h2 t-h1+s2)`.
- Builtin `Risch`CoupledDECancelTan[b0,b2,c1,c2,t,deriv,n] -> {q1,q2}`
  (PROTECTED | READPROTECTED, docstring). Reproduces **Example 8.4.1** exactly
  → `{t-1, 2x}`, every intermediate `(z,s,c,d,h)` matching the book.
- Shared helpers extracted (`rc_split_gaussian`, `rc_base_var`), reused by
  `risch_coupled_desystem` + the hypertangent-reduced path.
- 5 new test groups: Example 8.4.1 + recursion exact; constructed-solution
  diff-back against the SYMMETRIC system `[[b0-nηt,-b2],[b2,b0-nηt]]` (forced by
  the complex-scalar structure of eq. 8.14 — verified 3 ways) over `tan(x)`,
  `tan(2x)`, `tan(3x)`, `tan(x/2)`, `n=0..4`, zero & nonzero `b0`; `Ei` decline;
  malformed / non-hypertangent NULL paths. Clean `-Wall -Wextra`, all risch +
  transcendental suites green, valgrind footprint == `Sin[1.0]` baseline.
Scope: base field `k = C(x)` (single tangent tower). Deeper `k` carrying more
tangent monomials needs the full CoupledDESystem dispatcher (mutual recursion) —
follow-up (see Next).

### Track 8 — retire TrigToExp for real tangent (Bronstein §5.10) — Phase 1 (DONE)
Real hypertangent integration framework so `tan` integrates directly (staying
real) instead of via `TrigToExp`→exp machinery→`ExpToTrig` (which strands
`Tan`/`Tanh` at an I-laden `I x − Log[1+E^(2Ix)]` no simplifier collapses).
- `Risch`ResidueReduce[h,t,deriv] -> {g2,β}` (§5.6 residue criterion over k(t),
  `src/calculus/risch_hypertangent.c`): LazardRiobooTrager via the existing
  `Integrate`TranscendentalLogPart`, lower-vars gate; β=False on non-constant
  residue (h non-elementary, Thm 5.6.1).
- `Risch`IntegrateHypertangent[f,t,deriv] -> {g,β}` (§5.10 p.172 driver):
  composes HermiteReduce → ResidueReduce → IntegrateHypertangentReduced →
  IntegrateHypertangentPolynomial exactly per the book box (`p=h-D[g2]+r`,
  `p-D[q1]`, `+c Log(t²+1)` when Dc=0). Reproduces Ex 5.10.1 `{tan(x),False}`,
  5.10.2 `{0,False}`, 5.10.3 `β=True`; genuine-log + non-constant-residue cases.
- 6 new test groups in `tests/test_risch_hypertangent.c` (book values, RR
  residue/diff-back/decline, constructed `f=D[g0]+base` β=True over tan(x)/tan(2x)/
  tan(x/2), robustness). Clean `-Wall -Wextra`, all suites green, leak == baseline.

### Track 8 Phase 2 — wire IntegrateHypertangent into the live integrator (DONE)
- `rt_hypertangent_case` (`src/calculus/integrate_risch_transcendental.c`),
  dispatched in `rt_transcendental_case` just BEFORE `rt_trig_frontend`: finds a
  single `Tan[u]` (u rational in x), substitutes `t=Tan[u]`, builds `Dt=u'(t²+1)`,
  runs `Risch`IntegrateHypertangent`, integrates the leftover `C(x)` element,
  back-subs `t→Tan[u]`, collapses `Log[1+Tan²]→-2Log[Cos]`, diff-back gated.
- Real forms now: `∫Tan→-Log[Cos]`, `∫Tan²→Tan-x`, `∫Tan⁴→x+Tan³/3-Tan`,
  `∫Tan[2x]→-Log[Cos[2x]]/2`, `∫x Tan²→½(-x²+2Log Cos+2x Tan)`, `∫2x Tan[x²]→
  -Log[Cos[x²]]`. Non-elem (`x Tan`, `x² Tan²`) proven & left unevaluated. Tanh
  (special t²-1) stays on TrigToExp.
- Normal-pole gate: defers denominators whose normal part (after stripping t²+1)
  has degree ≥2 to the fallback — irrational-algebraic residues swell the
  `p=h-D[g2]+r` reconciliation over the residue field (driver perf follow-up).
- 4 new test groups (clean forms by value; real-ness via `Position[g,Complex[__]]`
  + diff-back; powers 1-6, scalings, nonlinear args, special/linear-normal-pole
  rationals, mixed x·Tan^even; non-elem unevaluated; Tanh guard). Suites green,
  leak == baseline.

### Track 8 Phase 3 — Cot + Tanh variants (DONE)
- Cot: reuses `IntegrateHypertangent` (special t²+1, η=−u'); wiring only
  (cosmetic `Log[1+Cot²]→−2Log[Sin]`). `∫Cot→Log[Sin]`, `∫Cot²→−Cot−x`.
- Tanh: new `risch_integrate_hypertanh` (special t²−1 SPLITS → reduced case
  decouples into two real Risch DEs `DP−2mηP=a+b`, `DQ+2mηQ=a−b`, no ℂ(i)).
  Refactored `IntegrateHypertangent`→shared `rh_integrate_core` param'd by
  reduced/poly sub-builtins + special poly. New builtins `Risch`IntegrateHypertanh`
  `…Reduced` `…Polynomial`. `∫Tanh→Log[Cosh]`, `∫Tanh²→x−Tanh`, `∫Tanh⁴→…`.
- Live `rt_hypertangent_case`→family dispatcher + `rt_hypertan_family` core over
  Tan/Cot/Tanh. Non-elem (`x Tanh`, `x Tanh³`) left unevaluated. Sech/Csch/Coth
  out of scope (separate heads). 4 new test groups; suites green; leak==baseline.

### Track 8 Phase 4 — Coth (DONE); Weierstrass tangent path (REVERTED)
- `Coth` added: same hyperbolic monomial as Tanh (Dt=u'(1-t²), special t²-1),
  reuses the hypertanh driver; cosmetic `Log[1-Coth²]→-2Log[Sinh]` (sign folded,
  diff-back gated). `∫Coth→Log[Sinh]`, `∫Coth²→x-Coth`. `rt_hypertan_family`
  refactored to (subrule, bval, carg) so Tan/Cot/Tanh/Coth share one core.
- Weierstrass `t=Tan[x/2]` front-end for rational circular trig (Sec^n, Csc^n)
  was implemented and REVERTED: driver integrates the forms fast in isolation
  (∫Csc→Log[Tan[x/2]]) but has a perf pathology on HIGH-MULTIPLICITY SPLIT poles
  — `(1-t²)³` (Sec³) blows up the Hermite reconciliation while `t³` (Csc³) is
  fine. Made the suite 10× slower + 7 FAILs → backed out (with the gate
  relaxation, exp-route-in-trig_frontend, and numeric/rational verify helpers).

### Next (blocked on driver perf)
- **Driver perf on high-multiplicity split normal poles** — the prerequisite for
  BOTH relaxing the normal-pole gate AND the Weierstrass front-end (Sec^n/Csc^n
  real via the tangent path). Root cause: Hermite reduction / the `p=h-D[g2]+r`
  reconciliation swells for `(t-a)^m (t-b)^m` denominators (m≥2). Fix this, then
  re-add the gate relaxation (Factor-all-linear) + Weierstrass.
- Driver perf for irrational-residue normal poles (deg ≥2 irreducible).
- `Sech`/`Csch` (not rational-of-single-kernel — need Sec/Csc or Weierstrass).
- Full CoupledDESystem dispatcher for deep multi-tangent towers.
- Audit finding A2-2 (residue criterion `r_s` logs) — the last MAJOR completeness
  gap; requires the transcendental driver to carry a partial log part + an
  unintegrated `r_n` remainder (architectural, roadmap P3).
