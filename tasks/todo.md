# Goursat cube-root third-kind (p=2/3) + general-affine, FLINT-accelerated

Goal: make these close under `Method -> "GoursatAlgebraic"`:
- (a) `Integrate[((x-1)^2(x+1))^(1/3)/x^2, x]`  — genus-0 double-root radicand, double pole
- (b) `Integrate[(1-k x)/((1+(k-2)x)(x(1-x)(1-k x))^(2/3)), x]` — third-kind, symbolic k

## Diagnosis (done)
Both are THIRD-KIND (F pole off the branch locus). Eigendescent correctly obstructs
(H0 genuinely != 0). `goursat_cubic_thirdkind` was p=1/3 only. Native canonic over the
descent tower Q(k)(k^(1/3), R(t0)^(1/3), omega) does not terminate -> route through the
FLINT algebraic-field engine (`flint_algebraic_field_normalize` = rigorous zero test,
`flint_algebraic_field_canonical` = RootReduce). FLINT engine confirmed to handle
symbol-radicand cube roots + roots of unity together. Constant radicands (numeric k =
3^(1/3)) are DEFERRED by the engine -> target symbolic k (what the user wants).
diff_back_ok is already numeric (pins params + x) so it is NOT a tower bottleneck.

## Plan
- [ ] P1. Extend `goursat_cubic_thirdkind` to p=2/3: coefficient direction omega^(a*j),
      a=3-pnum; matched Y-slot Y^(3-pnum); vanishing slots = the other two;
      C = F*N/(R*p_match). [math verified in experiment]
- [ ] P2. General-affine offset: t0 = F's single non-branch pole; if N(t0)=0 keep m=0
      (fast, preserves existing p=1/3 cases); else m = R(t0)^(1/3) - kappa t0, a = kappa t+m,
      N = R - a^3. [math verified: m_off = 3^(1/3)+2(-1)^(1/3) for case b]
- [ ] P3. FLINT routing: new static helper `tk_is_zero_tower(e)` = try
      `flint_algebraic_field_normalize` first (->0 => zero), fall back to canonic+is_zero.
      Use for the two vanishing-slot tests. Route the C computation (F*N/(R*p_match))
      through `flint_algebraic_field_canonical` when a genuine radical generator is
      present; keep native canonic otherwise / when FLINT absent.
- [ ] P4. Test (b) symbolic k end-to-end; confirm diff_back_ok closes; valgrind clean.
- [ ] P5. Case (a): eigendescent double-root path in `goursat_cubic` (accept dR=3,nr=2 ->
      2 finite roots + Infinity). Verify it closes (genus 0, kappa=1). If third-kind
      route also needed for the double pole, assess.
- [ ] P6. Regression: existing p=1/3 third-kind + full goursat/rat/simp suites; add
      GOURSAT_EXERCISES cases; docs (algebra.md) + changelog; memory update.

## Review (2026-07-04)
- P1/P2/P3 IMPLEMENTED and verified sound: p=2/3 coefficient direction + matched slot,
  general-affine m = R(t0)^(1/3)-kappa t0, and FLINT-routed zero test / canonic
  (tk_is_zero_tower via flint_algebraic_field_normalize; canonic_tower via
  flint_algebraic_field_canonical). FLINT routing SOLVES the tower non-termination:
  native canonic did not finish at 300s/pass; FLINT-routed slots return in <1s and are
  numerically correct (pzeroA~8e-17, pzeroB~3e-19, pmatch~0.13).
- BLOCKER (mathematical, not code): for case (b) the scale C = F*N/(R*p_match) is NOT
  t-free (Craw@t=1/7=0.686 vs @t=2/9=0.701) -> |D[G]-f|=7.24. A single linear-argument
  log family cannot capture the integrand: F = -k/(k-2) + pole-part, whose constant part
  integrates to an ALGEBRAIC term and pole-part to a Log (two eigencomponents).
- The reference's constructive algorithm (preprint 30042026B Cor 4.8(i) / CubicTest23) is
  the eigendescent substitution x=z^3 -> rational J1,J2, valid only when H~0 == 0. Case (b)
  has H~0 != 0 (measured 0.19-0.81i). Preprint L838-840 EXPLICITLY DEFERS the H~0 != 0
  rational-F (third-kind) case to Risch-Trager-Bronstein; the two-component split is an
  UNPROVEN assertion in cube_roots.tex with no criterion/formula. Case (a) additionally
  has a repeated root (violates every theorem's simple-root hypothesis).
- CONCLUSION: neither example is closeable by the reference's constructive method; both
  need a research-level cube-root third-kind Risch-Trager-Bronstein (or the unconstructed
  split). All experimental edits REVERTED; working tree clean. Awaiting user direction.
- Reusable: tk_is_zero_tower / canonic_tower (route Goursat tower canonic/zero-test through
  the committed FLINT field engine) is a proven pattern worth landing independently.

## RootReduce WL-faithfulness (G1–G5) — 2026-07-04 — DONE

Plan: /Users/user/.claude/plans/glowing-petting-pnueli.md

- [x] G1 constant algebraic numbers → Root/quadratic/rational (FLINT qqbar, new src/poly/flint_qqbar.c)
- [x] G2 nested constant radicals reduce ((Sqrt[18]+Sqrt[27])/Sqrt[5+2Sqrt[6]] → 3)
- [x] G3 RootReduce::argx / RootReduce::mtd diagnostics
- [x] G4 thread over Equal/inequalities/logic; exact qqbar decision for binary (in)equalities
- [x] G5 Method → Automatic/Recursive/NumberField distinct paths (identical canonical result)
- [x] src/rootreduce.c dispatcher; parametric towers still route to flint_algebraic_field_canonical
- [x] tests/test_rootreduce.c extended (G1–G5); all pass
- [x] docs/spec/builtins/algebra.md + changelog 2026-06-29.md updated

Review:
- WL-faithful Root index (reals asc, then non-reals) — matched WL exactly on every
  reference example, including the degree-6 (2^(1/3)+Sqrt[5]) and degree-15 Root-product.
- Fixed a real bug found while wiring Root-object input: build_fmpz_poly reused a temp
  polynomial across Plus/Times terms without zeroing in the variable branch → wrong minpoly
  (t^3+t+c → 2t^3+…) and a qqbar_add blow-up that manifested as a hang on the c→1 example.
- Behaviour change (more WL-faithful): degree≥3 constant results are Root objects, not cubic
  radicals; updated the one existing test that expected 1/(1+2^(1/3)+2^(2/3)) → -1+2^(1/3).
- Memory: paired init/clear audit clean; no my-code frame in valgrind definitely-lost.
