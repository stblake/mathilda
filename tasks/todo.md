# Thue-equation solver (general-degree, Baker + LLL) — todo

Plan: `/Users/user/.claude/plans/we-need-to-implement-agile-koala.md`
Contract: **provably complete, or DECLINE** (never a guessed answer).

---

## COMPLETION-PLAN M1 — Voronoi units, cubic (§3.1 of `docs/design/thue_completion_plan.md`)  ✅ DONE (2026-08-19)

**Goal.** Clear the Gate-2 declines where the fundamental unit's coordinates exceed
the coefficient-box search. Scope (decided from PARI): **rank-1 complex cubics only**
— the totally-real "Thomas family" declines are all *reducible* (x=−1 always a root),
i.e. correct declines, not targets; genuine small-reg totally-real cubics are already
box-found. Rank-2 Buchmann has **no validatable benchmark target** → deferred.

**Targets (monogenic, currently DECLINE → should become CORRECT):**
`Q(∛15)` reg 9.69 (fu coord 30 > box 12), `Q(∛42)` reg 11.06 (coord 42),
`Q(∛97)` reg 49.49 (~10-digit coords), `Q(∛41)` reg 56.29 (**12-digit coords**).
Non-monogenic d10/d12/d20/d999 stay declined (correct; that's M3 Round-2).

**Architecture leverage:** `nf_fundamental_units` is propose-then-certify. Voronoi
only PROPOSES the unit's integer coords; exact p-saturation + `|N|=1` is the
correctness boundary. A wrong proposal → certifier declines (safe), never wrong.

**Algorithm (rank-1 complex cubic, self-contained, no ideal-HNF subsystem):**
walk the chain of *second-kind adjacent minima* θ₀=1, θ₁, … of O_K, where
θ' = argmin σ₁ over {φ∈O_K : σ₁(φ)>σ₁(θ), |σ₂(φ)|<|σ₂(θ)|}; stop at the first
θ_k (k≥1) with |N(θ_k)|=1 → that is ε. Per step: LLL-reduce the fractional
ideal θ⁻¹O_K (well-conditioned neighbor search) via `lll_reduce_q` with the
transform recovered from an appended identity block (metric scaled by 2^P so the
identity passively records T ∈ GL₃(ℤ)); the reduced-basis O_K preimages g_m have
coords = rows of T; enumerate ψ=Σ eₘ gₘ over small e (|eₘ|≤3), pick the valid
neighbor with minimal σ₁. Track θ_k exactly in mpz coords (FLINT absorbs the
≤12-digit values); embeddings/norm via existing `nf_embed_int`/`nf_norm_int`.
arb precision adaptive (bump on inconclusive strict compares / as σ₁ grows).

- [x] `src/numbertheory/nfvoronoi.c` + decl in `numberfield_internal.h`
- [x] Refactor cert into shared helper (`cert_saturate`); wire Voronoi fallback in `nf_fundamental_units`; **certifier now mpz end-to-end** (24-digit units)
- [x] Add `nfvoronoi.c` to `tests/CMakeLists.txt` COMMON_SRC
- [x] Tests: `test_numberfield.c` (regulators vs PARI 6-digit, exact |N|=1), `test_thue.c` (rigorous x³−{15,41,42,97}y³=±1 → complete set) — ALL PASS
- [x] Benchmark 88: **CORRECT 48→56, 0 WRONG/0 CRASH**; check-c99 ✓, leaks=0, USE_FLINT=0 degrade compiles
- [x] Docs: solutions-of-equations.md, weekly changelog, plan docs + SOLVE_INTEGERS.md

## COMPLETION-PLAN M3 — Round-2 maximal order + O_K-basis unit search (§3.2)  ✅ DONE (2026-08-20)

**Scope (user-chosen):** general Round 2 (Pohst–Zassenhaus, any degree); cubics + quartics.
Targets (non-monogenic, currently DECLINE → CORRECT): complex cubics
`x³−{10,12,17,19,20}y³=1` (index 2–3; d17→{(1,0),(18,7)}, d19→{(-8,-3),(1,0)},
d20→{(-19,-7),(1,0)}), the Dedekind cubic `t³−t²−2t−8` (index 2), and
non-monogenic quartics `Q(d^{1/4})`.

**Key representation (minimizes regression risk):** O_K = `(1/D)·L`, L an integer
HNF lattice with rows B[i] in the θ-power basis; a field element is `(1/D)·v`,
`v = Σ c_i B[i] ∈ Z[θ]`. Then a **unit** `u=(1/D)v` has `|N(u)|=1 ⟺ |N(v)|=Dⁿ`,
and `embed(u)=embed(v)/D` — so the existing `nf_norm_int`/`nf_embed_int` on the
integer vector `v` are reused, only the target changes (`Dⁿ` not `1`). **Monogenic
= B=I, D=1 → byte-for-byte the current path.** The solvethue reconstruction stays
in θ-coords (via `coords_from_nfelem`) and is **unchanged**; only the unit search
walks the lattice L instead of Z[θ].

- [x] `nf_round2_maximal_order` (`src/numbertheory/nfround2.c`, FLINT matrices): p-radical
      Frobenius kernel, ring of multipliers via HNF, iterate; (W,D)+index+d_K. Used FLINT
      nmod_mat/fmpz_mat/fmpq_mat (no hand-rolled F_p nullspace needed)
- [x] O_K basis (W,D)+index+d_K in the struct; Gate 1 runs Round 2 on non-monogenic (monogenic unchanged)
- [x] Unit search over L (test |N(v)|=Dⁿ); Voronoi walks Z[θ] (result scaled D·w); saturation with D⁻¹
- [x] solvethue: build unit nf_elem as (1/D)·v; **Baker bound unchanged** (unit /D cancels in log-ratios)
- [x] Safety guard: M2 general-m gated to monogenic (nf_ok_index==1)
- [x] Validated: 130-case PARI grid **0 WRONG**; benchmark 88 CORRECT 65→81, 0 WRONG/0 CRASH;
      test_numberfield.c (index+d_K vs PARI incl. multi-round quartics), test_thue.c; c99, leaks=0, USE_FLINT=0

### Review — M3 (Round-2 maximal order + O_K-basis unit search)
- **Delivered:** non-monogenic `Solve[x^n-d·y^n==±1, Integers]` — `Q(∛{10,12,17,19,20})`,
  the Dedekind cubic, `Q(d^{1/4})` (incl. index-16 `Q(12^{1/4})`). Benchmark 88 65→81.
- **Round 2 (Cohen 6.1.3)** on FLINT matrices — the research agent's spec (validated on the
  Dedekind cubic end-to-end) was followed exactly. Key gotcha: `H_I` singular mod p ⇒ the
  β-change-of-basis is exact over Q, reduced mod p after. Frobenius exponent `p^k≥n`.
- **The low-risk representation:** O_K = `(1/D)L`; unit `u=(1/D)v` ⇒ `|N(u)|=1⟺|N(v)|=Dⁿ`,
  `embed(u)=embed(v)/D`, and the Baker bound's log-ratios cancel /D. So the reconstruction and
  bound are unchanged; only the unit-search lattice + norm target change. Monogenic = identity.
- **False alarm caught:** truncated THUE_DEBUG made x⁴-12 look under-maximized (index 4); the
  direct round2 test proved index 16 (correct). Lesson: validate the algorithm via a direct
  unit test, not truncated Solve-path debug.
- **Files:** `src/numbertheory/nfround2.c` (new), `numberfield.{c,h}`+internal (Round-2 wiring,
  O_K struct, accessors), `nfunits.{c,h}` (O_K-lattice search, D-aware saturation), `solvethue.c`
  ((1/D)v build + monogenic guard), tests, docs.
- **Follow-ons:** O_K-Voronoi (large-reg non-monogenic quartics), general m × non-monogenic (M2×M3).

---

## COMPLETION-PLAN M2 — General m (|m|≠1) via μ-enumeration (§3.3)  ✅ DONE (2026-08-19)

- [x] `thue_norm_reps_cubic11`: enumerate norm-m reps (canonical box from fundamental
      domain + norm constraint via inverse Vandermonde; keep N==m; over-cover safe)
- [x] μ-aware `thue_exponent_bound`: δ += log|μ^(k)/μ^(j)| looped over μ; C4=/μ_, Y2p=*μ_+,
      V0 += μ-height — all OVER-estimates. NULL/0 → |m|=1 byte-for-byte unchanged
- [x] Per-μ enumeration loop in `thue_enumerate`; dropped the |m|=1 gate (rank-1 only)
- [x] Validation: 270-case PARI grid **0 WRONG**; benchmark 88 CORRECT 56→65, 0 WRONG/0 CRASH;
      test_thue.c (7 M2 cases); c99, leaks=0, USE_FLINT=0 degrade; no |m|=1 regression

### Review — M2 (general m via μ-enumeration)
- **Delivered:** `Solve[x³−d·y³==m, {x,y}, Integers]` for `|m|≠1` over rank-1 complex
  cubics — complete set or proven `{}`. `x³−2y³={2,3,10}` solve; `={4,5,9,73,100}` → `{}`.
- **Key idea:** `N(x−θy)=F(x,y)=m` ⇒ `β=μ·unit` with μ a norm-m rep. Enumerate μ (bounded
  box), reuse the whole unit-exponent engine per μ with the μ-ratio added to the linear form.
- **Safety discipline:** every μ-constant is an OVER-estimate — a too-large bound is safe
  (reduction shrinks it / box-check declines), only under-estimation misses solutions (WRONG).
  Validated by the 270-case PARI grid (the only way to trust the bound).
- **The two test breakages fixed:** the *bounded* + *rigorous* entries previously DECLINED
  |m|≠1; two old test asserts said so. Now they solve — updated the asserts (correct new behavior).
- **Scope:** rank-2 totally-real |m|≠1 (cyclic-cubic-m2) deferred to M2b (needs the 2-D
  fundamental-domain box). Declines safely.
- **Files:** `src/solvethue.c` (μ-enum + μ-aware bound + per-μ loop), `tests/test_thue.c`, docs.

### Review — M1 (Voronoi units, cubic)
- **Delivered:** `Solve[x³−d·y³==±1, {x,y}, Integers]` returns the complete set
  for d=15,41,42,97 (was DECLINE). Any monogenic complex cubic whose fundamental
  unit outgrows the coefficient box is now solved.
- **Key insight:** the engine is *propose-then-certify*, so Voronoi need only
  PROPOSE the unit; exact p-saturation + |N|=1 is the correctness boundary → a
  wrong proposal DECLINEs, never a wrong answer. Kept the box as the fast path;
  Voronoi is fallback-only.
- **The subtle bug fixed mid-build:** the second-kind adjacent minimum is *not* a
  short Minkowski vector (it minimises σ₁ under an anisotropic |σ₂|<1 box), so
  "LLL then small combos" found nothing. Fix: rescale the σ₁ column by 1/U=2^-u
  and grow U until the box captures a neighbour (isotropising the search).
- **Scope call:** rank-2 totally-real cubics (Buchmann 2-D) have no validatable
  benchmark target (Thomas family reducible; small-reg cases already box-found),
  so deferred to M4 — avoids shipping unvalidatable geometry.
- **Files:** `src/numbertheory/nfvoronoi.c` (new), `nfunits.c` (mpz certifier +
  fallback), `numberfield_internal.h` (decl), `tests/{test_numberfield,test_thue}.c`,
  `tests/CMakeLists.txt`, docs.

---

## M1 — Number-field layer + Gate 1 (maximal order, monogenic-first)  ✅ DONE
- [x] Read exact signatures
- [x] Export `lll_reduce_q` real-lattice wrapper in `src/linalg/latticereduce.c` + `linalg.h`
- [x] Export `facint_factor_complete` in `src/facint.c` + `facint.h`
- [x] `src/numbertheory/numberfield.{c,h}`: field setup, disc, Dedekind (mod-p factor radical), monogenic certification (Gate 1), exact Sturm signature
- [x] Unit test `tests/test_numberfield.c`: 3 monogenic fields accept; non-monogenic/reducible/non-monic decline; Dedekind + signature verified — ALL PASS
- [x] Build clean (GCC 16, FLINT 3.6). Fixed small-char radical bug (mod-p factorization, not gcd(f,f'))

## M2 — Unit-group engine + Gate 2 (p-saturation)  ✅ DONE
- [x] `src/numbertheory/nfunits.{c,h}`: small-norm unit search, log-embedding (acb), greedy independent set, regulator
- [x] p-saturation certification via mod-p character-matrix RANK (rank==r ⇒ saturated, unconditional); DECLINE if not reached (no fragile enlargement)
- [x] Validated: regulators match LMFDB — ℚ(∛2)=1.347377 (rank 1), cyclic cubic=0.849287 (rank 2), ℚ(2^1/4)=2.158001 (rank 2)

## M3 — Thue engine + dispatch  ✅ DONE
- [x] `src/solvethue.{c,h}`: reduce to unit eqn, enumerate exponents, reconstruct (x,y), verify exactly. Degree-generic.
- [x] Wire `si_solve_thue` into `solveint_solve_integer`; expose `Solve`ThueSolveForm` test builtin
- [x] **Rigorous Baker(Waldschmidt)+de-Weger bound** (`thue_exponent_bound`): C1..C6 constants, Waldschmidt K3, de Weger LLL reduction (Prop 3.2), Q-dependent case (iii) via relation-detection + L-trick. arb/acb @1600 bits.
- [x] Downloaded the reference papers to `docs/references/thue/` + `ALGORITHM_NOTES.md`; implemented against the explicit constants.
- [x] **Solve[…,Integers] now returns the complete sets** for all three targets + more (Thomas 9 solns, x^3-7y^3=1). Non-monogenic / |m|≠1 decline safely.
- [x] `tests/test_thue.c` rigorous-path + reconstruction tests; 0 leaks (MSL); c99 clean

## M4 — Degree-4  ✅ came free (engine is degree-generic; x^4-2y^4=-1 solves)

## M5 — hardening/docs  ✅ DONE
- [x] Decline paths, degrade build, check-c99, docs (spec §6-F, changelog), no leaks
- [ ] (optional follow-on) Held-out gate extension; general a0/|m|; non-monogenic Round-2

## Review
- **DELIVERED (end-to-end):** `Solve[F(x,y)==m && Element[{x,y},Integers], {x,y}, Integers]`
  returns the COMPLETE finite solution set for irreducible monic |m|=1 forms over monogenic
  fields, via the genuine Tzanakis–de Weger algorithm (Baker's linear forms in logs + LLL
  bound reduction). All three user targets solve; validated against known/brute-force sets.
- **Number-field layer:** Gate 1 (Dedekind maximal-order cert), Gate 2 (fundamental units +
  regulator via p-saturation — regulators match LMFDB exactly).
- **Contract:** provably complete or safe DECLINE (non-monogenic, |m|≠1, |a0|≠1, degree/precision
  out of reach). Never a guessed/incomplete answer.
- **Key implementation lessons:** (1) Dedekind radical from mod-p factorization, not gcd(f,f')
  (small-char); (2) real-case linear-form coeffs need FULL arb precision (double loses all
  digits past ~16, but c0~10^60); (3) subfield units → Q-dependent linear form → degenerate
  lattice → handle via relation-detection + L-trick (paper's case iii).
- **Files:** src/numbertheory/{numberfield,nfunits}.{c,h}, numberfield_internal.h;
  src/solvethue.{c,h}; si_solve_thue + include in src/solveint.c; solvethue_init in core.c;
  exports lll_reduce_q (linalg), facint_factor_complete; tests/{test_numberfield,test_thue}.c;
  docs/references/thue/ (papers + ALGORITHM_NOTES.md).

## Review
- (to be filled at the end)
