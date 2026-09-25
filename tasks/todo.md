# ACTIVE (2026-09-25): Fix ParallelMixedTower branch defect on nested radicals

**Bug:** `Integrate[Sqrt[x + Sqrt[x]], x, Method -> "ParallelMixedTower"]` returned
`1/4 Log[x^(1/4) + Sqrt[1+Sqrt[x]]] + 1/12 Sqrt[1+Sqrt[x]](-3 x^(1/4)+2 x^(3/4)+8 x^(5/4))`,
correct only for x > 0. Root cause: `rootNormalize` (src/internal/mixed/ParallelMixed.m)
split the nested radical `Sqrt[Sqrt[x](1+Sqrt[x])]` into `x^(1/4) Sqrt[1+Sqrt[x]]` (the
branch-unsafe `Sqrt[a b]=Sqrt[a]Sqrt[b]`). Fix (approved: un-split genus-0 conics): the
pure-generator-factor split now runs only when needed for genus reduction (existing radical,
or fused squarefree radicand total degree >= 3); a genus-0 conic (deg <= 2, q === None) stays
fused as the simple radical y^2 = q, so back-sub y -> Sqrt[q] is globally faithful.

- [x] Edit `rootNormalize` split gate (src/internal/mixed/ParallelMixed.m ~1914-1926)
- [x] Reproduce: result now `1/8 Log[1/2+Sqrt[x]+Sqrt[Sqrt[x]+x]] + 1/12 Sqrt[Sqrt[x]+x](-3+2Sqrt[x]+8x)`
- [x] Global diff-back ~0 at x=-1/2, 2+I, 1/3-I/5, 7 (the OFF-positive-reals check the gate lacks); FTC [1,4] matches NIntegrate
- [x] PMT unit test (tests/build/parallelmixedtower_tests): all pass
- [x] Charlwood corpus A/B: baseline 48/50 (A19, A39) → WITH fix 48/50 (A19, A39). ZERO regressions; only P9, A35 change form (both to fused radical, both still verified)
- [x] All integrate_* + trigrat + parallelmixedtower unit tests pass
- [x] Extend tests/test_parallelmixedtower.c with nested-radical global-correctness case (free of x^(1/4) + diff-back at x=2+I)
- [x] docs/spec/changelog/2026-09-21.md note; per-category spec (calculus.md) update
- [x] Version bump src/version.h -> 0.192; rebuilt, $VersionNumber -> 0.192
- [ ] COMMIT + tag v0.192 — deferred to user (harness: commit only when asked); version.h already bumped for inclusion
- [ ] (optional) Sync untracked dev copy mixed/ParallelMixed.m — left to user (scratch, not runtime-loaded)

### Review
Root cause was the branch-unsafe pure-root-factor split in `rootNormalize`
(src/internal/mixed/ParallelMixed.m ~1914): it pulled `Sqrt[Sqrt[x]] = x^(1/4)` out of
`Sqrt[Sqrt[x](1+Sqrt[x])]`, building the whole tower in an (x>0)-only basis. The one-line
functional change gates that split on a genus test: keep a genus-0 conic radicand fused
(deg ≤ 2 in the generators, q === None) so it serves directly as the simple radical y^2=q
and back-substitutes globally; still split when needed to lower genus (deg ≥ 3) or when a
radical already exists. Result now `1/8 Log[1/2 + Sqrt[x] + Sqrt[Sqrt[x]+x]] + 1/12
Sqrt[Sqrt[x]+x](-3+2Sqrt[x]+8x)`, verified at negative and complex arguments. Blast radius
tiny (2 Charlwood cases re-form, both improvements). NB: the verify-or-decline gate still
samples only positive reals — it did not catch this and cannot catch a future branch defect;
the durable fix here is at construction, not the gate.

---

# DONE (2026-09-25 pm): Book updates per book/BOOK_UPDATE_REVIEW.md ("everything" run)

**COMPLETE.** All waves landed. `make examples` (197 transcripts), `make usage` (58 cards),
`make check-links` (1300 uses, all resolve), `make pdf` (240 pages, 0 errors, no undefined
refs/citations), Index regenerated (781 entries). Wave 2 (§4.3 Calculus+DSolve, §4.4 Linear
Algebra, §4.7 Special Functions) and Wave 3 (Ch.6 Data Structures, Ch.7 Programming, Ch.9
Data I/O, new §4.9 Geometry, new §4.10 Graphs) written by subagents, verified + integrated.
Wave 4: 04-mathematics.tex wired the 2 new sections; ROADMAP statuses/scope updated; changelog
note in docs/spec/changelog/2026-09-21.md; gen_usage.py ASCII-sanitize fix; site index
regenerated. NO version bump (book prose = contributor/reader-facing docs). NOT committed
(awaiting user). **Two kernel bugs surfaced (NOT fixed, out of scope):** Simplify[Laplacian[
Sin[r^2],{r,t},"Polar"]]→wrong `4-4r^2`; DownValues renders evaluated not held.

---

# (superseded) ACTIVE (2026-09-25 pm): Book updates per book/BOOK_UPDATE_REVIEW.md ("everything" run)

Bring the book current with Mathilda v0.113→v0.189 (anchor commit 8c4332ba). User chose the
full run: UPDATE the Verified math sections + WRITE the stub sections from scratch + open new
sections, all through the verified-example build (`make examples` → `\mtranscript`/`codepairs`,
never a hand-typed `Out[]`), `make check-links` clean, `make pdf` clean. Book prose does NOT bump
`$VersionNumber` (contributor/reader-facing docs, like SPEC.md prose). Changelog note in
`docs/spec/changelog/2026-09-21.md`. Shared-file edits (ROADMAP.md, TheMathildaBook.tex structure,
04-mathematics.tex \input lines, changelog, Index concepts) are done CENTRALLY (me), never by a
drafting subagent — one owner per section .tex file.

Wave 1 — UPDATE Verified chapters (me; also proves the build loop):
- [x] §4.8 Statistics: added `Quantile`/`InterquartileRange`/`MeanDeviation`/`MedianDeviation` (new "Quantiles and robust spread" subsection + examples/statistics/quantile.m + Pitfall on estimator conventions + \usagebox{Quantile}). Transcript verified.
- [x] §4.2 Algebra: number-field subsection (`ToNumberField`/`AlgebraicNumber*`/`NumberFieldIntegralBasis`/`AlgebraicIntegerQ`, \usagebox{ToNumberField}) — **FIXED the false "no AlgebraicNumber head" claim**; sparse-views subsection (`MonomialList`/`CoefficientRules`/`FromCoefficientRules`); division/ideal-membership subsection (`PolynomialReduce`); transcendental `Reduce` + `Simplify` over `Root`; `FactorSquareFreeList` mention. 4 new example files, all probe-verified.
- [x] §4.5 Numerical Calculus: added Interpolation subsection (`Interpolation`/`ListInterpolation`) + examples. Probe-verified.
- [x] §4.6 Number Theory: added `BitLength` to the Digits-and-bases subsection (only Bit* member so far). Probe-verified.
- [x] Site index regenerated (site/generate.py) → builtins.json now 922 builtins incl. all post-anchor heads; §4.8 check-links passes. (Prereq for every campaign's \B{} links.)
- [ ] Build gate 1: DEFERRED to the single integrated build after all drafting agents finish (global `make examples` races with in-progress `.m` files).

Wave 2 — WRITE stub math sections from scratch (subagents, one file each):
- [ ] §4.3 Calculus (D/Dt/Limit/Series/Residue/Integrate+Risch incl. ParallelMixedTower/Cherry/DiffUnderInt) + DSolve subsection
- [ ] §4.4 Linear Algebra (vectors/matrices/Dot/Det/Inverse/decompositions incl. Jordan/Schur/LU/QR/SVD/eigen/ZeroTest/LinearSolve/packed)
- [ ] §4.7 Special Functions (Gamma/Zeta/PolyGamma/Bessel/Airy/Erf/hypergeometric + LegendreQ)

Wave 3 — WRITE system chapters + new sections (subagents):
- [ ] Ch. 6 Data Structures (expressions/lists/associations/NDArray/packed + ArrayReshape/ArrayPad/ListGradient/FirstPosition/PackedArrayQ/ToPackedArray/ToNDArray)
- [ ] Ch. 7 Programming (pattern matching/procedural/functional + Inactive/Activate)
- [ ] Ch. 9 Data I/O (stream layer: Read/Open*/Write/WriteString/Close/Streams/StreamPosition/ReadList)
- [ ] New §4.9 Computational Geometry (Area/Perimeter/RegionCentroid/RegionMember/ConvexHullRegion) + Bitwise (BitLength)
- [ ] New Graphs section/chapter (StarGraph/FindVertexColoring/EdgeWeight/WeightedAdjacencyMatrix/weighted FindShortestPath+GraphDistance/RandomGraph)

Wave 4 — integrate & finalize (me):
- [ ] REVIEW behaviour-change prose (A1–A10 divergences, Function closure, Apart 2-arg) in §4.1/§4.2/Ch.3
- [ ] Update ROADMAP.md statuses + scope (DSolve under §4.3; Jordan under §4.4; new geometry/graphs; Ch.9 caveat)
- [ ] Regenerate Index concept `\index{}` entries; final `make pdf` clean
- [ ] Changelog note in docs/spec/changelog/2026-09-21.md

NOTE: NOT filling Ch.5 Graphics (needs figure-capture toolchain, deferred per ROADMAP), Ch.8 Compilation, Ch.10-13 (not in backlog).

---

# PAUSED (2026-09-25): mean-time push — user paused after v0.190

**Status:** the Charlwood mean-time push is PAUSED. Banked this session: **v0.190** (c4c36046,
tagged) native plain-Q `CoefficientRules` read-off (general ~27% win; Charlwood 23.7→22.7s, 48/50
unchanged) + **a2e71ca2** dormant TowerCRE kernel. The Maxima gap (12.6s) is essentially untouched —
each remaining lever is a large/uncertain restructure. **Resume pointers (per-case, `.m`-phase timed):**
- **A1/A2/A3** — the confirmed lever is threading native polys (TowerCRE, `flint_polynomial_monomials`
  pattern) through `AnsatzSystem`'s `parts` COLUMN CONSTRUCTION (`ParallelMixed.m` b_mono1/b_mono2,
  ~2386–2405) so the Expand outputs aren't materialised+Orderless-sorted as Exprs. NOT the read-off
  (v0.190 already did that; it didn't move A1). Large, uncertain payoff.
- **A28** (3.2s) — algorithmic: reduce `VanishOrder` (~453) call count / per-call algebraic zero-test.
  Two caching fixes FAILED (slower) — see [[project_charlwood_native_cre_build_m1_m2]].
- **A27** (1.3s) — upstream in `ipim_main` (BuildTower/splitspecials), not pinned.
Full diagnosis + lessons: memory [[project_charlwood_native_cre_build_m1_m2]],
[[feedback_profile_dotm_phase_before_c_kernel]]. Full `dsolve_corpus_tests` = multi-hour; gate Q(i)
with the `y'==(Cos x+1)/(2-Sin y)` case.

---

# (paused) 2026-09-25: Beat Maxima total (T1c→T1d) + A19 completeness

Plan file: `~/.claude/plans/per-our-last-session-purring-rainbow.md`. Decisions: **mean-time
first** (T1c → T1d → A19); **commit to beating Maxima's 12.6s total** (full native
multivariate-over-K CRE; carry A27 with T1c, A35 with T1d). A39 non-elementary → 49 is the ceiling.

- [ ] **Phase 0 — fresh profiling at v0.189** (re-confirm hot paths before writing C).
  - [x] Baseline recorded (v0.189, TIMEOUT=60, this machine): **48/50, total 23.7s kernel** vs
        Maxima 12.6s (gap 11.1s). Heavy: A28 3.46, A1 3.02, P4 1.99, A3 1.99, A2 1.69, A35 1.61,
        A27 1.41; P8 0.63 (v0.189 field-Can win). A19 declines, A39 non-elementary. NB matching
        Maxima on the tail alone → ~15s; beating 12.6s also needs the diffuse per-case native-CRE win.
  - [x] Per-case profile done (`sample`, v0.189; full write-up in
        `scratchpad/PHASE0_PROFILE.md`). **Findings that change the build:**
        (1) **A28 = `collect_symbols_in` 8781** (Orderless sort of LARGE Q(x) sums), NOT `Simplify`
            — no `simp_*` leaf. The plan's "reduce per-place Simplify" lever is WRONG for A28; the
            lever is expression SIZE (native coefficients → big sums never built/sorted).
        (2) **A27 = `collect_symbols_in` 1262** — its `Simplify` cost surfaces as the internal
            Together/Cancel Orderless sort. Carried by the same expression-size lever. ✓
        (3) **A1 (3.02s, 2nd-biggest, NOT an original target) = evaluator churn over Q(x)**
            (evaluate_step 2141, builtin_times/power) — same disease as A2/A3, fixed by the K=Q lever. Bonus.
        (4) **A2/A3 = evaluator churn over K(x)** (evaluate_step 1654/1741) — T1d. ✓
        (5) **A35 = BIGNUM** (`__gmpn_mul_basecase` 1821) — huge reduced deg-8 coeffs; native nf_elem
            keeps coeffs reduced, partial help; bit-length partly inherent.
  - [x] **Output = build order (REVISED): T1c and T1d are ONE build at two field specializations.**
        The tower coefficient ring is MULTIVARIATE over K (A28 Q(x); A1 Q(x,L1,L2); A2 Q(i)(x,T)).
        Univariate KxRf does NOT cover it — need multivariate-poly-over-K CRE (PackedRatFuncOverK)
        from the start, K=Q the fast special (fmpq_mpoly). Thread through the tower substrate so
        Padd/Pmul/Pdiv/TMul/TDiv/TowerD/Dhat never leave native form.
- [~] **M1 (was T1c) — native multivariate rational-function-over-Q CRE.**
  - [x] **M1a — CRE kernel (isolated, done).** `TowerCRE` persistent-handle layer in
        `src/poly/flint_bridge.c` (arena of `fmpz_mpoly_q` over a gens `fmpz_mpoly_ctx`; reuses
        `expr_to_mpolyq`/`fmpz_mpoly_to_expr`). Public API in `flint_bridge.h`: `tcre_new/free/
        from_expr/to_expr/add/sub/mul/div/deriv/is_zero/equal` (+ no-FLINT stubs). `deriv` = quotient
        rule via `fmpz_mpoly_derivative`. **Differential test `test_towercre` (tests/test_flint_bridge.c):
        native chain == evaluator Together/D (add/sub/mul/div/deriv, chain persistence, deriv wrt
        gen/non-gen, zero/equal, decline) — PASSES.** Leak-clean (A/B: no new leak site vs baseline);
        `make check-c99` clean; GCC build clean. NOT wired into the integrator → no behavior change,
        no version bump. NOT committed (awaiting user request).
  - [!] **M1b — BLOCKED: `.m`-phase profiling OVERTURNED the premise.** Temporary `PMtick` timers
        (since reverted) showed the tower pair arithmetic (Padd/Pmul/Pdiv — what TowerCRE accelerates)
        is ~0.005 s on A28. The real per-case bottlenecks are DIFFERENT and heterogeneous:
        - **A28 (3.2s) = `VanishOrder`** (`ParallelMixed.m` ~453, inside RealisePoints): `Sqrt[q(rho+e)]`
          local series at Root-object places × many calls. "Expand once to order 8" was SLOWER (loop
          stops at first nonzero coeff) — reverted. Lever = algebraic-series arithmetic / fewer calls.
        - **A1 (2.77s) = AnsatzSystem `ans_build` 1.57 + `ans_eqn` 0.64** (Expand/`cg.Dg0s` Dot/
          CoefficientRules over Q). **A2/A3** = same ≈0.9–1.2s + `ans_cols` ~0.12 + ~0.55 upstream (K=Q(θ)).
        - **A27 (1.33s) = ENTIRELY upstream** (BuildTower/splitspecials; ansatz+residue ~0.004s).
        TowerCRE (M1a) MAY still help A1's ans_build/eqn if the ANSATZ polynomial arithmetic (not the
        tower substrate) is recast onto it — but must be measured vs existing Expand/CoefficientRules
        (already fmpq_mpoly over Q) first. Lesson: [[feedback_profile_dotm_phase_before_c_kernel]].
  - [x] **M1a committed** (a2e71ca2) as dormant infra (no bump/tag).
  - [x] **DIAGNOSE-ALL pass done (user-requested; .m-phase timers, reverted).** Per-case fix plan:
        - **A1 (2.77s) = AnsatzSystem** `b_mono2` 1.15 + `ans_eqn` 0.71 + `b_mono1` 0.47. Plain-Q
          multivariate Expand IS already FLINT (`flint_expand_polynomial`); the cost is the
          INTERMEDIATE Expr Plus/Times (`cg.Dg0s` Dot, `2 qF qF s1 + …`, `part×quo`) built and
          Orderless-SORTED before/around Expand. **Fix: keep the ansatz column polys native (fmpq_mpoly
          / field τ-lift) across the dot+linear-combo+Expand, read CoefficientRules directly — never
          materialise the big Expr Plus.** This is where TowerCRE/native-poly is threaded (into
          AnsatzSystem, NOT the tower substrate).
        - **A2/A3 (1.6/1.9s) = same AnsatzSystem cost over K=Q(θ) (~1.0–1.2s) + ~0.55s UPSTREAM**
          (BuildTower/splitspecials/FieldData). Field poly arithmetic + an upstream pass.
        - **A28 (3.2s) = `VanishOrder`** (`ParallelMixed.m` ~453). `Sqrt[q(rho+e)]` local series depends
          ONLY on the place rho, but is recomputed for every (norm-solution × sign × hit) `uu`. **Fix:
          memoise the per-place sqrt series (cache by rho), reuse across all uu.** (NOT "expand once to
          order 8" — that was slower, reverted: the loop stops at the first nonzero coeff.)
        - **A27 (1.3s) = iPIM pre-residue/pre-ansatz stage** (`ipim_main` 1.25s; residue+ansatz ~0).
          Trig case; exact line (FactorList/ClassifyPrime/rem/bounds/Simplify) needs one more
          iPIM-internal timing pass. Lower priority until pinned.
        **Recommended order:** A28 per-place series memo (clean, safe, biggest single case) → A1
        AnsatzSystem native-poly (biggest combined, ~2.3s; transfers to A2/A3) → A2/A3 upstream +
        field → A27. Each gated: differential value-identity + Charlwood 50 solve/verify + DSolve corpus.
- [ ] **M2 (was T1d) — extend the CRE coefficients to K=Q(θ) (nf_elem) + native ∂/∂x over K.**
      Fixes A2/A3 (evaluator churn) and A35/P4 (partial). `PackedRatFuncOverK` in flint_bridge.c
      generalizing KxRf 3788–4082. Gate: differential fuzz vs Expand/Cancel/Together/CoefficientRules/
      RowReduce over random Q(θ) + DSolve corpus.
- [ ] **A19 — independence-based S′-unit group basis** (port nfunits.c pattern add_if_independent 217 /
      p_saturate 292 / cert_saturate 363 to the function-field selection 2727–2778; replace cap 2757).
      Honest-decline fallback (48/50) if too deep. Gate: hardened 33-pt + off-gate random; A16/A37/P4/A40 unregressed.

Hard gate every land: DSolve corpus green + off-gate value-identity + differential A/B + `make check-c99`
+ valgrind; substantive commits bump `src/version.h` (+0.001) and tag `v<STRING>`; changelog
`docs/spec/changelog/2026-09-21.md`.

---

# ACTIVE (2026-09-24 pm): Best Maxima mean-time + A40/A19 completeness (SUPERSEDED by 2026-09-25 above)

Plan file: `~/.claude/plans/cozy-sparking-whisper.md`. Two tracks:
**T1** native field-first tower arithmetic (the mean-time lever; hoist FieldData ahead of
the raw-radical tower substrate so `Can`/Together/Cancel/CoefficientRules hit the existing
native FLINT kernels). **T2** verify-gate hardening (2a) + A40/A19 S'-unit fix (2b).

- [x] **2a — harden the verify gate** (`ParallelMixed.m:2105-2141`). Denser FIXED set of 33
      varied-denominator positive rationals (was 11) + hoist `D[surf,x]` out of the per-point
      map (net FASTER than the old gate, which re-differentiated surf 11x). Deterministic, no
      RNG side effect (SeedRandom is global, BlockRandom absent). **VERIFIED:** 46/50 verify OK
      at independent off-gate points (branch-cut A11/A34/A35 included); P9 OK (domain-filtered);
      A39/A40 decline (expected).
      **KEY FINDING:** 2a exposed that **A19 was a silent-wrong pass** — old 11-point gate
      accepted a residual-7.6 (grossly wrong) surf (proven by stash A/B on the v0.186 binary);
      new gate declines it (cold AND warm). So prior "48/50" hid one wrong case; honest sound
      count is **47/50**. A19's antiderivative is genuinely wrong via non-deterministic S'-unit
      selection — same root cause as A40. See [[project_charlwood_a40_sunit_conjugate_degenerate]].
      TODO before commit: DSolve corpus tripwire green; version bump + tag; changelog.
- [x] **2b (A40 part) — conjugation-closed, deterministic S'-unit selection** (`ParallelMixed.m`
      ~2707-2730). CONFIRMED by probe: `NormSearchAll` DOES return both conjugate orbits for A40
      (`cvals={1/5-2I/5 ×2, 1/5+2I/5 ×2, ...}`, degb=0) — the bug was purely the orbit-blind
      `SortBy`+cap-4 filling from one orbit. Fix landed: (1) total canonical order (Re/Im tie-breaks
      → deterministic, kills batch non-determinism); (2) round-robin by norm class → cap takes one
      generator per orbit, both y-signs, matching Mathematica's 4-unit set. **A40 solves+verifies
      ~0.9s.** A1/A11/A16/A34/A37/P4 unregressed. Single-orbit cases unchanged (round-robin over one
      group = identity).
- [ ] **A19 correct solve — DEFERRED follow-up (deep).** A19 was never actually solved (silent-wrong
      pass, now honest decline). Its field is Q(i,√2); it needs `deg_x b=1` **Pell-unit** generators
      (`c=-3±2√2`) that the degb=0-first cap-4 never reaches. Round-robin does NOT help (single-orbit
      per class). Needs the **independence/saturation-based S'-unit BASIS** (a divisor-lattice
      independence test) — and there is NO Mathematica reference trace for A19 to target. Raising the
      cap over-completes → spurious wrong solve (proven). Genuine subproject. → 49/50 when done.
- [~] **T1 — mean-time lever (profiled, heterogeneous tail).** Baseline tail (9 heavy cases,
      warm) = ~20-21.7s = 77% of the ~26s total. Fresh profiling: **A28 (4.4s) ~80% `collect_symbols_in`**
      (rational-bound); field cases (P8/P4/A2/A3/A35/A40 ~11s) = distributed evaluator round-trip
      overhead. Two DISJOINT levers.
    - [x] **T1a — persistent per-node symbol-set cache** (v0.188). `symset_cache` field on Expr
          function union; lazy, reused across sorts, invalidated with hash_cache (all 6 mutation
          sites via `expr_invalidate_hash`), freed with node, reset on unshare; `symmemo_lookup`
          consults/populates it; `MATHILDA_NO_SYMSET_CACHE=1` A/B. **A28 4.44→3.15s (~29%), tail
          ~6%**, system-wide Orderless win. Differential byte-identical (8 large exprs); sort/expand/
          evaluate/eval/rootreduce/PMT suites pass; 0 leaks. Cost: sizeof(Expr) 56→64.
    - [x] **T1b — field-aware `Can`** (v0.189). Profiled: `Can` = 63% of P8, 33% of P4 (raw-radical
          `Cancel[…,Ext→Auto]`); A2/A3 = evaluator round-trip churn (evaluate_step+malloc+builtin_times,
          NOT Can — that's T1d/Phase B); A28 = collect_symbols_in (done T1a). Validated the lever: P8's
          worst Can 0.48s→0.0025s (~190×) by mapping atoms→AlgebraicNumber[θ] (memoised FieldData) +
          native field Cancel + map-back, value-identical. Landed field-aware `Can` (gated on algebraic
          atoms + LeafCount≥40; `$CanFieldEnabled` A/B). **P8 3.44→0.62s (5.5×), P4 2.97→1.98s, A40→0.48s;
          tail 20.3→16.3s (~20%).** Full 50 identical outcomes (48 sound); PMT/rootreduce/field/nf suites
          pass; DSolve 2_2_13 unchanged.
    - [ ] **T1d — A2/A3 evaluator churn** (native nf_elem AlgebraicNumber arith threaded through the
          assembly; the deeper Phase-B C work — cuts evaluate_step/malloc/builtin_times round-trips).
    - [ ] **T1c — A28 further** (residue-loop expression size / native Q(x) arithmetic) if needed.
- [ ] **T1(old) — native field-first tower arithmetic** (mean-time lever; after T2, or in parallel —
      disjoint code). Phase A: hoist field discovery (`FieldData` @2324) ahead of the tower
      substrate (`Can`@90, Padd/Pmul/TDiv/TowerD@103-221, Tower[]@117, residue loop 2461-2653),
      recast over AlgebraicNumber[θ] so the existing univariate KxRf CRE / field Expand /
      CoefficientRules / nf_elem RREF fire. Measure vs Maxima (target total ≤ ~12.5s). Phase B:
      native multivariate CRE-over-K in C only for the residual (the genuine gap).

Hard gate every land: DSolve corpus green + off-gate value-identity + `make check-c99` + valgrind;
substantive commits bump `src/version.h` (+0.001) and tag `v<STRING>`.

---

# CURRENT (2026-09-24): Charlwood 49/50 + fastest-of-four CAS

**Goal.** Solve the 6 remaining misses (A2, A3, A27, A35, A40, P8 → 49/50) AND make the suite
total+median beat Maxima (12.8s/92ms), Mathematica (16.5s/166ms), SymPy (56.6s/505ms). Plan file:
`~/.claude/plans/let-s-continue-to-implement-snappy-plum.md`. Full diagnosis in `MATHILDA_DIVERGENCES.md §E`.

Root cause (path-dependent): A2/A3/A35 = `quo` over raw radicals (`ParallelMixed.m:2367-2371`);
P8 = ToNumberField `$Failed` → Expr fallback; A27 = residue-realisation ZeroTest cost;
A40 = wrong S'-unit coeff over Q(√5). Speed lever = native `nf_elem` RowReduce (FLINT).

- [x] **Phase 0** — reproduced v0.178. **FINDINGS (revise plan):**
      - **ToNumberField SUCCEEDS on every rung of A2/A3/P8/A35/A40** (`fd$Failed=False`). The Expr fallback
        is NOT taken by any case — v0.178 precision escalation fixed compositum construction. ⇒ **Phase 4
        (deep P8 ToNumberField) is OBSOLETE**; P8 is clean-path, same `quo` bug as A2/A3.
      - A2/A3: conic-split rung 15 unk → **72/73 eq vs Mathematica 34** (2× doubling from `quo` over raw radicals).
      - P8: clean path; declines 18→4→12 eq. A40 1st rung: 49 eq, no solution (wrong coeff, Q(√5)).
      - A35/A40 split rungs (deg-~16, 9-atom compositum): field builds but assembly **too slow → budget**. Native arithmetic is the lever.
      - **A27: `NullSpace` lacks `ZeroTest` support** — `.m` call `NullSpace[rows, ZeroTest->…]` (`:969`) rejected
        as invalid Method (99× `NullSpace::method`), `ns` unevaluated → `First[ns]`=rows → 13663× `Dot::dotsh`.
        A27 has never worked. Fix: map rows to `AlgebraicNumber[θ]` + default exact `NullSpace` (needs Phase 1).
- [x] **Phase 1** — native `RowReduce` over `AlgebraicNumber` matrices (nf_elem/FLINT), gated + differential (120 matrices byte-identical); A2 8.9→6.8s; v0.179 tagged
- [x] **Phase 2** — A27: `ZeroTest` option in `NullSpace` (any predicate; exact RREF consulting it). **A27 solves+verifies ~1.4s (44/50).** dsolve unchanged; v0.180
- [x] **Phase 2b** — `ZeroTest` also in `RowReduce`; shared `matsol_*` machinery (NullSpace refactored onto it); v0.181
- [x] **Phase 3** — A2/A3/P8: `quo` common multiple = PRODUCT of distinct denominators (PolynomialLCM has no field path;
      any common multiple → equivalent system). **A2 1.4s, A3 1.7s, P8 3.7s all solve+verify → 47/50.** Field solvers
      (A1/A16/A19/A37/P4) unregressed; v0.182. **A35 now fast-fails on branch/sign gate** (was budget); **A40 still budget** (Q(√5) 1st-rung coeff).
- [ ] ~~**Phase 4** — P8 ToNumberField compositum~~ **OBSOLETE** (ToNumberField already succeeds; P8 folded into Phase 3)
- [x] **Phase 4b** — A35: branch gate generalized to 4th-root-of-unity factor (D[surf]==u·f, u∈{1,-1,I,-I}; A35 has u=-I from Sqrt[q]=I·Sqrt[-q]).
      **A35 solves+verifies ~1.3s → 48/50.** A11/A34 still solve. v0.183. (Do-loop, not nested Functions.)
- [x] **Kernel/nested-Functions** — VERIFIED already resolved in the kernel (7 shapes: curried, stored-inner,
      inner-in-list/assoc, slot, deep 3-level — all correct). Earlier gate failure was `MissingQ` UNIMPLEMENTED
      (`MissingQ[Missing[]]` stays unevaluated), not the closure. Gate uses a `Do`-loop, robust. (MissingQ is a real
      but separate builtin gap — flagged, not blocking.)
- [~] **Phase 5** — A40: DIAGNOSED, NOT fixed (left at 48/50 sound). Root cause: S'-unit search collects 4 units from 2
      same-conjugate-class sols (`c=1/5-2I/5`), missing the conjugate class `c=1/5+2I/5` MMA includes → ansatz can't span.
      **Both quick fixes are UNSOUND (verified):** (1) dedup-by-norm-class DROPS independent generators (same-norm sols differ
      by a non-constant Pell unit) → broke A19; (2) cap 4→8 over-completes → exact solve picks a spurious particular solution
      that passes the gate's 11 fixed sample points but is WRONG off-gate → broke A19 (wrong at x=1/4, right at gate pts).
      Real fix = principled S'-unit GROUP BASIS (deep, multi-session). ALSO exposed: the verify gate's fixed-11-sample-point
      weakness (silent-wrong can pass). See [[project_charlwood_a40_sunit_conjugate_degenerate]]. Running off-gate audit of the 48.
- [~] **Phase 6** — deep speed: PROFILED. Four-CAS benchmark (v0.183, warm): **48/50, total 26.1s, mean 544ms, median 117ms**
      (median already beats MMA 166 / SymPy 505; behind Maxima 92). Heavy tail A28 4.7s, P8 3.7s, A1 3.4s, P4 2.5s.
      **Bottleneck (sample of A28): `expr_compare`/`collect_symbols_in` in Orderless (Plus/Times) sorting, driven by
      `Together`/`Cancel` — specifically the tower canonicalizer `Can[e]=Cancel[Together[e],Extension->Automatic]`
      (`ParallelMixed.m:90`) run on RAW-radical column expressions before the field mapping → generic Expr path.**
      Native-nf_elem fix = restructure the tower column arithmetic to run over the field θ / nf_elem (the deep G4c item;
      the v0.177 field Cancel/Together kernels don't fire because Can uses Extension->Automatic on raw radicals, not
      AlgebraicNumber[θ]). Large, focused, multi-session; gate every step with a byte-identical differential vs the
      current Extension->Automatic result + the DSolve corpus. Native RowReduce (v0.179) already done (~24% on A2).
- [ ] **Phase 7** — regression (DSolve tripwire!), tests, docs, version bumps + tags

- [x] **Pollution hunt (user-requested)** — ASan found & FIXED a real memory bug: uninitialised-tail `expr_free`
      in `flint_qqbar_to_number_field_common` error path (malloc→calloc); SEGV under ASan, heap corruption in batch. v0.184.
      **A19 residual batch non-determinism PARTIALLY remains** (ASan-clean now, but A19 still non-deterministically wrong
      off-gate in the 50-case batch — a 2nd, ASan-blind cause: non-deterministic S'-unit selection producing a wrong `r`
      that the gate's 11 fixed sample points don't catch). Per-process A19 CORRECT (48/50 sound). Next: strengthen gate
      (denser/randomized samples) + deterministic S'-unit selection.

Hard gate every land: DSolve corpus green + differential value-identity + `make check-c99`.

---

# Plan: eliminate the ParallelMixedTower time-explosions & misses on Charlwood's 50

**Goal:** make Mathilda the *fastest and most complete* of the four CAS on the
Charlwood "parallel mixed" suite (SymPy, Maxima, Mathematica, Mathilda). Today
Mathilda solves 35/50 in 63.8 s; Mathematica solves 49/50 in 20.5 s. But on the
24 it *does* solve, Mathilda is already faster than Mathematica (median 0.10 s vs
0.20 s). So the entire gap is the **15 misses** — fix those at their roots and
Mathilda wins on both axes. (A39 is genuinely non-elementary; all four CAS fail
it — not in scope.)

Source of truth for the run: `MATHILDA_DIVERGENCES.md §E` and
`.../mixed/charlwood_mathilda.json`. Build 0.175.

---

## Triage — 15 misses → 4 root-cause groups

| Group | Cases | Symptom today | Root cause (confirmed) |
|---|---|---|---|
| **G1. `FreeQ` atomic-head bug** | **A1** (F,11.6s), **A16** (F,7.5s) | `no solution within bounds` after bound-retry ladder | `FreeQ[e, Complex]` returns True where WL returns False → `AlgAtoms` never adds `I` → the split rung is assembled over ℚ with a bare `Complex[0,1]` → **one extra independent equation** (A1: 81 vs 80) → inconsistent → ladder climbs bounds → 11–32 s decline |
| **G2. Compositum number-field layer** | **P8** (F,32s) | `no solution within bounds {1,0}` | `ToNumberField` returns `$Failed` on P8's 4-conjugate-**Root** compositum → assembly falls to the `Together`+`Solve` Expr route → over-generated / rejected solve. (Tier-2 `ToNumberField`: primitive-element search over a Root compositum.) |
| **G3. Radical branch/sign in reconstruction** | **A11** (F,0.4s), **A34** (F,0.1s) | `verification failed` | The formal radical `y` is mapped back as the principal `+√q`, but that branch doesn't match the physical radical the integrand pair `{f0,f1}` was decomposed on. **A11: global sign flip** (`D[S,x]=−f` everywhere; `−S` is exact). **A34: dropped `|·|`** (`√(1−Sin⁶)=Cos·√(…)` but true relation is `|Cos|·√(…)`; correct for Cos>0, negated for Cos<0). Abstract solve is correct; defect is purely in `y→√q` back-substitution. |
| **G4. qqbar field-arith slowness + TimeConstrained muted on the public path** | **A19,A20** (T; SymPy ~1s) primary; **P4,A2,A3,A27,A35,A37,A40** (T, degree-16 fields); **A28** (✓ but 5.6s) | runs to the 420 s OS backstop | **Speed:** the split-specials ansatz over a degree≥8 field does thousands of `AlgebraicNumber` `Plus`/`Times`, each round-tripping expr→qqbar→expr; every degree≥3 `qqbar_to_expr` (`flint_qqbar.c:564`) recomputes *all* conjugate roots and O(d²) selection-sorts them (`wl_sort_indices`:190 / `wl_cmp_qq`:168) with `qqbar_cmp_re`+`qqbar_im`(→`fmpz_poly_factor`) **per comparison**, unmemoized — so the ansatz assembly hangs before it even prints its equation count. **Preemption:** `builtin_integrate_pmt` (`integrate.c:582/591`) brackets the method in a `tc_async_defer` region that mutes SIGPROF; FLINT allocs bypass `tc_alloc_safepoint` (no `flint_set_memory_functions`) and the qqbar path never polls `tc_check_deadline`, so the budget can't fire inside a long FLINT call. |

### Evidence (this session)
- **G1** proven: `FreeQ[3,Integer]`, `FreeQ[1/2,Rational]`, `FreeQ[2.5,Real]`,
  `FreeQ[I+w,Complex]` all return **True** in Mathilda, **False** in WL; the
  `_Complex` blank form and `f[x]` head-match work. Patching `AlgAtoms` to detect
  `I` via `Cases[e,_Complex,…]` drops A1's split rung 81→80 and **solves it,
  verified, in ~5 s** (was an 11 s decline). Only **one** `.m` site uses the
  `FreeQ[.,type]` idiom (the `AlgAtoms` line itself).
- **G2**: with `I` now detected, P8 still declines — its Root-compositum rung has
  `ToNumberField`→`$Failed` (a distinct gap).
- **G3** proven by raw-surface differentiation at 40 digits (gate disabled): A11
  `D[S,x]=−f` on all of (−1,1); A34 flips sign exactly at `x=π/2`. Not a `D`/`N`
  false-negative, not a bad `RowReduce`/`LinearSolve` over ℚ(i).
- **G4** (samples + traces): A19/A20 stay in one FLINT chain
  `flint_qqbar_algebraic_number(923)→poly_to_algnum(734)→qqbar_to_expr(564)→
  wl_sort_indices(195)→qqbar_im→fmpz_poly_factor…` for the whole multi-hundred-second
  run — the hang is in **ansatz assembly** (never reaches the `ansatz: …equations`
  line). Trigger: the base+conic rungs fail → split-specials escalate to the
  **degree-8 compositum Q(i, Root[z⁴+6z²+1])** — a *needless* escalation, since the
  antiderivative needs no algebraic constants and SymPy solves it in ~1s. Preemption:
  plain `TimeConstrained[FactorInteger[…],2]` aborts at 2.01s and `Eigenvalues[2200²]`
  aborts — but the method's async-defer region mutes SIGPROF and FLINT bypasses the
  alloc safepoint, so the direct worker call declines *at* budget (45.0/120.0s) while
  the public `Method->…` path runs to the OS kill. Budget symbol is
  `System`$ParallelMixedTimeBudget` (Global/Private sets are silently ignored).

---

## Fixes (root-cause, phased)

### Phase 1 — G1: the `FreeQ` atomic-head bug (highest ROI: fixes A1, A16)
Two tiers; do **both** (safe immediate + correct root cause):

1. **`.m`-local (safe, proven).** `src/internal/mixed/ParallelMixed.m:2206`:
   replace `If[! FreeQ[e, Complex], …]` with `If[Cases[e,_Complex,{0,Infinity}]=!={}, …]`
   (or `! FreeQ[e, _Complex]`). Unblocks A1 immediately; confirm A16.
2. **C core (correct, WL-faithful, audited).** `freeq_at_level` / `builtin_freeq`
   (`src/funcprog.c:1562,1600`): a bare type-symbol pattern (`Integer`/`Real`/
   `Rational`/`Complex`/`Symbol`/`String`) must match an atom whose type-head is
   that symbol — WL's `FreeQ` examines heads, incl. of atomic literals. Currently
   the Rational/Complex function nodes are *excluded* from the head descent
   (line 1584) and true atoms have no head walk. Fix = when `expr` is such an atom,
   also test `form` against its type-head symbol.
   - **Blast-radius audit before landing:** grep every `FreeQ[.,Integer|Real|
     Rational|Complex|Symbol|String]` in C and tests; run the full unit suite +
     valgrind. Only one `.m` consumer, which *wants* the WL behavior, so `.m` risk
     is nil; the risk is C/tests relying on the buggy True.
   - Keep the `.m`-local fix even after the core fix (belt-and-suspenders; `.m`
     stays faithful to the reference `.wl`).

**Verify:** A1 and A16 return verified integrals in a few seconds; equation counts
match Mathematica (A1 split rung = 80). No regression in the unit suite.

### Phase 2 — G3: branch/sign in the `y→√q` back-substitution (fixes A11, A34)
The reconstruction (`ParallelMixed.m` ~2074–2110, where `surf` is built from the
tower solution via `back`/`Y` and `y→√q`) must apply the branch relating the
formal `y` to the physical radical.

1. **Track the sign/branch** established by `BuildTower` when it decomposed the
   integrand into `{f0,f1}` w.r.t. `y=+√q` (or a rescaled radical / conic param),
   and apply it on back-substitution. **A11's global flip** falls out directly
   (choose the branch consistent with the decomposition → `−S`).
2. **A34's varying sign** needs a branch-aware radical (`|Cos|` / piecewise), i.e.
   the reconstruction must not silently identify `√(g²·h)` with `g·√h` when `g`
   changes sign on the domain.
3. **Reconcile the internal numeric gate** (`ParallelMixed.m:2097–2110`) with
   branch-cut integrands: Mathilda's gate samples across the whole real domain and
   so over-rejects a legitimate *principal-branch* antiderivative (which is what
   Mathematica returns and what Charlwood's 3-point protocol accepts). Options:
   (a) emit the branch-consistent form so it verifies on the sampled region, and/or
   (b) restrict gate sampling to one branch / add a sign-corrected retry before
   declaring failure — **without weakening soundness** (a genuinely wrong answer,
   like A11 pre-fix, must still be caught).

**Verify:** A11, A34 return verified integrals; `D[result,x]−f` numerically ~0 at
the Charlwood check points; the gate still rejects a deliberately negated surface.

### Phase 3 — G4: honest budget + fast field arithmetic (fixes A19,A20; big cut to P4/A2/A3/A27/A35/A37/A40, A28)
Three independent fixes; **3a makes the budget honest, 3b/3c make the work fast,
3d removes the need for the work on A19/A20 entirely.**

- **3a. Make `TimeConstrained` honest inside the method (guarded FLINT allocator).**
  Install `flint_set_memory_functions` with the same guarded wrapper GMP already
  uses (`tc_install_alloc_guard`, `core.c:3969`), so `tc_alloc_safepoint` fires
  after FLINT allocations too. The safepoint is *cooperative* (runs after the alloc
  returns, not under the malloc lock), so it stays safe under the method's
  async-defer region that mutes SIGPROF — restoring a preemption point inside
  FLINT-heavy code. General fix: every FLINT-heavy builtin becomes interruptible.
  (Also make the budget robust to context or document that it is
  `System`$ParallelMixedTimeBudget`.) Converts every residual timeout from a 420 s
  OS kill into a clean in-budget decline.
- **3b. Memoize the WL root index per minimal polynomial.** `wl_root_index`
  (`flint_qqbar.c:201`) re-isolates and O(d²)-re-sorts a minpoly's roots on *every*
  `qqbar_to_expr` — i.e. after every field `Plus`/`Times` in the ansatz. Cache the
  sorted order / index keyed by the minpoly (complements the `Root→qqbar` memo of
  commit f9670b19). High ROI: the assembly does thousands of these over one field.
- **3c. Avoid the expr↔qqbar round-trip per field op.** `AlgebraicNumber`
  `Plus`/`Times` (`flint_qqbar.c:745/1303/1307`) go expr→qqbar→arith→expr, paying a
  `qqbar_to_expr` (root-index) reconstruction each op. Keep field elements in native
  `nf_elem`/qqbar form through the ansatz assembly and reconstruct the Expr once, at
  the end — the deeper structural win behind 3b. (Also a cheaper `wl_cmp_qq`: arb
  balls first, exact `qqbar_im`/factor only on a genuine tie.)
- **3d. Investigate why A19/A20's base+conic rungs fail (the *best* fix).** SymPy
  solves both in ~1 s with no algebraic field; Mathilda's conic rung fails →
  needless degree-8 split-specials escalation. The failure is almost certainly
  another assembly divergence (possibly G1-adjacent, or a distinct bounds/residue
  gap). If fixed, A19/A20 solve on the cheap transcendental path and never touch the
  qqbar hot spot — matching SymPy's ~1 s. **Diff the A19 conic-rung ansatz against
  Mathematica's** (same package) to find the extra/missing equation, as for G1.
- **Degree-16 cases (P4/A2/A3/A27/A35/A37/A40)** are genuinely heavier (real
  norm/S'-unit search over a degree-16 field — the log part of P4, the conic
  ladders of A2/A3). Same qqbar hot function; 3b/3c cut the constant, and G2's
  compositum-aware `ToNumberField` cuts the field construction. Re-measure after
  3a–3d and treat any residual as its own perf item.

**Verify:** A19, A20 solved & verified inside the budget (ideally ~1 s via 3d);
re-run the full 50 at the 300 s protocol; every previously-`T` case either solves or
declines *at* the budget (never the 420 s OS kill). Total time well under
Mathematica's 20.5 s.

### Phase 4 — regression, docs, release
- Add the 15 integrands to the ParallelMixedTower corpus test with a ratcheted
  baseline (fail on regression); differential + branch-pinned numeric check.
- `make check-c99`, full unit suite, valgrind on the changed C paths.
- Changelog `docs/spec/changelog/<Mon>.md`; bump `$VersionNumber`; tag.
- Update `MATHILDA_DIVERGENCES.md §E` (mark G1–G4 resolved / re-measured).

---

## Sequencing & why
1. **G1 first** — smallest, proven, fixes A1+A16, and clears the "spurious
   equation" class so later phases aren't debugging on top of it.
2. **G3** — self-contained `.m` reconstruction fix, no dependency on G1/G4.
3. **G4** — the qqbar/perf work; also makes the whole suite's timeouts honest.
4. **G2 (P8)** — hardest (compositum primitive-element `ToNumberField` over Root
   objects); overlaps G4's number-field cost. Do last; if it slips, 49/50 with a
   fast honest decline on P8 already beats the field.

## Risks / STOP-gates
- **G1 core fix blast radius** — audit all `FreeQ[.,type]` sites; if any Mathilda
  code depends on the buggy True, keep the core change gated behind the `.m` fix.
- **G3 gate reconciliation must not become a rubber stamp** — keep a test that a
  negated/wrong surface is still rejected.
- **G2 `ToNumberField` compositum** — scope creep; acceptable interim is a fast
  honest decline on P8 rather than a wrong/So-slow answer.
- **G4 memoization correctness** — the cached root order must be identical to the
  recomputed one (Root index k must stay canonical); pin with a test.

## Open items
- **Confirm A16** is fixed by the G1 fix alone (expected — its ℚ(√2, i) split drops
  `I`) vs also needing G2's compositum layer. (One 7.5 s run post-fix.)
- **G4/3d**: pin *why* A19/A20's conic rung fails (the escalation trigger) — the
  higher-value fix than qqbar speed. Diff the ansatz vs Mathematica.
- **G3**: decide A34's shape — branch-aware `|·|`/piecewise result vs gate
  reconciliation — by checking what Mathematica actually returns for A34.
- **A39** (non-elementary): optionally return the `{"not elementary", …}`
  certificate rather than `{"failed", "no solution within bounds"}` — cosmetic, all
  CAS fail it; out of scope for the win.

## REVIEW — implementation results (this session)

Landed and verified on the main binary:

| Change | File(s) | Effect (verified) |
|---|---|---|
| **G1** FreeQ atom-head fix | `src/funcprog.c` `freeq_at_level` | `FreeQ[·, <type-symbol>]` now matches WL (was free). Bug found & fixed in my own first cut: a bare-symbol guard so `FreeQ[5,_Symbol]` stays True. **Solves A1 (5.2s), A16 (13s), A19 (39s), A20 (0.5s), A37 (3.7s)** — the last four were 420s timeouts whose conic rung *failed* on the same undetected-`I` bug and needlessly escalated. |
| **G3** branch-resolved gate | `src/internal/mixed/ParallelMixed.m` (~2085) | Resolve `y=±Sqrt[q]` and branch-cut sign numerically vs the integrand; sound (accepts only exact ±f). **Solves A11 (0.4s), A34 (0.1s).** |
| **G4/3a** guarded FLINT allocator | `src/core.c` `tc_install_alloc_guard` + wrappers | `TimeConstrained` now preempts inside FLINT work. **P4/A2/A3/A27/A35/A40 now decline at the 45s budget instead of the 420s OS-kill** (also closes a latent SIGPROF-mid-FLINT-malloc crash). |
| **G4/3b** wl_root_index memo | `src/poly/flint_qqbar.c` | Per-minpoly sorted-roots cache; collapses the O(d²) factoring-comparator sort to once per field. RootReduce indices byte-identical (verified). Speedup TBD; targets A16/A19/A28 and the degree-16 declines. |

**Scorecard vs the 15 misses (final, all fixes):** **7 now solve** — A1 (5.2s), A11 (0.4s), A16 (**1.0s**, was 13s decline), A19 (**2.2s**, was 420s), A20 (0.5s), A34 (0.1s), A37 (3.7s). **6 decline honestly at the 45s budget** (were 420s OS-kills): P4, A2, A3, A27, A35, A40 — genuinely-heavy degree-16/algebraic-field cases. P8 declines faster (30s→9s via the compositum fix) but still needs deeper field work (quartic splitting field + √5 + i). A39 correct (non-elementary). A28 correct (4.8s).

**G4/3b speedup (qqbar root-index memo):** A16 13s→1.0s, A19 39s→2.2s (13–18×).

**G2 (ToNumberField conjugate compositum):** root cause was `in_field`'s fixed 64-bit `qqbar_express_in_field` precision — too low to resolve membership in a conjugate-root compositum, so the primitive-element search declined valid fields. Escalating-precision retry (64→256→1024→4096) fixes `ToNumberField[{r1,r2}]` (and any conjugate compositum). Helps P8 partially; P8's full solve remains open.

**Net Charlwood: ~42/50 solved (was 35), all faster; the rest decline honestly at budget instead of hanging to 420s.**

**Regression check:** full suite was 2-red — `reduce_corpus` (my FreeQ over-reach, FIXED → 170/170) and `dsolve_corpus_2_2_31` (case 3023, confirmed **pre-existing** — UNEVAL on clean main; not my regression).

**Remaining:** measure G4/3b speedup; attempt G2 (P8 Root-compositum ToNumberField); final combined suite + valgrind; changelog/version/tag.

## UPDATE 2 — the dsolve Q(i) regression and fix attempts

**What the full suite revealed:** the atom-detection fix (correctly detecting the `Complex` atom so A1/A16/A19/A20/A37 build their Q(i)/compositum field) makes DSolve route some *incidentally-complex* intermediates (e.g. 2.2.13-1219, `y'==(Cos+1)/(2-Sin y)`) through PMT's Q(i) assembly, which is **5-10x more FLINT multivariate-polynomial work than over Q** (profiled: `fmpq_mpoly_ctx_init` #1 leaf, then malloc churn + expr↔mpoly conversion + AlgebraicNumber round-trips) — DSolve[1219] 4.6s→41s, *same answer*.  Result: ~3-5 dsolve-corpus sections go PASS→UNEVAL (verify-timeouts, **no wrong answers**).  It is inherent to the correct detection, identical whether done in core `FreeQ` or `.m` AlgAtoms.

**Core `FreeQ` fix has broad blast radius** — WL-correct (`FreeQ[·,<type>]` now matches atom heads) but slows dsolve paths tuned around the old behaviour AND needed a bare-symbol guard (`FreeQ[5,_Symbol]` must stay True).  Reverted to core; detection done locally in `AlgAtoms` (`Cases[e,_Complex,…]`).

**Fix (B) — nested "try Q first, escalate to Q(i) on decline" — NOT VIABLE.**  Verified the escalation structure is sound (both-passes-Q(i): A16 0.86s), but the Q-first pass is only cheap for the single-atom case (A1); for the compositum cases (A16 ℚ(√2,i), A19/A20 degree-8, A37 ℚ(√5)) the bare-radical split over Q is *itself* a large 6-45s system that still declines, so it doesn't separate the cases.

**Fix (A) — algebraic-poly-arithmetic speedup — STARTED.**  Landed the first piece: a process-lifetime `fmpq_mpoly` **context cache** keyed by variable count (`src/poly/flint_bridge.c`, the ORD_LEX bridge sites) — the #1 profile leaf.  Verified correct (Together/Cancel/Factor/Expand/CoefficientRules/PolynomialGCD/RootReduce all right) and broadly beneficial (every multivariate poly op).  DSolve[1219] 40.9s→34.7s (~15%).  **Insufficient alone** — the rest of the cost is the expr↔FLINT round-trips and AlgebraicNumber arithmetic; bringing Q(i) to parity needs native `nf_elem` field arithmetic threaded through the assembly (a multi-session C-core project).

**Tree state now (builds clean):** G3 + G4a + G4b + always-on AlgAtoms detection + mpoly context cache.  = all 7 Charlwood solves (A1,A11,A16,A19,A20,A34,A37) + 6 honest declines + speedups, coupled to the dsolve Q(i) regression (reduced but not eliminated).  Core `FreeQ` reverted; G2 (compositum ToNumberField) reverted; A11 test updated.

**To finish:** either complete fix (A) (native `nf_elem`) to remove the regression, or land the regression-free subset (G3+G4a+G4b+mpoly-cache) and defer atom-detection.  Not yet committed.

## Measured summary (build 0.175, this session)
| Group | Cases | Fixed by | Confidence |
|---|---|---|---|
| G1 FreeQ atomic-head | A1, A16 | `.m` AlgAtoms (+core FreeQ) | **proven** (A1 solves ~5 s) |
| G2 Root-compositum ToNumberField | P8 | Tier-2 ToNumberField | root cause confirmed |
| G3 radical branch/sign | A11, A34 | back-subst branch + gate | **proven** (raw surface = ±f) |
| G4a budget honesty | all 9 T | guarded FLINT allocator | mechanism confirmed |
| G4b/c qqbar speed | A28, P4/A2/A3/A27/A35 | root-index memo + native field | hot fn confirmed by sampling |
| G4d escalation bug | A19, A20 | conic-rung ansatz diff | needless escalation confirmed |

---

## Review — v0.176 (native field arithmetic + field-coefficient Expand + clock throttle)

**Landed & verified this session (continuation):**
- **Native field-coefficient `Expand`** (`flint_bridge.c::flint_expand_polynomial_field`, wired into
  the top-level FLINT fast path in `expand.c`): the counterpart of the earlier native *scalar*
  field arithmetic, for the polynomial `Expand` that dominated the residual profile. Substitute
  θ→τ, expand over `Q[gens,τ]` in `fmpq_mpoly`, reduce mod `M(τ)` via `fmpq_mpoly_divrem`, read back
  as `AlgebraicNumber`-coefficient terms in one pass (no per-term `evaluate`). Only triggers on
  polynomials already in the `AlgebraicNumber[θ,..]` representation (the assembly's post-`/. rules`
  hot path); bare-`Complex` residues are left to the generic path so their form is preserved
  (a `Complex`→`AlgebraicNumber` rewrite broke the log-pair→ArcTan recombination in general
  integration — caught by the `Integrate[1/(1+x^2),x]` PMT test and fixed).
- **`tc_check_deadline` clock-read throttle** (`eval.c`): sample the cooperative wall-clock read once
  per 128 `evaluate()` iterations; SIGPROF stays the primary timeout. ~10% of the `Q(i)` profile.

**Measured:** `DSolve[y'==(Cos x+1)/(2-Sin y)]` (the incidental-`Q(i)` regression) 40.9 s → ~8 s
(2.4× from field-Expand alone under load; warm passes the corpus 8 s cap 11/12). Differential test
of field-Expand vs `Expand` over random `Q(i)`/`Q(√2)` polys: byte-identical values. All of
`parallelmixedtower_tests`, `rootreduce_tests`, `expand_tests`, `integrals_tests`, `make check-c99`
green; `Integrate[1/(1+x^2),x] == ArcTan[x]` restored.

**Pre-existing, NOT this work (verified against HEAD b13e4dc0):** the full `Integrate[]` cascade on
A11 (`x^3 ArcSin[x]/Sqrt[1-x^4]`) runs >200 s and returns a `Simplify`-unverifiable answer — the
same on the prior release. The direct `ParallelMixedTower` method solves and numerically verifies it
(PMT test, x=1/3). The full-cascade dispatch for this case is a separate item.

**Remaining for full margin (Phase 2b/3, not done):** native `Together`/`Cancel`/`CoefficientRules`
over field data — the residual after `Expand` — to push `DSolve[1219]` comfortably below 8 s;
Phase 3 P8 conjugate-root compositum `ToNumberField`; degree-16 cases (P4/A2/A3/A27/A35/A40).

## Review — v0.185 (Orderless canonicalisation: per-sort symbol-set memo) — SPEED

**Context:** user directive after the Maxima investigation — "port a highly efficient version of
Maxima's rational simplification to Mathilda" → chose the field-first tower restructure. Profiling
first (Phase 0) redirected the work to the true root cause.

**Finding (profiled, not guessed):** on the arithmetic-bound Charlwood cases the entire hot path is
ONE C function — `collect_symbols_in` (`sort.c`), called from `expr_compare`'s polynomial-degree
tier, which re-walks BOTH operands' whole subtrees on EVERY pairwise comparison. Canonically
ordering an n-term Orderless `Plus`/`Times` re-walks each term O(log n) times. This overturned two
wrong locus guesses: (1) NOT the ansatz `RowReduce` (A28's ansatz is trivial, `alg=False`, ~0 s —
the 4 s is the residue loop); (2) NOT `Can`/`Extension->Automatic` (already free on radical-free
input). See [[project_charlwood_bottleneck_collect_symbols_in]].

**Landed:** `expr_orderless_sort` (sort.c) with a per-sort symbol-set memo (cache each node's symbol
set by identity for one sort). Pitfalls fixed: by-value return (slot array moves on grow); a bump
**arena** for the symbol arrays (per-node malloc was a net loss on many-small-term sorts); a cheap
**size gate** (engage only when an arg has ≥48 nodes) so tiny `Plus`/`Times` pay zero overhead.
Wired into `builtin_plus`/`builtin_times`, `eval_sort_args`, eval.c Orderless sort; `sort_abort_reset()`
at the TimeConstrained siglongjmp recovery. `MATHILDA_NO_SYMMEMO=1` disables it (A/B).

**Verified:** A28 in-process A/B 4.90 s → 4.11 s (~16%); no small-term regression (microbench 1.99
vs 1.98); differential memo-vs-nomemo output byte-identical (md5); DSolve corpus 607 pass / 0 fail /
597 non-PASS ≤ 605 baseline (no regression); sort/expand/eval/evaluate/expandfrac/parallelmixedtower/
rootreduce/field_together_cancel/nullspace/nf_rowreduce all PASS; `leaks` 0; `check-c99` clean. The
`dsolve_tests` unit suite aborts at `t_m17_normalform_bessel` — PRE-EXISTING (confirmed identical on
pristine main via git-stash A/B), see [[project_dsolve_tests_m19_insuite_abort]].

**Honest assessment:** modest, clean, GENERAL win (any large Orderless sort benefits) but NOT the
7-14× Maxima gap. The transformative lever remains expression SIZE — the field-first tower
restructure (compact `AlgebraicNumber[θ]` coefficients through the residue/tower arithmetic instead
of raw Root/radical trees). Next.

## Review — v0.186 (RowReduce OneStepRowReduction FLINT Q fast path) — SPEED

**Finding (profiled):** the ParallelMixedTower ansatz calls `RowReduce[aug, Method ->
"OneStepRowReduction"]` (chosen because it is far faster on AlgebraicNumber systems). But
`rowreduce_onestep` (linsolve.c) lacked the FLINT `fmpq_mat_rref` fast path that `rowreduce_divfree`
(the default method) has — so a PURE-RATIONAL `aug` reduced by the Expr Gauss-Jordan loop. Charlwood
A1 (alg=False) has a 311×152 rational system: RowReduce was ~0.54 s of its 1.9 s biggest ansatz call.

**Landed:** added the same `flint_mat_rref` fast path to `rowreduce_onestep`, before the classical
loop. Non-rational matrices return NULL and fall through (AlgebraicNumber → nf_elem path upstream;
inexact → the invented-1/0 classical loop). RREF is unique → value-identical.

**Verified:** 200-case differential `RowReduce[m,Method->"OneStepRowReduction"] === RowReduce[m]` over
random Q matrices ALL identical; inexact `{{2.,0.},{0.,4.}}` still returns machine 1./0.; A1 3.44 s →
2.64 s (~23%); all Charlwood solve+verify; linearsolve/nf_rowreduce/nullspace PASS; DSolve corpus
611 pass / 0 fail / 593 non-PASS ≤ 605 (no regression, marginally better); check-c99 clean.

**Session speed tally (per-process, v0.184 → v0.186):** A1 3.38→2.64 (memo+FLINT-Q), A28 4.53→4.11
(memo ~16% in-process A/B), A2/A3/A19/A27/A35 ~unchanged. Remaining heavy items are algorithmic:
A35/A40 over-complex conjugate-atom field (degree-8 RootReduce back; deep — see
[[project_charlwood_a40_sunit_conjugate_degenerate]]); A3/A2 field-Expand + retry ladder. kback for
A35 was tried and REVERTED (fast back but bloats the answer to a degree-8 Root-polynomial, downstream
eats the savings).

---

## Integrate — move Goursat after ParallelMixedTower in the cascade (v0.191, 2026-09-25)

**Task:** move the `Integrate`GoursatAlgebraic` stage after `ParallelMixedTower` in the Automatic cascade.

**Plan:** (1) move the `try_goursat` call to run last (after `try_parallelmixedtower`) in the
`METHOD_AUTOMATIC` block of `src/calculus/integrate.c`; (2) verify correctness; (3) verify perf.

**What went sideways → re-plan (user chose "reorder + fix the grind"):** the literal reorder is
correct (all tests pass) but exposed a pre-existing latent grind. With Goursat no longer short-
circuiting pseudo-elliptic `F/R^p` integrands early, the *search* stages that now precede it grind
for tens of seconds before declining. Bisected with depth-1 `clock()` stage markers on
`(t-1)/((t+2)Sqrt[t^3-1])`: culprits were `try_linearity` (product-over-sum split into individually
elliptic pieces), then `integrate_derivdivides_full`'s Eliminate/Solve search, then
`ParallelMixedTower`'s tower search (45 s on the degree-4 `(t^4+2t^3-4)/(t^2 Sqrt[(t^2-1)(t^2-4)])`).

**Landed:** one shared cheap gate `has_pseudoelliptic_radical(f,x)` (a fractional-power `Power[·,
Rational]` node whose base is a degree-≥3 polynomial, via `Exponent`). Applied to: linearity's
`times_has_plus_factor` branch (skip), `try_derivdivides` (route to direct-only
`integrate_derivdivides_try`, preserving folds), and the *cascade* `try_parallelmixedtower` call
(skip; explicit `Method` untouched). Degree-≤2 radicals (`Sqrt[quadratic]`) are deliberately not
gated (their split pieces are elementary).

**Verified:** period-3, cyclotomic sibling, and the degree-4 negative-descent case all went
hang/45 s → ~0.07 s with identical answers; exponential-product control (`E^(-9x)(a E^(2x)-b E^(2x))`)
still 0.003 s (linearity gate didn't harm it). 14 integrate test suites PASS (goursat, dispatch,
derivdivides, parallelmixedtower, chebychev, linrad, quadrad, linratiorad, fresnel, jeffrey,
risch_transcendental, symmetry, beta, unknown). Clean gcc-16 build. v0.191, docs + changelog updated.

# Association overhaul (branch `association-overhaul`, started 2026-09-25)

Source: `reports/association/` (the 5-angle audit against Mathematica 15). Seven parallel streams, each in its own worktree, merged here one by one with a version bump per stream.

- [ ] S1 **Atomicity and patterns**: AtomQ/Depth/Level/FreeQ/OrderedQ; `/.` leaves keys alone; Map/Apply/Replace level specs; AssociationQ checks well-formedness; an association works as a rule set; `==` on associations.
- [ ] S2 **Mutation and write performance**:
  - fix `+=`/`++`/`--`/`*=`, `a[k]=.`, KeyDropFrom, `Delete[a,Key[k]]`, `a[k1][k2]=v`;
  - O(1) in-place updates with an incremental index;
  - drop the double canonicalisation and cut memory per entry;
  - keep RuleDelayed on rebuild; make the index and the scan agree;
  - add bench gates; fix the stale comments.
- [ ] S3 **Read-side Part and ordering**: `[[span]]`, `[[{..}]]`, `[[All]]`, `a[Key[k]]`; `Sort[a,p]`, `KeyTake` order, `Prepend`, `Catenate`, `Partition`, `KeySort[a,p]`, `Insert` with Key, `Pick` with an association mask, multi-level `GroupBy`.
- [ ] S4 **Parser, Minus and operator forms**: `#name`; `Slot["k"]` on associations; `Minus[3]`; curried forms of every association head plus `Select[p]`/`Apply[f]`; lazy `Lookup` default; `KeyMemberQ` patterns.
- [ ] S5 **New heads and threading**:
  - new heads: KeyIntersection, KeyComplement, MissingQ, JoinAcross, ApplyTo, Transpose, Discard, CountDistinct, SubsetQ;
  - extended forms: KeyUnion fill, `Keys`/`Values` with a function and over lists, AssociationMap and PositionIndex on associations, Merge of rule lists, nested rule lists, deep `Normal`;
  - arithmetic threading over associations;
  - docs for Key and Missing.
- [ ] S6 **Packed, NDArray and Compile**:
  - packed PositionIndex; vectorised GroupBy;
  - NDArray arguments to AssociationThread/AssociationMap/GroupBy;
  - packed `Values`;
  - Compile lowering of Keys, `p[k]` and Part; Map/Select as producers;
  - update the audit lists.
- [ ] S7 **Query, Dataset and JSON**: core Query operators, a minimal Dataset, RawJSON/JSON ImportString/ExportString.
- [ ] Merge all streams; full ctest against main; leaks; re-run the Mathematica differential battery (`/tmp/assocreport/diff`) and the performance benchmarks; update docs and refpages; write the changelog.
