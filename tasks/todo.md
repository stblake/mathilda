# Bronstein integration campaign (started 2026-11-12)

Living plan: src/internal/bronstein/BRONSTEIN_PLAN.md (phases, design, progress log).

- [x] 0. Bookkeeping (plan copied, this section)
- [x] 1. Spikes and substrate: fieldFn (Trager Field.m), context bridge, Risch`RischDECertified, Risch` builtins on a theta symbol
- [x] 2. Tower.m, Poly.m, Curve.m
- [x] 3. Hermite.m, LimitedIntegrate.m, AlgebraicPart.m (Special.m folded in)
- [x] 4. Puiseux.m, Special.m (exp), Residues.m
- [x] 5. Divisors.m, Specialize.m
- [x] 6. LogPart.m
- [x] 7. Pipeline.m + Method -> "Bronstein" (v0.6518)
- [x] 8. Trager + Bronstein cascade stages (v0.6519; budget profile still to rerun)
- [ ] 9 / 9b. FriCAS stress; integral generator + benchmark 256
- [ ] 10-11. Shape B
- [ ] 13. Close-out

# DSolve PDE corpus + roadmap (started 2026-11-11)

Plan: ~/.claude/plans/let-s-review-all-the-eventual-honey.md. Roadmap: DSOLVE_PDE_PLAN.md.
Corpus: DSolve_test_status/pde/ (Abbasi "PDE in CAS", 2011 records; README.md, STATUS.md).

- [x] Import: raw file + lookup.tsv (2012/2012 rows; Mma 1547, Maple 1802 match the site);
  tools/pde_corpus/import_abbasi_pde.py -> 14 DE_PDE_*.m, 2011 converted, 1 refused (pde-392,
  undefined `bc` in the source)
- [x] Verifier prelude + selftest (16/16, incl. refusing wrong transport/wave/series answers)
- [x] Runner test_dsolve_pde_corpus.c, shared ../corpus_fork.h (ODE runner refactored onto it)
- [x] Classifier tools/pde_corpus/classify_pde_corpus.m -> classes.tsv; report tool with method buckets
- [x] Wrong answers fixed: v0.6497 D[f,{x}]; v0.6498 wave/heat dropped forcing; v0.6499 fos
  Integrate-wrt-expression + PossibleZeroQ value-keyed atoms
- [x] Post-fix baseline (v0.6499): 371/2011 verified (18.4%; Mma 76.9%, Maple 89.6%), FAIL 0;
  14 ctest gates dsolve_pde_corpus_<chap>_tests (nightly; nonpass = measured + max(2, timeouts), FAIL 0)
- [x] valgrind cross-section (valgrind_dsolve_only.m, 50 records, NO TimeConstrained):
  definitely lost 17,464 B (1x) -> 19,096 B (3x), essentially flat; PDE path clean
  (dsolve_parse failure leak fixed v0.6501). Residual per-call growth is in core engines:
  - [ ] ReplaceRepeated: do_replace_all < apply_replace_repeated_nested (~13 KB/pass)
  - [ ] Together: rt_expand_exp_sums < rat_canon_build (~13 KB/pass)
  - [ ] replace_bindings rb_rec < do_replace_all (~12 KB/pass)
  - [ ] zero_test zt_exactify_reals make_rational_mpz (small)
  - note: under TimeConstrained every abort abandons its nodes by design (time-and-date.md);
    the verifier run (valgrind_sample.m) measures mostly that, not method leaks
- [x] v0.6501: TimeConstrained timeout restores Begin/BeginPackage stack + $Context/$ContextPath
- [ ] Waves W0..W9 (DSOLVE_PDE_PLAN.md §2), each to the §3 definition of done
- [x] v0.6502: fit rebases indefinite integrals (Duhamel form; no Integrate[., 0]);
  LinearConstantCoefficients no longer declines such IVPs; g_nh_active registered for TC unwind
- [x] **dsolve_tests intermittent SIGALRM in t_m26_general_forcing** -- FIXED v0.6504.
  It was not a stall: `sample` of a hung process showed 93% of samples in eval_dep_forget
  (Module-local reclamation scanning the whole dependency log, O(n^2) under one long
  dispatch frame; the log size depends on earlier tests, hence the heisenbug).  The log
  now chains each symbol's entries, and removal is unconditional, which also closes a
  use-after-free (a mark reset to 0 by a binding bump left dangling entries).
  6 parallel runs: all 400 s alarms before, 152 s pass after.
- Pre-existing ODE gate failures (identical at v0.6496 record by record, so not this session):
  2.2.1 (28, 32, 35), 2.2.16 (1505/1508/1573 UNEVAL, **1582 and 1592 FAIL** -- 1592 is the
  separable IVP y' == 2 x y (1 + y^2), y(0) == 1 returning a 0/0-at-0 branch pair), 2.2.35,
  2.2.52, 2.2.65 (43 non-PASS vs baseline 34).
- Open (found here, out of scope):
  - [ ] D[Integrate[g[K[1]], {K[1], 1, x+y}], x] -> malformed List; D of Inactive[Integrate] with K[1] dummy
  - [ ] FirstCase[e, p, default, levelspec] (4-arg) unimplemented
  - [ ] ODE corpus 2.2.1 red at v0.6496 (2.2.1-28, 2.2.1-32 UNEVAL); Chini stress case
    tests/test_dsolve_stress.c:247 fails at v0.6496 -- both pre-date this session
  - [ ] ODE prelude: implicit answers are "UNK"=PASS unverified; port the PDE implicit-function check

# Gradshteyn & Ryzhik integral review: corpus + triage + GR table (started 2026-11-09)

Plan: ~/.claude/plans/let-s-review-all-integrals-recursive-nygaard.md
Tooling: tools/gr/ (TRANSCRIBE.md contract, check_raw.wls, merge.py, verify_mathilda.py,
build_corpus.py; oracle = tools/mathar/verify.wls + instance.wls).
Buckets: Can-do / Algorithmic (ranked method families) / Table (lazy-loaded GR .m).

- [x] Tooling: transcription contract, syntax checker, merge, Mathilda NIntegrate cross-check, corpus builder
- [ ] Pilot: 2.0-2.12 (p112-120), 3.38 (p395-401), 6.51 (p707-713) -> calibrate yield/accuracy/cost
- [ ] Chapter batches: Ch2 (p112-295) -> Ch3-4 (p296-667) -> Ch5 (p668-679) -> Ch6-7 (p680-907) -> Ch17 (p1156-1189)
- [ ] Spot audit: 2% of pages re-transcribed blind, diff < 1%
- [ ] Mathilda Integrate run (S/N) -> Can-do; WRONG/CRASH/NONNUMERIC bug list fixed
- [ ] NIntegrate DISAGREE list triaged/fixed
- [ ] Triage: docs/design/gr_integrate_triage.md (algorithmic families ranked; table bucket)
- [x] Corpus snapshot committed to tools/gr/corpus/ (5636662df, 2026-11-10; was only in a session scratchpad)
- [x] GR table stage (Method -> "GRTable") -- v0.6493, LAST stage of both cascades, 3745 records.
  STATUS + NEXT STEPS: docs/design/gr_integrate_triage.md (sec. 6). Pending: condaudit re-run with the
  correction step, final corpus pass + tests/ corpora + nightly floors, ch3-7/17 triage, bugs (#132, NSum 6.666,
  Inactive leak, runner Abs check), other-Abs convention review.
  - Design (2026-11-10): engine src/internal/GRTable.m = unification table (x-free subterms -> Optional slots,
    params solved from slot equations, book conds via Refine, numeric self-check); records in
    src/internal/GRTables/{Indefinite,DefHalfLine,DefWholeLine,DefUnit,DefTrig,DefOther}.m from
    tools/gr/build_table.py (verified + S-mode not CORRECT in a baseline run with an EMPTY table).
    C: integrate.c try_grtable after CRCTable (indefinite), try_grtable_definite before Roach-early/NL.
  - [x] engine + generator + C wiring; smoke tests pass (Beta Mellin symbolic+numeric, ArcTan instance)
  - [ ] re-verify NOHEAD records (special functions now exist), finish ch34 pycheck
  - [ ] baseline Integrate runs (empty table) over GR2..GR17 corpora -> run*.jsonl
  - [ ] build table from baseline; second run with table; WRONG/CRASH = 0
  - [ ] triage doc for ch3-7, ctest gr_* nightly, docs, changelog, version bump + tag
- [ ] ctest gr_* corpora (nightly), docs, changelog, version bump + tag
- [ ] AFTER all dev campaigns: G&R integrals as a benchmarks/ gap experiment (user, 2026-11-09)
- Status 2026-11-09: transcription DONE (104 batches, ch2-7 + ch17; gap audit clean -- remaining
  holes are 7th-ed numbering: 3.248.5, 3.851.1-4, 6.290.1-2). ch2 verified: Mathilda 1686/1803,
  Python agrees on 1622 both-verified; 63 both-mismatch (book/reading errors), 40 conflicts blocked.
- Mathilda bugs found by the pipeline (fix phase; minimal repros in scratchpad/gr/bugs.m, bug4.m):
  - [ ] B1 D[Sum[x^k/k,{k,1,n}],x] -> Derivative[1,0][Sum] with leaked index k (D not threaded over Sum, symbolic limit)
  - [ ] B2 TimeConstrained cannot interrupt Sum closed-form/Simplify/FLINT mpoly gcd (G&R 2.673: >25 min past a 300 s budget)
  - [ ] B3 NIntegrate silent WRONG: Int_0^inf J_{7/4}(7x/4) J_{3/4}(7x/4) = 0.630 (true 1/(2b)=0.2857), no message
  - [ ] B4 NIntegrate SEGFAULT: Int_0^inf I_{5/6}(x/5) J_{3/8}(5x/8) K_{5/6}(x/5)
  - [ ] B5 N of infinite Sum of BesselJ stays unevaluated (no NSum fallback in N)
  - [ ] B6 N[Derivative[0,1][BesselJ][3/8,5/8]] unevaluated (order-derivative of Bessel)
  - [ ] B7 Integrate WRONG (ch2 corpus): 2.558.4d -> 0; 2.423.39 Cosh Coth^2; 2.462.3b/10a/10b instances; 2.619.1/2/11 spurious I
  - [ ] B8 Integrate -> Indeterminate on numeric instances of (a+bx+cx^2)^(3/2)/x^3 family (2.263.4, 2.267.6, 2.269.6)
  - [ ] B9 NONNUMERIC closed forms under own Assumptions (2.246c, 2.445.2, 2.451.3a, 2.452.2b, 2.584.92/96)
  - ch2 Integrate baseline (pre-fix, contended CPU): S 1027/1621, N 1314/1621; triage: CAN-DO 1020, INSTANCE-ONLY 289,
    elliptic 192, risch 50, symbolic-order 30, liouvillian 21, BUG 19
  - [x] B1,B3-B8 fixed, merged to main + tagged/pushed v0.6470-v0.6474 (2026-11-09); B2 not reproducible once harness holds values; B9 open (ArcTan/ArcTanh choice ignores Assumptions); divergent Bessel products not detected (run to timeout)
  - [ ] Special-functions campaign RUNNING (worktree agent, from v0.6475): ParabolicCylinderD(109 G&R uses), GegenbauerC(47),
        Jacobi elliptic family(~130), HankelH1/H2 + spherical Bessel(47), EllipticTheta(12), Weierstrass(~15), BetaRegularized,
        BarnesG check, AppellF1/F4 -- all surfaces (machine/MPFR real+complex, ND/packed, Compile), >= reference/scipy speed, valgrind-clean
  - [x] Special fns families 1-7 done on branch worktree-agent-a2a1e3bfda298329b (v0.6475-v0.6481), unmerged; AppellF1-F4 done v0.6485; MPFR compile campaign done v0.6482-84 (151/151 real, 145 complex); full line on branch worktree-agent-a968fb29571eade76 awaiting user merge. TODO: InverseWeierstrassP, BetaRegularized 50-digit speed, Appell F2-F4 continuation, Plot MPFR autocompile
  - [ ] RUNNING (started on top of the 7-family branch, from v0.6482): arbitrary-precision
        lowering campaign -- every special function (existing + new) must lower in Compile[] AND auto-compilation
        at WorkingPrecision > MachinePrecision (real+complex), reusing its arb/MPFR evaluator; add an MPFR mode to
        check-compile-coverage (ratchet). Found: at WP30 NIntegrate falls back to the interpreter for all new heads
        and for Erf, AiryAi, Hypergeometric1F1, EllipticE/K, MathieuC; only Bessel J/Y + elementary lower.
        Machine Compile[] + machine autocompile already work for all new heads (verified with MATHILDA_COMPILE_DIAG=1).
  - [ ] Transcription remaps (main session): MathieuCe/Se -> MathieuC/S[MathieuCharacteristicA/B]; LobachevskyL -> Clausen/PolyLog;
        LommelS1 -> 1F2 form; ThetaN -> EllipticTheta (heads with no builtin: Lommel S/V, Volterra, MacRobert E, Schlaefli stay TABLE/no-head)
  - [ ] Missing heads -> separate special-functions campaign (user): ParabolicCylinderD, HankelH1/H2,
        SphericalBesselJ/Y, GegenbauerC, JacobiSN/CN/DN, AppellF1-F4, BetaRegularized, BarnesG,
        Weierstrass P/Zeta/Sigma, Lommel s/S, theta, Lobachevsky, Volterra nu/mu, MacRobert E

# SIAM 100-Digit Challenge benchmark + holdouts + fixes (started 2026-11-07)

Plan: ~/.claude/plans/pasted-content-id-9fcd-similar-to-fluttering-bentley.md
Register: benchmarks/SIAM100.md. Scorer: benchmarks/siam100_table.py.
Target: on every case, Mathilda's correct digits >= the best of WL / scipy+mpmath /
Maxima / FriCAS, and time < 1.5x the fastest correct system; machine + wp50 tiers.

- [x] Phase 0: harness -- FriCAS column in run_all.py (.input halves, VALUE lines kept),
      siam100_common.{m,py,mac,input}, siam100_cap.lisp, siam100_table.py
- [x] C1 (exp 166) encoded in all five systems, reference by two methods
- [x] Phase 1-2: encode C2-C10, D1-D20 (exps 167-195) with references (weak: C3 40d, D8 5d, D14 12d, D15 8d, D18 7d, D20 7d)
- [x] Baseline run of 166-195, status table in SIAM100.md (benchmarks/siam100_results/baseline_v0.6377.md)
- [x] Phase 3: holdout A written + committed (e1ad3fd0)
- [x] Holdout A pre-fix run DONE: Mathilda 262/600 on target (benchmarks/SIAM100_HOLDOUT_A.md, siam100_results/holdoutA_v0.6377.md) (v0.6377 binary, running; load peaked >30 from holdout-B refs + fix builds -> DIGITS valid, TIMINGS approximate; if speed verdicts matter, re-time from a `git worktree add <dir> v0.6377` build on a quiet machine)
- [x] Phase 4: sealed holdout B, 10 per problem, branch siam100-holdout-b (worktree ../Mathilda-holdoutB, exps 226-255, 6 agents writing; WL verify after scratchpad/WL_FREE), baselines only
- [x] Phase 6 (2026-11-09): final run v0.6460. Originals 69/76, holdout A 506/600, **holdout B 395/600 (single run; now spent)**. Verdict and queue in benchmarks/SIAM100.md; details in benchmarks/SIAM100_HOLDOUT_B_RESULTS_v0.6460.md. Next queue: D17 {0,0} wrong (20), D10 FoldList speed, C6/D8 multi-D NIntegrate timeouts, D14 FEM operator forms, D2 NSum accgl, D18/D19/D20 unevaluated, C3/C5 speed. Any fix needs a holdout C to re-measure.
- [x] Post-verdict trivial fixes v0.6461-67 (pushed 2026-11-09): Element of exact constants, Lookup[a, {Key[k], ...}], GroebnerBasis Method -> "F4"/"Buchberger", Fold/FoldList compiled rows/vectors/exact folds (D10 1684 -> 60 ms), RandomChoice packed array cells, NSum step/digit-length summands (D2 bd2.1/.2/.8/.9 15/50 d)
- [ ] Follow-up queue (non-trivial; holdout B is spent, so these are measured on 166-255 as training data and need a holdout C):
  - [ ] D17: RegionBoundary and RegionDistance absent (252 guarded to "unevaluated")
  - [ ] D14/D18: 3-D and nonlinear FEM (NDSolve/NDEigenvalues on 3-D regions; nonlinear operators)
  - [ ] C6/D8: multi-dimensional NIntegrate timeouts (MPFR cubature; >3-D deterministic rules)
  - [ ] hd10.1: RandomChoice / Log element pipelines slower than numpy (machine kernels)
  - [ ] Packed arrays inside a plain List materialise when ragged/mixed (transparency-gate design; needs a redesign decision)
  - [ ] D2: bd2.7 nested digit sums 0 d; nested cases ending at 0-1 d now ~300-450 ms (original D2 1.6 -> 303 ms) -- honest but slow
  - [ ] C3/C5 speed; D19/D20 remainder (see Phase 5 open lines)
- [x] Phase 5: fix waves (W wrong / X absent / D fewer digits / S slower)
  - [x] v0.6378 NDSolve BVP by shooting (silent wrong answer on two-point conditions)
  - [x] v0.6379 WhenEvent (all integration loops; d12 12d/274ms, d12.wp50 50d, hd12 20/20), v0.6380 NDSolveValue, v0.6381 unlisted unknowns (d13.wp50 50d)
  - [x] v0.6382-86 linalg: MPFR LU honours TimeConstrained; sparse band solver (hc7.2 7.3 s -> 2 ms); exact sparse fill-column wrong answer; MPFR sparse by iterative refinement (hc7 wp50 all 50d, fastest); Norm[m,2] MPFR
  - [x] v0.6387-90 FindRoot vector eqs/MPFR systems/backtracking + get_approx_mpfr bignum Rational (D20/hd20 fixed); minimisers converge to WP (C4 wp50 15->58d); Brent residual test (d13 machine 12d); penalty scaling (D19 no longer returns start)
  - [ ] open: D19/hd19 488 s (FD gradients x 64 vars); D19 wp50 constraints at MPFR; C5/hc5 minimax (Max[] objective); hd4.9 NSum noise; hd3.3-10 unscored; hd15.7 MOL PDE at WP60; C3 Table/N build 1.5 s (speed agent)
  - [x] v0.6391-93 NSum direct summation of fast-decaying series + alternating underflow + EM far-probe hang (c10.wp50 timeout -> 0.33 ms; c8.wp50 2354 -> 1.4 ms), saturated partial sums (hc8.7), Series keeps SeriesData fields exact (hd20.9)
  - [ ] open (NSum agent): MPFR compiler lacks With[] and Sech (hc10.9 wp50 903 ms); hd20.9 wp50 Series 59 s; hd4.4 Moebius-li series ~4d; hc8.7-10 machine 1d short in FindRoot
  - [x] v0.6394-97 FoldList compiled over exact elements (d10 1677 -> 62 ms), packed RandomChoice (hd10 330 -> 13 ms), Partition overhang forms (d18 105 s -> 0.46 s), Table exact-rational kernel (c3 447 -> 70 ms); full unit tier 620/620 at v0.6397
  - [ ] open (speed agent): hc3.10 transcendental Table body 34 ms; d10.wp50 bignum loop 3.7x slower, 2d vs 3d
  - [x] v0.6398-6400 NIntegrate: oscillatory endpoint rule at WP (C1 wp50 0->50d/0.36 s; hc1 wp50 50d), 2-D point singularities (c6 7->11d, hc6 10-12d), MPFR lobes with inexact params (c9 wp50 timeout -> 48d)
  - [x] v0.6401-6407 MangoldtLambda, CubeRoot, RiemannR (+ derivative registry deriv.h, v0.6407), RotationMatrix, ArcLength/RegionMeasure/Annulus, WaveletPhi/Psi; full tier 625/625
  - [ ] open: hc1.7 Cos[Exp[1/x]] (needs u=e^(1/x)); hc6 machine 1-2d short of scipy, hc6.7 3-D; hc6 wp50 MPFR cubature; NIntegrate::accgl prints ~1e-13 for Automatic goals; hd4.7 NSum of Lambda(k)/k^s; d16 NDSolveValue heat PDE (pdeic); NDEigen*/FEM (192,193,219,222,223)
  - [x] v0.6408-11 NIntegrate parametric oscillatory correctness (C9 machine 2->10d), SQP for linear constraints + MPFR linear constraints (D19 488->22.5 s, wp50 50d), cbrt(eps) central step, minimax epigraph (C5 7.6->14d, wp50 50d; hc5 all correct)
  - [ ] open (optimiser): hc9 machine 6-9d (NIntegrate default-goal accuracy); hd19.10 wp50 18d; nonlinear constraints at WP
  - [x] v0.6412-17 NDSolve PDE variable order (d16 pdeic), exact MOL at WP (hd15.7 32d), FEM: NDEigenvalues/NDEigensystem/NDSolve on planar regions (P2, Delaunay, graded), NIntegrate over a mesh, ContourPlot over a region; d14 8d (WL 3), d18 Koch 0.024 off (WL 0.052), d17 5.2d (WL 0)
  - [ ] open (FEM): 10-digit FEM targets need P3+/Richardson; 191 Daubechies data (25-pt default grid vs t~5e-7); hd16 rough data 1-3d (grid from initial data); FEM is machine-only; hd14.6/.10 corner predicate; hd17 crossing nodal lines ~3d
  - [x] wave 2 agents (done): absent heads (RiemannR, MangoldtLambda, CubeRoot, RotationMatrix, ArcLength, Annulus, WaveletPhi/Psi); evaluator speed (Table/N, FoldList, RandomChoice pipelines, D18 polygon build)
  - note: D15 blow-up digits = default PrecisionGoal semantics (not a defect; reference worse); hd15.7 MOL PDE at WP60 returns nothing (open)
  - [ ] wave 1 agents: NDSolve events/NDSolveValue/unknowns/blow-up; NIntegrate WP60 wrong + 2-D digits + MPFR speed; FindRoot systno/non-real + optimiser MPFR stall
- [ ] Post-fix Mathilda run of 166-225 at v0.6411 (scratchpad/s100post, other systems from archived JSONs) -> next fix queue
- [x] Holdout B 25/30 families committed on siam100-holdout-b (3f7fa420), .m halves WL-verified (8 coding fixes; no truth disputed). Note: WL C4 (229) column untrustworthy (NMinimize under TimeConstrained returns box corners after earlier calls); run wolframscript output through a pipe, not `> file` (corrupts)
- [ ] Holdout B C6-C10 (231-235): author still writing; then WL check
- [ ] Holdout B baselines of WL/Py/Maxima/FriCAS (quiet machine), then Mathilda ONCE at the final version
- [x] v0.6418 #129 Linux startup crash (fortified realpath) + check-c99 rule; v0.6419-21 last digits (QK15 constants, tanh-sinh endpoint panel, FindMinimum Armijo); v0.6422 USE_MPFR=0 + LAPACK-less builds; v0.6423-31 speed (packed N, BigInt borrow, Surd index, NIntegrate quiet guard, MPFR Power, compiled WP FindMinimum, flat Plus/Times, quadratic NRoots, machine Series)
- [x] v0.6432 FIXED: NIntegrate numericalized a nested symbolic endpoint (Log[2 + Sqrt[3]]) at machine precision under WorkingPrecision -> 60: NIntegrate[1, {x, -L, L}, WP->60] - 2 N[L,60] = 2.7e-16, result labelled `60 (hd7.8 wp50 31d, pre-existing since v0.6377). Sqrt[2], Log[3] endpoints fine.
- [ ] open: C2 machine 1.7-2.5x (2x2 Det via LAPACK, Select/Tuples); D10 wp50 2.2x (interpreted Fold steps); hd20 machine 1.7-6.6x (exact Table/Series before N); N of an exact array now PACKED (reference: unpacked) -- decide
- [x] Linux CI: Outer vs Table 1-ulp on glibc = interpreter/ND Log used Re clog(x+0I) (glibc: log1p formula) vs VM log -> one helper mth_log_machine_z (v0.6434); basin_hopping seed fragility (918/1000 -> 999/1000, v0.6435); USE_MPFR=0 fm_mpfr_prog_lookup guard (719e5d10)
- [x] x^n surfaces: one double-double machine-power helper for scalar/ND/Compile (v0.6436; 100% correctly rounded; Compile scalar loop +5-30 ns per power on non-FMA x86, packed 2-5x faster); complex overflow -> Overflow[] fallback
- [x] Benchmark fairness: C6/hc6/hd15 .m machine tiers now ask PrecisionGoal 12 (PDE 10) like .py/.mac (5e6147e9); WL+Mathilda must be re-run on 171, 201, 220
- [x] c6.machine PG 12: nested QAGS v0.6433 (19.5 s -> 0.45 s)
- [x] Wave 2026-11-08 (v0.6436-v0.6458): x^n surfaces; KKT refinement for constraints at WP (hd19.10 18 -> 50 d); MPFR compiler With + Sech etc. (hc10.9 OK); Zeta' numeric + NSum MangoldtLambda Dirichlet (hd4.7 1 -> 15/50 d); FindRoot extra Newton step (hc8 +1 d); NIntegrate phase-coordinate endpoint (hc1.7 12/50 d), accgl goal, iterated WP; FEM path integrals, Dirichlet closure, NDEigenvalues Richardson, MOL auto grid + no unresolved answer; Det/Dot small, Table trig kernel (hc3.10 43 -> 2 ms); Together inexact exponents opaque (was silently exact); one numeric inexact Series path (hd20.9 wp50 56 s -> 0.3 s)
- [ ] open after wave: C2 machine ~2x Py (Select over vector tuples: compiled-predicate path is scalar-only); D10 wp50 (bignum loop compiler); hd20.1-8 machine 2-4x Py (exact Series by design); hc6 wp50 timeouts (needs 2-D MPFR method; hc6.1 48 d in 173 s); hd4.4 Moebius-li tail asymptotics; hd19.5/.7 machine 1-2 d (8-digit goals); FEM functionals 7-9 d (P3 elements); hd14.8 9 -> 7 d; hd17 crossing nodal lines; hd16.3/.9 boundary-incompatible data; wp50 MOL rows 25-point grid (explicit MPFR integrator); koch[4] mesher "boundary vertex could not be inserted"; Together[x^2.0 + 1/x] = (1.0 + x^3)/x (reference keeps x^3.) and Cancel[(x^2.0-1)/(x-1)] = 1.0 + x (reference leaves it); hc9 wp50 46-49 d + hc9.8 machine (optimiser agent running)
- [x] Holdout B 30/30 families committed (568b39a1 C6-C10), all .m halves WL-verified, no truth disputes; run_all pipes wolframscript stdout (lost 16/20 lines of 234 to a file)
- [x] Phase 6: final run, holdout B once, verdict (v0.6460; see line above)

Findings so far:
- C1 wp50: NIntegrate[Cos[Log[x]/x]/x, {x,0,1}, WorkingPrecision->60] returns -1.186 (WRONG) in 13.4 s;
  WL gets 9 digits in 0.18 s.
- C2: machine rows 1.8-2.2x slower than numpy (NSolve per bounce); wp50 best (47d).
- C3: dense SVD/Norm 3.1x slower than numpy at n=4096, 12x at n=500; wp50-n100 best.
- C4 wp50: NMinimize[..., WorkingPrecision->60] returns only machine digits (15 vs mpmath 50).
- C5 machine: 7d vs scipy 11d (NMinimize of a minimax); wp50 timeout (all).
- C6 machine: 2-D NIntegrate around corner singularity only ~7d (scipy 12d); wp50 timeout (WL 7d).
- C3 truth only 40 digits (Strebel via Trefethen page); c3.wp50 need 38.
- D16: Mathilda ABSENT WaveletPhi/DaubechiesWavelet (WL has them; others skip).
- D17: Mathilda ABSENT NDEigensystem/NDEigenvalues/ArcLength (WL only).
- D18: Mathilda ABSENT NDEigenvalues, RotationMatrix; ~105 s burned building a symbolic polygon before "unevaluated".
- D19: FindMinimum with constraints + black-box objective makes no progress (WRONG, returns ~start);
  wp50: "general constraints at WP > Machine not supported" (WL 10d, Maxima cobyla 14d machine, scipy 11d).
- D20: FindRoot::nlnum "non-real f during iteration" on an exact-rational real polynomial at wp50 (WRONG/fail; WL/Maxima/FriCAS 50d).
  machine robbins row 14d vs 15d best, 4.5x slower than scipy.
- D18 truth 7 digits (Banjai), D20 truth 7 digits consensus, D16/D17/D19 well defined (D17 exact = 2).
- D12/D13: NDSolve SILENTLY IGNORES WhenEvent (integrates to tmax, no message) -> wrong answers; D13 NDSolve::nrnum.
  NDSolveValue ABSENT. (Highest priority: silent wrong answer.)
- D15: NDSolve PDE ignores WorkingPrecision silently (machine result); machine 3.5x slower than WL (both 3d).
- D14: no FEM (DirichletCondition/NeumannValue ABSENT) -- WL only, 2-3 digits.
- D12 truth 70d, D13 68d, D14 12d consensus, D15 8d consensus.
- C7: wp50 LinearSolve on 20000x20000 SparseArray at 60 digits UNEVALUATED (WL 72 s); n=2000 wp50 97x slower than WL;
  machine 12d vs 16d (but 250x faster than WL).
- C8: NSum-in-FindRoot fine (2x slower than Maxima machine, 24x at wp50); inline NSum with symbolic t inside FindRoot HANGS;
  NDSolve cannot take a PDE (unevaluated).
- C9: NMaximize of NIntegrate: 3d vs WL 7d; wp50 FindMaximum TIMEOUT (WL 26d).
- C10: NSum sech series wp50 TIMEOUT (mpmath 1.7 ms); machine 10x slower than mpmath.
- D1: everyone wrong (digit-indicator NSum).
- D2: everyone wrong (slow series); Mathilda warns accgl but returns 60-digit wrong value.
- D3: machine 11d correct (WL wrong local min) but 29 s; wp50 FindMinimum 197 vars TIMEOUT (600 s).
- D4: RiemannR ABSENT (WL 55d builtin); FindRoot machine returns start point silently (no warning).
- D5: all optimisers stall at the nonsmooth kink (everyone 0d); wp50 Mathilda timeout.
- D6: machine 13d vs 16d (scipy/Maxima); 0.09 ms timing suspicious (warm-up cache?).
- D7: machine 11d vs 15d; wp50 6.5x slower than mpmath; hundreds of spurious Surd::int messages; CubeRoot ABSENT.
- D8: 10-D NIntegrate unevaluated ("Monte-Carlo requires constant finite real bounds"); c2R2 1.65x slower than scipy.
- D9: machine MatrixPower 465x465 3.1x slower than numpy.
- D10: FoldList with pure function over 1e6 elements 14x slower than WL; wp50 bignum loop 5x slower than Python.
- D11: FindRoot rejects vector equation f[th,T]=={4,1} (FindRoot::systno); scalar form diverges (T->58.7) -> nlnum.
- D8 truth 5d consensus (MC); D9 exact rational; D10 64d consensus; D7/D11 70d+.
Holdout A findings (dev set):
- hd20.9: Series of an inexact expression numericalises SeriesData fields (point/orders become 0.0, 5.0) -> Normal fails -> WRONG.
- hd20.1-9 wp50: FindRoot WorkingPrecision->60 "non-real f" (same as d20) -- unevaluated.
- hd19: FindMinimum + LinearSolve black box with box+volume constraints wrong on 9/10 machine rows; WP>machine equality constraints unsupported.
- hd16-18: NDSolveValue(PDE)/NDEigenvalues/NDEigensystem/Annulus/ArcLength ABSENT (WL only, 0-5 digits).
- LinearSolve badly scaled triangular: spurious "no solution"/badly conditioned warnings, 6-10 digits lost.
- hd6: wp50 NIntegrate 2.5-10x slower than mpmath; machine typically 14 vs 15 digits.
- hd7: FindRoot WP60 returns complex roots with tiny Im -> Surd::preal noise; machine 11d vs 13-15; hd7.8 wp50 31d.
- hd8: 3-D machine NIntegrate timeout (WL 8d).
- hd9: machine MatrixPower 3-5x slower than numpy.
- hd10: FoldList pure-function loop 7-12x slower than WL; Exp[Mean[Log[RandomChoice[..,10^6]]]] 47x slower.
- hc1 (196): wp50 NIntegrate WRONG on hc1.5/.9/.10 (0d, 18-55 s; WL 9/40/3); hc1.6/.8 2/6d vs WL 8-12/13-44; hc1.1-4 best.
- hc2 (197): all correct, fastest at wp50; machine 1.6-6.1x slower than numpy (per-bounce NSolve overhead).
- hc3 (198): all digits, fastest at wp50; machine Norm/SVD 5-12.5x slower than numpy (all 10 rows).
- hc4 (199): machine correct & fastest on all 10; wp50 NMinimize never refines past ~16 digits (mpmath polish 50).
- hc5 (200): NMinimize of discrete minimax only 6-8d machine (WL 13-15, scipy 11-12); wp50 timeouts.
- TOOLING: wolframscript whole-file runs of 199/200 drop CHECK lines intermittently (license contention?) -- watch in run_all.
- hd1 (206): all systems fail (digit-pattern sums); Mathilda wp50 timeouts; IntegerString[n,3] summand unevaluated in NSum.
- hd2 (207): nobody reaches 10d; Mathilda mostly matches/beats others.
- hd3 (208): machine global min found on all 10 (11-12d, others 15); FindMinimum WP60 returns only ~12d (WL 50 on 3).
- hd4 (209): RiemannR ABSENT, MangoldtLambda ABSENT; FindRoot "non-real f" again (hd4.9 machine); NSum Gram series wp50 timeout.
- hd5 (210): Norm[A, 2] UNEVALUATED on arbitrary-precision matrices (blocks every wp50 row); Frobenius-norm FindMinimum 0d (others 10-13).
- The earlier baseline kill (exit 144) was an agent's pkill -f "Mathilda -file"; baseline was resumed for 191-195, data complete.
- hd11-15 (216-220): WhenEvent ABSENT (all hd12 fail); NDSolveValue ABSENT;
  ** NDSolve two-point BVP {-u''==Exp[2x], u[0]==0, u'[1]==0} returns u(1)=-6.6e-23 SILENTLY (true 2.097) ** -- top-priority silent wrong answer;
  Bratu BVP NDSolve::underdet; NDSolve[sys, x, ...] with more unknowns than listed -> nrnum (WL accepts);
  FindRoot[f[th,T]=={a,b}, {{th,..},{T,..}}] opaque NumericQ fn -> systno (all hd11);
  blow-up endpoint 8d machine (scipy 13), WP60 only 32-33d; hd13 machine/wp excellent & fastest.
- hc6-10 (201-205):
  ** WP60 sparse LinearSolve (parallel MPFR LU) IGNORES TimeConstrained (>12 min twice in full-file runs; 60 s alone) ** -- abort-check gap in the parallel LU;
  hc6 machine 2-D NIntegrate 7-9d vs Py/Maxima 11-13; wp50 all timeout (nobody correct);
  hc7 machine sparse solve 13-4100x slower (hc7.2 5.4 s vs 1.3 ms pentadiagonal!), wp50 unevaluated (hc7.1) / timeout (5) / ~100x slower (4);
  hc8 NSum over BesselJZero non-convergent (hc8.7); wp50 NSum eigen-series 73-297x slower than mpmath; machine 2.7-13.5x slower;
  hc9 FindMaximum WP60 stalls at 32-39d (mpmath findroot 50); NMaximize∘NIntegrate wp50 timeouts;
  hc10 2-D NSum wp50 no answer (hc10.3/.8); hc10.6 wp50 212x slower than Maxima.
  Tooling: wolframscript drops later lines after NMaximize overruns TimeConstrained in a file (abort escapes vrun?) -- WL C9 column under-reports;
  201/202 halves edited after full-file run -> re-run in the holdout-A pass.

# ISSAC '97 Systems Challenge benchmark + fixes (from 2026-11-06, v0.6225)

Plan: ~/.claude/plans/let-s-review-the-mathematica-hazy-bumblebee.md
Benchmark: benchmarks/126-135 (issac97), register + status in benchmarks/ISSAC97.md (the source of truth; update its status column per fix).
- [x] Step 1: experiments 126-135, run_all bench-timeout + group H, first run both columns
- [~] Tier 1 wrong answers: W5 NSum (v0.6226); W1-W4 NIntegrate, W6 LinearSolve (agents)
- [~] Tier 2 options/args: A1-A9 (agents: FindRoot, NDSolve/FindMinimum, linalg)
- [ ] Tier 3 absences: X1 PolyLog Integrate, X2 Limit (agent), X3 LUBackSubstitution (agent), X4 GaussianQuadratureWeights, X5 FunctionInterpolation
- [~] Holdout: experiments 136-145, 100 problems (agent); every fix re-checked against its family
- [~] W17 TOP PRIORITY (user): PossibleZeroQ[Pi/3^700] / Pi/3^700 == 0 True -- zero test magnitudes in doubles; fix + audit every double-magnitude site (agent)
- [~] GROEBNER (user, 2026-11-06: "we absolutely need highly efficient Groebner bases"): state-of-the-art exact engine over Q -- F4 mod p (GBLA-style elimination), multi-modular CRT + rational reconstruction with a verified lift, DRL -> FGLM/sparse FGLM for zero-dim lex; standard families katsura/cyclic/noon/eco/reimer + dense random, DRL and lex, timed vs baseline in benchmarks/122 (agent)
- [ ] TARGET (user, 2026-11-06): all 240 holdout-A cases CORRECT and FASTER in Mathilda than the baseline (not just within 1.5x): wave 2 = wrong/incomplete, wave 3 = every SLOWER case by area
- [~] Holdout B (sealed, agent authoring; branch kept off main until the final run): once the holdout-A fixes land, a fresh 100-problem holdout written by an agent that sees neither the fixes nor holdout A; run on final main only, never fixed against in this campaign (holdout A became training data once its failures were fixed)
- [ ] Before final holdout-A scoring: re-run the baseline column for 142 and 145 (their h7 N rows and h10.7 changed after the baseline was taken)
- [x] Before push: nightly-tier check (2026-11-07: no campaign regression; Eigensystem nonsymmetric eigenvectors (pre-existing, exposed at v0.6252) fixed v0.6340; bench_pack/bench_eval rescaled to the faster calibration; dsolve stress, solve_corpus, dsolve_corpus failures identical at v0.6225) -- -- dsolve_stress / dsolve_m5_stress / dsolve_m12_stress stop at their first assertion on main (v0.6316); compare with a v0.6225 build of the same targets; fix any regression from this campaign
- [x] PUSHED 2026-11-07: main 8edb9c4d..2c780d44 and tags v0.6226-v0.6342 to origin
- [x] Wave 4 (user, 2026-11-07) done v0.6343-v0.6373; holdout C (sealed) at v0.6373: Mathilda 228/249 vs baseline 231/249, faster in 185 of 215 both-correct; open cases in benchmarks/ISSAC97.md --: holdout B gaps (B becomes training data) + sealed holdout C to re-measure overfitting
- [x] FINAL: user asked to push commits + tags once everything lands (git fetch --tags; check no remote tag collides; git push origin main --follow-tags)
- [~] Tier 4 bottlenecks: B1-B3 Series/Residue (agent), B4-B5 MPFR linalg (agent), B6 RootCount, B7 NSolve deg 1000, B8 inexact GB, B9, B10

# Roach algorithms: Integrate[..., Method -> "Roach"] (from 2026-11-05, v0.6189)

Plan: ~/.claude/plans/pasted-content-id-62fd-solving-integral-lexical-sparrow.md
Papers: Roach ISSAC'96 "Hypergeometric Function Representations", ISSAC'97 "Meijer G Function Representations", Quantum 2008 (all in ~/Documents/Research/post_phd_research/algebraic_integration/research_papers_and_books/).
Baseline v0.6189: paper integrals 7/30 correct, 2 WRONG (#25, #26); pFq values 2/11; '96 gallery 0/13; '97 gallery 0/13 closed.

- [x] Phase 0 (v0.6197): WRONG #25/#26 from LinearRadicals (get_degree_poly on unexpanded radicand); #3 abort = Limit`Scoped re-scoping recursion. Fixed + tests; ctest 600/602 (meminfo, parallelmixedtower: load flakes, pass alone)
- [x] Phase 1 (v0.6203): pFq engine (src/special_functions/hrep_*.c, companion-basis module, FLINT), seed certs, Deriv formula, FunctionExpand hook
- [x] Phase 2: done on branch feat-struvel-anger-kelvin 7241091b (agent); merge pending
- [~] Phase 3 (v0.6203): Lerch Phi + roots of unity, integral-formula generated certs, Euler/Pfaff/Kummer, unit theorems + limit, 2-param families. OPEN: PFD duplication, Struve L/Kelvin certs (after Phase 2 merge), 2F1(1,1;2/3;8/9) prints with (-1)^(1/3)
- [x] Phase 4 (v0.6203, v0.6205): Slater -> pFq engine; G-module walk incl. log case; Struve L certs + order reduction
- [x] Phase 5 (v0.6206-07): Method->"Roach" in src/internal/RoachIntegrate.m (PBM table, one/two-G, full 2.24.1 conditions, Log^k by d/ds, Lerch term, E^-x / Tan substitutions, IBP, numeric check). OPEN: Cauchy PV theorem (#20), G^{4,3}_{4,4} at 1 (#18), #9, #24; public MeijerGReduce builtin not exposed
- [x] Phase 6 (v0.6211): Integrate`Elliptic (agent), all 7 Roach elliptic integrals
- [~] Phase 7 (v0.6209, v0.6211): Roach + elliptic as late fallbacks in the Automatic cascade; corpora unchanged vs main (main has pre-existing WRONG: M6869/M6895/M8659/M8723, extreme 31/32 -- agent fixing). Default Integrate: 25/30 Roach integrals (7 before). OPEN: benchmark, RoachCorpus.m

# FeynCalc round 14: last failures + speed (from 2026-11-04, v0.6166)

Status v0.6176: per-dir 13391/13412, 2 SKIP as UNSUPPORTED (parallel subkernels; user decision), 19 left. Fixed this round: Replace whole-expression (v0.6174), #126 quotient follow-up (v0.6175), numloop seed + Collect head (v0.6176: TID-ID50, Solve3-ID2).
Remaining, by cause (bisect reports 2026-11-04; repros in scratchpad agentG7/, agentS2/, agentX/):
- A Simplify[e, sum == 1] per-subterm reduction (simp_builtins.c simp_reduce_by_equations): Prepare ID11/13/14
- B nested Collect + Factor-oriented coefficients (try_collect_per_variable; simp_split_additive n<4 gate + raw-symbol linking): Parametrize ID6/20/11, FourDivergence-ID18, PauliTrickAnticommutator-ID3 (3 u.v - D u.v -> -((-3+D) u.v)), DiracTrace-ID71 (component split, first-term rule)
- E fold monomial into power of same product (simp_power.c): Parametrize ID9
- H common base with integer-differing exponents, shared Together/Simplify helper: Parametrize ID44, ID10 (with N)
- N D[f, {x, n>=2}] keeps x-free factors outside (deriv.c higher_order_partial): Parametrize ID10
- Q quotient orientation: expanded denominator leading positive, sign into numerator sum or outside (FCDiffEqChangeVariables ID1-3)
- P partial monomial content from Collect[_, m] (AddScalingParameter-ID6): only full gcd
- M match_flat_leftover element-driven subset order (match.c ~3139): DiracOrder-ID54
- R cubic/radical conventions (solvepoly.c solve_cubic_radical; power.c bigint square part, trial bound ~65521; rationalised quadratic root): SharedToolsFRH-ID2, Factor3-ID14

Agent findings to act on (repros in scratchpad agentW/agentT/agentP):
- Matcher (W, M): Orderless/Flat first-match order differs -- reference assigns elements to slots first, bare blanks take leftovers, nested Orderless arg orders tried before moving elements; rest_. (Optional) is enumerated first because the reorder is skipped when any_opt (match.c ~1882-1940, ~2190-2520). Repros: pc[L[i4],L[i5]] hd[1,L[i5]] hd[3,L[i4]] UpValue; ReplaceList[ch[e,j] ch[x,i] Aa[i,j], r_ ch[a_,ii_] Aa[ii_,jj_] :> ii] = {j,i} ref; Replace[s[1] s[2] s[3] s[4], rest_. s[a_] :> {rest,a}] = {s[2]s[3]s[4],1} ref. Fixes TID-ID65 (WRONG answer: drops D*SPLRD via LCPC PairContract), DiracOrder-ID54, PauliChainJoin-ID19, SUNSimplify-ID83/84. Also Replace at level 0 matches a subset of a Flat head (Replace[a+b+c, a+b -> x] gives c + x; ref unchanged) -- wrong.
- T: (1) Factor's sign-rearrangement veto on product ties (simp_search.c ~782) -- reference takes lower LeafCount on SC tie for products; (2) per-coefficient factoring after Collect (try_collect_per_variable); (3) simp_pipeline_rational ties take Together's form; (4) simp_split_multiplicative lifts 1/96 + Apart (SUNSimplify-ID23).
- S: FourDivergence-ID18 Simplify collect over Sqrt[m^2 - sp[q,q]] + sign out of coefficient; Cardano form (positive radicands, real root first) for Solve/ToRadicals (SharedToolsFRH-ID2, Factor3-ID14); Sqrt[27067195991679469617841] = 104729 Sqrt[2467796558401] (square part beyond trial division); TID-ID50 0 -> 0. inside FDS (ExplicitLorentzIndex[0.]); FCLoopFromGLI-ID29#2 + FCLoopFindIntegralMappings-ID14 need $ParallelizeFeynCalc / LaunchKernels.
- Follow-ups: SeriesCoefficient[FresnelS[x], {x,0,n}] returns unflattened Times; Block[{Times}] reference still folds numbers; FortranForm narrow-width break choice; Unique[{"a","b"}] counter per element.

Plan: ~/.claude/plans/let-s-continue-our-feyncalc-wondrous-castle.md
Reference run saved (local) in packages/feyncalc/harness/ref_wl/.

- [x] P0 per-dir run at v0.6166: 13371/13412, no losses vs BASELINE_PASS (+Factor3-ID10, SharedObjectsBehavior-ID201); geomean 1.42x (1.51x at r13b; uniform ~10% drift, the same binary varies ~10% run to run -- interleave A/B for speed claims)
- [ ] W  (NEXT) DiracTrace-ID46 done v0.6168 (Factor sign on expansion); ID70 is equal (dummy labels), TID-ID65 (WRONG answer), DiracOrder-ID54, PauliChainJoin-ID19, SUNSimplify ID83/84: matcher Orderless/Optional binding order (src/match.c)
- [ ] G7 FCFeynmanParametrize ID6/9/10/11/20/44, FCFeynmanPrepare ID11/13/14 (causes A, B, E, H)
- [ ] T  form ties (agentT report 2026-11-04: four Simplify mechanisms, not one tie rule -- see review below): DiracTrace-ID71, FCDiffEqChangeVariables ID1-3, PauliTrickAnticommutator ID1/3, SUNSimplify-ID23, AddScalingParameter-ID6
- [ ] M  SUNSimplify ID83/84 (Flat+Orderless rest_ binding order)
- [x] P  FCToTeXPreviewTermOrder ID1-4: InputForm nested-paren printer v0.6169 + Block clears attributes/builtins v0.6170
- [x] F  Write2 ID1-3, 5-7: FortranForm/CForm printer + OpenWrite/ToString FormatType+PageWidth + MantissaExponent exact exponent ; v0.6167 (+6). Follow-up: reference line-break choice at narrow widths is global (breaks before a parenthesised factor), ours greedy
- [ ] S  singletons: FourDivergence-ID18, Solve3-ID2, SharedToolsFRH-ID2, Factor3-ID14, TID-ID50 (MSG), FCLoopFromGLI-ID29#2, FCLoopFindIntegralMappings-ID14 (LaunchKernels)
- [ ] SP speed: profile done (agentSP, scratchpad agentSP/agg_*.txt): interpreter overhead dominates (eval 25%, malloc 15%, matcher 12%, node new/free 10%, rb_rec 5%, Set value walks 4%, canonical order 5%). Done: lazy has_return v0.6171. Next, ranked:
  1. rb_rec/replace_bindings substitution skeleton per rule (match.c ~2640) + ov_inject_rec (options_builtin.c ~411): 5-8% corpus
  2. eval_mentions_symbols_of on every Set (eval.c ~3630): needs a per-node symbol summary
  3. allocation: expr_new_function two mallocs (node + args), real_value_cmp mpz inits: 8-12%
  4. sort.c symbol_cmp strrchr per compare -> cache short-name key; rational fast path
  5. Isolate: stored values re-evaluated when extracted via DownValues (0.5x file)
  6. apply_up_values per-arg def lookups -> per-outer-head bit; evaluator arg-classification bitmask
  7. Collect2-ID11 bulk Expand of a 100k-term sum (expand_range_first)

# MathieuS[a, q, z] + MathieuSPrime (v0.6166, 2026-11-04)

Plan: /Users/user/.claude/plans/pasted-content-id-0290-following-on-golden-cookie.md

- [x] 1. Machine kernel: odd flag, sine sums, parity at integer nu, imag_kind_s, Taylor arm at small |z|
- [x] 2. Arb kernel: odd flag, Clenshaw/Horner sine forms, imaginary-axis parity fix (MathieuCPrime bug)
- [x] 3. Bridge, builtins, exact/parity rules, sym_names, docstrings, D/Series, ND kernels, AWARE, Compile
- [x] 4. tests/test_mathieu.c: 10 MathieuS groups + MathieuC regressions
- [x] 5. Benchmark rows + numeric_sweep probe; docs, changelog, version

### Review
- The reference shares MathieuC's coefficients for S: C + I S is Floquet (checked to 30 digits).
  It is WRONG for complex a (ODE residual about 3) and at a = a_n (residual -4.5), where we return 0.
- Found and fixed in the shared plan:
  - b_{2j} gluing divided by the zero centre coefficient; replaced by fold-based glue_integer;
  - inverse iteration was singular at the snapped nu;
  - Newton's noise floor near band edges;
  - nu tie for a < 0, which flips the sign of S;
  - truncated complex-z sums at the window edge;
  - CPrime zeroed on the imaginary axis.
- Speed:
  - 10^4 digits: 0.15 s wall / 0.24 s CPU (reference 0.27 s), from two worker threads plus a top-down Newton schedule;
  - 500 digits: 1.2 ms (reference 5.7 ms);
  - se_2 over 10^6: 20 ms (scipy 765 ms);
  - generic S over 10^6: 19 ms (reference about 126 s).
- Gates:
  - c99, messages, packed-aware and compile-coverage pass;
  - ctest -LE nightly passes 589/589;
  - valgrind definitely-lost equals the startup baseline (13,440 B) for 1x and 6x workloads.

# MathieuC[a, q, z] + MathieuCPrime (from v0.6164, 2026-11-04)

Plan: /Users/user/.claude/plans/pasted-content-id-c7dd-in-src-special-f-velvet-hinton.md
Normalization (verified vs reference): generic nu -> sum c_k cos((nu+2k)z), sum c_k^2 = 1;
integer nu -> DLMF (int ce^2 = pi).

- [x] 1. Reference grid (signs, branches, near-characteristic, exact cases)
- [x] 2. Machine helpers + (a,q) coefficient cache
- [x] 3. arb kernel in flint_num_bridge (precision ladder; 500 / 10^4 digits)
- [x] 4. mathieu.c builtins, exact/parity rules, attrs, docstrings, sym_names
- [x] 5. ND kernels, AWARE, Compile, D, Series
- [x] 6. tests/test_mathieu.c + valgrind
- [x] 7. Benchmark rows (vs wolframscript, scipy)
- [x] 8. Docs, changelog, version bump, commit + tag

### Review
- Normalization: matches the reference wherever nu is real. For complex nu, the user
  chose the analytic Sum c_k^2 = 1, which leaves a z-constant factor against the
  reference (ratio-constancy tests).
- Speed:
  - 500 digits: 1.5 ms (reference 16 ms); 10^4 digits: 0.40 s (reference 0.70 s).
  - 10^6 machine points: 15 ms. ce_2 over 10^6: 17 ms, against scipy 920 ms.
- Wrong reference value: at MathieuC[-1, 2, 5 I] the reference returns ~1e45; mpmath confirms ~1e-3.
- Gates: c99, messages, packed-aware, compile-coverage, nd-surfaces and array-exactness all
  pass. Valgrind of the REPL session equals the startup baseline.
- Also fixed: the Attributes print order (NHold* before NumericFunction).

# Issue #125: DifferentialForms package (from v0.6162, 2026-11-04)

Reporter: `<< DifferentialForms.m; d[x]` gave `Dt[x]`. Plan: corpus of every exported
function + usage examples + identities, compared with recorded answers, root-cause kernel fixes.

- [x] H1 harness: tools/check_differentialforms.py + tools/differentialforms/run_df_cases.m + tests/differentialforms/corpus.m (236 cases), `make check-differentialforms`; package local in packages/differentialforms/ (untracked)
- [x] F1 Literal = HoldPattern for matching (SYM_IS_HOLDPATTERN) -- the reported bug
- [x] F2 positional ReplacePart[expr, new, pos(, npos)] (Indices -> Basis/HodgeStar/...)
- [x] F3 named operator heads interned System` (\[Wedge] reached package's Wedge rules)
- [x] F4 Dt honours Constant (wrong answer: ball Stokes 4/3 Pi R^3 == 2 Pi^2 R); Attributes[{...}]
- [x] F5 Context["name"] via $ContextPath
- [x] F6 display: SequenceForm, Format[..., InputForm|TeXForm], HoldForm kept in InputForm, held 1*x, TeX strings
- [x] F7 tests/test_package_compat.c; docs; changelog; v0.6163
- [ ] Follow-up: notebook LaTeX channel (print_latex.c) ignores FormatValues -- the notebook shows altform[1, {x}] while the text payload shows (1) dx
- [ ] Follow-up: draft #125 reply (show user before posting)

## Review (v0.6163)
- Corpus 236/236 (recorded answers cached locally; 44 identities all True; 10 TensorProduct cases are hand-derived `expect`).
- Harness trap found: a Do loop variable named `c` in the driver captured the corpus's free symbol `c` (Do binds dynamically) -- driver locals are $-prefixed.
- Valgrind: definitely-lost identical (13,440 B / 420 blocks) for Print[1], 1x and 20x of the new paths.

# DSolve corpus verifier: PossibleZeroQ only (from v0.6160, 2026-11-03)

The prelude's dsResidVerdict accepted a branch when a majority of six real sample
points were small -- a numeric acceptance path (forbidden: PZQ is the verifier).
1582's explicit branches satisfy the ODE only on Re(x/(1-Cx)) > 0; PZQ says False.

- [x] V1 prelude: dsResidVerdict = TimeConstrained PossibleZeroQ on the symbolic residual (True OK / False BAD / undecided UNK); UNFIT only for IVPs
- [x] V2 re-run all 66 corpora, old vs new verifier, diff per-record verdicts
- [x] V3 triage each new FAIL: branch-restricted answer (DSolve bug) vs PZQ false negative (PZQ bug)
- [ ] V4 fix the separable inversion (1582 class) -- keep only branches PZQ verifies, else the implicit relation
- [ ] V5 further fixes by class; re-gate; STATUS/README/PLAN/changelog; bump + tag

## Review (v0.6161)
- Verifier: PZQ, then PZQ∘PowerExpand (raw, then radicands factored) per the user's branch policy; IC at a removable singularity by two-sided Limit; UNFIT only for IVPs; harness loads distributions.m (was missing: DiracDelta[-1] unevaluated).
- PZQ fixes: ladder resolve on exact point + top-two-rung trend + unresolvable point -> re-draw (1585); distribution stage zero_test_distrib.c (1506 class); trigexp FALSE untrusted with opaque generators (1118).
- Corpus: non-PASS 996 -> 1008 (27 P->F, 14 P->U, 3 U->F, 29 U->P). 7 gates raised to measured.
- Open (V4/V5 work queue): Re/Im of symbolic roots (1292/3571/5808); Arg in Integrate system answers (1439/2707/2761/2778/4550/4575); radical inversions (1582, 5016/5017/5047/5048, 5335/5343, 1592, 3056); 14 IVPs with unfitted C at a singular IC point; Limit[Sqrt[x^2]/x, x->0, Direction->1] = 1 (wrong).

# Nasser's 2026-11-03 report (7 DSolve items), from v0.6158

Plan: ~/.claude/plans/let-s-fix-the-following-sprightly-scott.md

- [x] W1 dsolve_parse flattens nested equation/condition lists (items 1-3, item 4 IVP)
- [x] W2 Integrate`QuadraticRadicals symbolic a: opaque Sqrt[a] in Euler I (dropped Log = wrong answer), expanded-radicand degree (110 s), closed form + Collect tidy instead of the Euler decline
- [x] W3 SymmetricSquare: expand exponents, distribute mu (x^a C1 + C2 + C3 x^-a)
- [x] W4 LowerDerivativeReduction: all branches, ds_sqrt_extract_squares (shared), Cauchy repeated integration, inert fallback (items 5, 6)
- [x] W5 separable inversion Solve muted (Solve::nongen)
- [x] W6 ds_tidy_integrand on held integrands (item 7)
- [x] tests (test_dsolve 4 new, test_integrate_quadrad 1 new), docs, changelog
- [x] corpus gates (66) vs STATUS baselines; commit v0.6159 (Integrate) + v0.6160 (DSolve), tags

## Review
- ctest -LE nightly: 587/587. check-messages, check-c99 green. leaks: 0.
- DSolve corpus: 65/66 gates green; §2.2.16 red with the same pre-existing 1508/1582 records (documented in DSOLVE_PLAN, red on main before).

# DSolve: the next 1000 of Nasser's corpus, §2.2.56–§2.2.65 (problems 5501–6500), from 2026-11-03, v0.6147

Plan: ~/.claude/plans/let-s-continue-our-implementation-sequential-nebula.md
Two batches of 500: A = M84–M88 (§2.2.56–60), B = M89–M93 (§2.2.61–65).
HEAD reference: per-record TSVs of all 61 older corpora on the v0.6147 harness (scratchpad), for A/B.

- [x] A1 corpora DE_examples_2256..2260.m generated (500 scalar, 0 IVP), parse OK, check-corpus-indvar green
- [x] A2 baseline: 394/500 (88/56/86/72/92), 0 FAIL; converter fixes (a2Csc, b1x glue, \operatorname juxtaposition, nested arbfun placeholder)
- [x] A3 triage: ChangeOfVariable families (lead), polynomial-coeff special forms + reductions + Riemann-P (agent), first-order polynomial-in-p (agent)
- [x] A4 fixes: IVP local branch v0.6148, ChangeOfVariable v0.6149, ParametricQuadrature & co v0.6150, SpecialFunctionTransform / KnownSolutionReduction / SF screens v0.6151–53, PolynomialGCD + Integrate v0.6154
- [x] A5 gates, reports, README/STATUS, DSOLVE_PLAN M84–M93, changelog, lessons
- [x] B1 corpora 2261..2265 (converter: \operatorname juxtaposition, nested arbfun placeholder, Sqrt shredding); baseline 311/500
- [x] B2–B5 nonlinear reductions v0.6155, RiemannP + RationalSubstitution v0.6156, Kovacic budget v0.6157, v0.6154 regressions v0.6158; gates and docs

## Review
- Batch A 388 → 458/500 (v0.6147 base), Batch B 311 → 385/500 (v0.6149 base): 699 → 843 / 1000, 0 FAIL in 66 corpora.
- Per-record A/B v0.6158 vs v0.6147: +342; the main FAIL (2.2.1-71) fixed; 4426/4727 flutter; 5227 open (correct implicit answer, harness verify budget).
- Open: 5227; §2.2.16 gate red before the campaign (1582 verifier); TC-abort leak of the DSolve problem struct; 5869 basis simplifier.

# FeynCalc speed campaign: Mathilda >= 10x the reference kernel (from 2026-11-01, v0.6117)

Plan: ~/.claude/plans/let-s-continue-our-dev-tranquil-bentley.md
Metric: the geometric mean, over .test files, of T_ref / T_mine, for ids passing in both runs.
Both kernels run per-dir with `python3 tools/check_feyncalc.py --per-dir --compare-timing <wl dir>`
(writes packages/feyncalc/harness/TIMING.md). Quick loop: `python3 tools/feyncalc/bench_fc.py --ref-cache <saved>`.
Profiling: `MATHILDA_DISPATCH_STATS=1` (per-head dispatch, distinct exprs, epoch resets), and `sample <pid>`.

- [x] P0.1 AbsoluteTime[] sub-second (it quantised every harness timing to whole seconds) -- v0.6118
- [x] P0.2 harness microsecond timing + --compare-timing / TIMING.md; bench_fc.m/.py -- v0.6118
- [x] P0.3 reference timing run (WL 281 s wall, 100%); baseline v0.6117: corpus 0.414x (files), micro 0.30x
- [x] P1.1 matcher by-pointer def/attribute lookups (symtab_node_def) -- 15-20%
- [x] P1.3 Orderless scratch once per call -- 15-25%
- [x] P1.x Block restore/entry precise invalidation (was a global epoch reset per Block exit)
- [x] P1.x dependency-recorded fixed points + no-fire memo (Pair: 98k dispatches over 106 distinct exprs) -- micro 0.30x -> 0.66x
- [x] full ctest + corpus per-dir correctness + timing after the evaluator round (v0.6119): pass set identical
      (13214), corpus 0.414x -> 0.734x, total compute 1740 s -> 419 s (ref 272 s), wall 1540 s -> 377 s at --jobs 10
- [x] system-wide strcmp -> SYM_ sweep: 1137 sites / 113 files (v0.6120); corpus geomean 1.07x on a quiet machine
- [ ] decide: 30 Global-named literals kept as strcmp (FractionBox, SqrtBox, PositiveIntegers, CompositeQ, ...)
- [ ] P1.2 pattern predicates per stored rule (pat_is_literal / pat_has_optional ~8%)
- [ ] re-profile: what dominates now (sunF, paveReduce, contractLC are the worst micro ratios)
- [x] P2 round 1 (v0.6121): per-slot dispatch keys (Pair 36 -> 1 try/call; DummyIndexFreeQ 8.1 -> 2.8 s),
      no foreign-call group spans under Flat (Contract-ID23 9.9 -> 0.49 s, ref 0.45; MatchQ[2 x y, y_^n_. x] now False),
      FreeQ pointer walk for atom/_h/Alternatives forms;
      corpus 1.072x -> 1.157x (total 315 -> 261 s, ref 272 s), micro warm 1.00x -> 1.09x
- [x] P2 FactorList2 Check (v0.6122): prime-power trial division in Simplify (33M divisions/pass) -> exact roots;
      Together/Cancel exponent generators when exponents differ by an integer. RegularizeDivergence-ID2 6.5 s -> 6.5 ms
- [ ] Together[1/(x^e - 1) + 1/(x^e + 1)] left uncombined (reference: 2 x^e/((-1 + x^e)(1 + x^e)))
- [ ] Together pulls generator content like the reference: x^e (1 + 2 x)/(1 + x)
- [x] P2 OptionValue injection rebuilt every bound argument (fresh nodes, no stamps): path copy (v0.6123);
      DummyIndexFreeQ-ID5 2.8 -> 0.88 s (ref 1.27), corpus 1.157x -> 1.296x, total 257 -> 206 s
- [x] P2 (v0.6124) Block/Dispatch/With pure + frame-scoped Block bindings + held args out of INERT (FCI in DiracChain's condition);
      UpValue outer-head key (SPD Set upvalues); Det minor memo. Corpus 1.296x -> 1.443x, total 206 -> 188 s
- [x] Inverse: num/den kept as quotient; small exact matrices via memoised adjugate (v0.6125); corpus 1.456x, total 184 s
- [ ] AbsoluteTime/N taint frames (FCPrint timing inside scoped-condition bodies); Module not pure (ModuleNumber)
- [ ] OptionValue::nodef not emitted for an unknown option
- [ ] memoized NonCommFreeQ repeat call on a 50k-term sum ~5 ms vs 20 us (literal DownValue lookup: hash + deep compare)
- [x] FindDivergences flake = Simplify not splitting (a b)^c for assumed-positive a (v0.6126); +3 never-passing ids
- [ ] Refine[(-la^4)^s, la > 0] -> (-1)^s la^(4 s); FullSimplify same split
- [ ] OptionValue in bound args still rewritten when the template mentions OptionValue (inject into template before binding)
- [ ] P2 TdecGrc-ID3 TIMEOUT (ref 29 s); FCFeynmanProjectiveQ-ID2 1.05 s vs 0.04

# DDI Vol. 1: close the remaining open cases (2026-10-31, from v0.6100; updated 2026-11-01 at v0.6112)

Plan: ~/.claude/plans/let-s-stress-test-both-synchronous-valiant.md (supersedes the P-items below where they overlap)
Standing rules: valid for all complex x; PossibleZeroQ over Complexes is the only verifier (a PZQ false negative = fix PZQ);
no time caps on an explicit Method; one commit + tag per fix; don't edit ParallelMixed.m or rebuild ./Mathilda during a corpus run.

**Where we stand (v0.6112, tests/data/ddi_vol1_report.json):** Trager or PMT decides 664 of 803 (FriCAS 528); PMT 564, Trager 415;
0 wrong. Known conflicts: inverse_function_3 (iii) BOOK-CONFLICT, multiple_radicals_0 (vi)-inst DISAGREE (FriCAS's "NONELEM" is an unevaluated return).
**139 open:** 71 PMT declines + 68 timeouts; book says 122 elem, 5 nonelem, 12 unlabelled; FriCAS decides 54 of them (50 ELEM), times out on 83.
90 are parametric, 64 of which have a decided -inst twin. By family: nested_radicals 71, generalised_gunther 16, chebychev 10,
inverse_function 8, euler 8, multiple_radicals 6, laurent_poly_subst 5, rational_subs 5, other 10.

- [x] W1 PolynomialLCM with parameter-rational coefficients (poly.c) -- 5 "solution does not verify" (v0.6101, +5)
- [x] W2 own-radical back-substitution for every back rule (generators, Euler w, Y, nested) -- 38 "verification failed" (v0.6102, +35)
- [x] W3 shared KummerRadicals.m; Trager splits composita; PMT pieces via Trager first (v0.6103, Trager +47, PMT +17); book conflict: multiple_radicals_0 (vi)
- [x] W4 kernel + conic-in-x^k substitution pre-pass (v0.6105, +35: 23 certificates, 12 answers); PZQ real-band screen v0.6104
  - [ ] AlgebraicPoint in eulerStep for a conic with no rational point (nested 13/17; old P3b)
  - [ ] nested 11 (x^4 conic with a pending outer root) still times out
- [x] W5 parameter constant field in PMT -- conic path v0.6106 (+30), FLINT RowReduce over Q(params)(radicals) v0.6107-08,
      parametric Solve fix v0.6109, norm-search logands over Q(params) v0.6110 (+34), PZQ machine screen in doubles v0.6111,
      monomial / two-point Kummer substitutions v0.6112 (+74). 490 -> 564 PMT since v0.6105.
  - [ ] 52 parametric timeouts left: profile (m = 4 Together blow-up in the ansatz, multmeth-m1, gg_*, nested_*) -- folds into W7
  - [ ] 17 "no solution within bounds" (chebychev 5, nested 10 parametric pieces): bounds / generic ansatz over Q(params)
  - [ ] 6 "certificate withheld: transcendental constants": a generic non-elementarity proof over Q(params) (needs every bound proved generically)
  - [ ] 4 "residue class lost: no places over the prime for a root of the residue" (parametric residues)
- [ ] W6 general-genus torsion realisation via a Trager adapter (5 "needs torsion realisation (milestone iii)" + 1 partially ramified prime);
      PMT {p, pts, taus} -> Trager divisor -> torsionOrder/jacOrder/principalGenerator -> log terms
- [ ] W7 timeouts: $PMProfile flag timing each iPIM stage; group the 68 timeouts (52 parametric, 16 numeric: nested 9 (i)-(v)-inst,
      nested 3/10/11/12/14/17/19, inverse_function 2 (v), timofeev 3 (v)) by dominant stage and fix by algorithm; Trager P4c (2 timeouts)
- [ ] W8 Trager over Q(params) (47 Trager declines "radicand is not a rational function over Q"; 6 "algebraic constant in the integrand")
- [ ] W9 general curves: primitive element + integral basis for F(x, y) = 0 (38 PMT "tower construction failed", 80 Trager "nested radicals":
      roots nested over non-conics such as Sqrt[x^4 + 1] inner; inverse_function 2, chebychev 2) -- separate plan
- [ ] W10 close-out: full rerun with FriCAS; regenerate docs/notes/ALGEBRAIC_INTEGRAL_PERF.md (tools/ddi_report.py) with before/after vs v0.6087
      and v0.6100 plus the Mathilda-vs-FriCAS timing table; build and send the PDF; review section here; memory update
      (PMT parameter field, substitution pre-pass, two-point Kummer)

Carried over from the earlier stress-test list (below), still open:
- [ ] P2b answers whose intermediate radicand lies ON its cut across a real interval are numerically fragile (multrad_0-ii-inst): prefer forms off the cut
- [ ] P2c power-free radicands before choosing m
- [ ] P2d remove the internal numeric fallbacks now that PZQ is the sole verifier (Trager spotCheck/numericZeroQ, PMT branchUnit)
- [ ] P2a complex-point verification in tools/trager_fricas_stress
- [ ] P3c partial ramification (if W6 does not cover it)

Housekeeping:
- [x] ctest -LE nightly: two failures predating this work (mixed_numericalize_tests, match_exhaustive_observe) -- triage
  (2026-11-02: both were cwd-relative .m paths assuming tests/build/; fixed in the test + absolute corpus paths for all 14 CMake corpus args; 4/4 pass via ctest)
- [x] MEMORY.md is over its size limit (24.7 KB of 24.4 KB): shorten index lines
  (2026-11-02: 26.5 KB -> 22.3 KB; five tooling lines moved to memory/index_tooling_verification.md)

# DDI Vol. 1 stress test: fix every bug/deficiency found (2026-10-30, v0.6087+)

Plan: ~/.claude/plans/let-s-stress-test-both-synchronous-valiant.md
Standard: antiderivatives valid for ALL complex x; one commit + tag per fix; no PMT time caps.

- [x] P0 baseline report docs/notes/ALGEBRAIC_INTEGRAL_PERF.md/.pdf (PMT rerun, reclassify), commit v0.6087
      (PMT budgets off on explicit Method, eulerStep Return-in-With fix, harness/corpus/report)
- [x] P1a Series: fixed pad loses k orders at a k-th order pole, short result returned silently (series.c)
- [x] P1b PMT gate on positive evidence (fac local, Hold/Temporary, reclim Check, C overflow flag)
- [x] P1c FLINT qqbar zero^(-n) guard + FLINT throw hook (no kernel abort)
- [x] P2a complex-point verification in harnesses (ddi_stress; trager_fricas_stress still TODO)
- [~] P2b drop branch-cut acceptance DONE (v0.6093 PZQ gate, x in Complexes); own-radical fallback DONE v0.6096 (+32 PMT); TODO: answers whose intermediate radicand lies ON its cut across a real interval are numerically fragile (multrad_0-ii-inst) -- prefer forms off the cut
- [ ] P2c power-free radicands before choosing m
- [x] P2d PossibleZeroQ is the sole verifier (PZQ deep-cancellation fixed v0.6091, harness PZQ-only; Trager/PMT internal numeric fallbacks still to remove): triage 77 PZQ=False verified answers (24 Trager, 53 PMT), fix PZQ per
      class, then remove numeric fallbacks (Trager spotCheck/numericZeroQ, PMT branchUnit, harness numeric acceptance)
- [x] P3a u = x^k pre-pass (v0.6095; its sub-answers were all rejected by the Temporary gate until v0.6098: +5)   - [ ] P3b eulerStep AlgebraicPoint fallback   - [ ] P3c partial ramification (if needed)
- [ ] P3d PMT timeouts: $PMProfile, fix dominant stages algorithmically
- [x] "tower construction failed" (151 DDI declines, all BuildTower::radicals): 70 root-on-conic, 41 flat products, 40 other nested
  - [x] v0.6097 Euler step before rootNormalize's generator split + RationalSqrt for parametric conics: 56 towers build, +19
  - [x] v0.6099 KummerSplit / KummerIntegrate for flat radical products: 22 of 41 decided (21 + 1 certificate)
  - [ ] nested 13/17: conic without rational point -> P3b (algebraic point in eulerStep)
  - [ ] parametric pieces declining "no solution within bounds" (chebychev 5 params, nested 10 params) -> P3d bounds
  - [ ] 40 roots nested on non-conics (x^4 + 1 inner): general curve, the Trager Phase-5 plan
- [x] P4a Trager routing (rational/number-field rational done v0.6094; radical + algebraic constant needs NF constant field -> with P4b)  - [ ] P4b Trager over Q(params)  - [ ] P4c Trager timeouts
- [ ] final full rerun, before/after report, PDF
- Follow-up plan: Trager for general algebraic function fields (nested / multiple radicals)

# Trager: remaining slowdowns via FLINT` polynomial-matrix kernels (2026-10-30, v0.6084)

Plan: ~/.claude/plans/let-s-look-into-and-happy-perlis.md

- [x] P1 .m fixes: gcd idealReduce over NF, hnfModDet in jacReduce (NF), pseudo-random twoGenerators,
      sequential mod-q check in torsionOrderSpec (A + 2B multiplications)
- [x] P2 FLINT`PolynomialHermiteForm (DeterminantMultiple / Products / Denominator), DetAdjugate (Products),
      WeakPopovForm, QuotientRemainder, ExtendedGCD, FieldGCD (modular over Q(theta)) + .m dispatch
- [x] torsion screen (N delta-bar at 3 further good primes) -- closed the order-21 cases
- [x] P3 residue-field tower (residueFieldTower / factorOverField / extendFieldByRoot)
- [x] user correction: no search caps (MaxFieldDegree, crude cut-off, MaxTorsionSearch, MaxPrimes, ShiftBound,
      HermiteMaxIter -> derived); growing BSGS
- [x] P4 tools/trager_profile.m; conditional L-polynomial step NOT needed (max former-timeout time 23 s)
- [x] P5 test_flint_polymat (push), trager_torsion_stress (nightly), suite extensions, corpus torsion-heavy,
      leaks clean; FriCAS stress rerun
- [x] P6 docs (flint.md, calculus.md), changelog, lessons, version 0.6084

# Trager: reduced divisor classes and performance vs FriCAS (2026-10-30, v0.6083)

- [x] Jacobian.m: single-place frame, Popov reduction, canonical keys, BSGS, compact generators
- [x] spec-based torsion (Pinf / R0 references), O(log) second-prime check
- [x] Hensel-lifted residue powers, HNF mod det, two-generator addends, gcd idealReduce
- [x] sum-form verification; real form with atom Together and machine sign sampling
- [x] tests: test_trager_jacobian.c (79), suites updated; FriCAS stress rerun
- Open: 10 stress timeouts (residue splitting fields of degree 6-20 force primes ~150-430);
  1/((x-2)(x^3+1)^(1/3)) and x/((x-2)(x^3+1)^(2/3)) (genus 1, cubic residue field, order 21).

# Trager (simple radicals) + Kauers log heuristic as .m Integrate methods (2026-10-29)

Plan: docs/notes/TRAGER_PLAN.md (copy of ~/.claude/plans/i-want-to-use-proud-sparrow.md)

- [x] T1 docs/notes/TRAGER_PLAN.md + todo; probe Modulus support of polynomial builtins
- [x] T2 Phase 0-1 port (Common, Validate, Reduce, Genus, Basis, Arithmetic) + test_trager_phase0, test_trager_basis_arith
- [x] T3 Shift + Hermite + test_trager_shift_hermite
- [x] T4 Residues + test_trager_residues
- [x] T5 Normalize/HNF (field-generic), Divisor, PrincipalGen + suites
- [x] T6 TorsionModP (Ch 6 §2-3) + test_trager_torsion
- [x] T7 TragerLogTerms, KauersLogTerms, LogToArcTan, Reassemble, Pipeline, Surface + suites
- [x] T8 integrate.c glue: Method -> "Trager" / "TragerKauers", docstrings, test_trager_integrate; version bump
- [x] T9 FriCAS stress harness + corpus (nightly), triage, BASELINE
- [x] T10 docs/spec calculus.md, changelog, review section, lessons

### Review (v0.6082)
- Design change from the plan (recorded in docs/notes/TRAGER_PLAN.md): the y^i components are integrated separately
  ([T] Appendix A). Residue divisors are then exact two-generator ideals, so there is no branch-place hack.
- Engines:
  - Trager: Q-basis over a small splitting field, else per Galois orbit over Q[Z]/(q), with fallback to the
    Q-basis construction.
  - TragerKauers: Groebner heuristic, verified exactly by traces.
- Tests:
  - 9 unit suites, about 1,950 assertions, all pass, each under 60 s.
  - FriCAS stress over 360 integrands: Trager 180 agree-elem / 132 agree-nonelem / 3 Mathilda-only /
    0 FriCAS-only / 0 wrong / 0 disagree; Kauers 58 FriCAS-only, its baseline.
- Fixed along the way:
  - Definite integrals with an ArcTan-of-quotient antiderivative now go through the Newton-Leibniz jump guard.
  - The LinearSolve FLINT-path leak.
- Open:
  - Large torsion powers over cubic or larger residue fields (1/((x-2)(x^3+1)^(1/3)), order 21).
  - Genus >= 2 torsion searches for index 5/6 radicals: needs divisor reduction, [T] Ch 6 Sec 3.
  - .m number-field arithmetic speed (an index-6 genus-0 case takes about 25 s).

# Polynomial ideals on FLINT: GCD, LCM, GroebnerBasis, PolynomialReduce, Eliminate, CAD (2026-10-28)

Plan: ~/.claude/plans/let-s-continue-our-performance-agile-creek.md

- [x] P0 benchmarks/122-polynomial-ideals (.m/.mac/.gp; four systems), baseline
- [x] P1 PolynomialGCD / PolynomialLCM (v0.6077): one FLINT front for every arity, Modulus multivariate (both heads),
      zero arguments, E^(k u) generators; keep the LCM fold shape (fewer-term quotient)
- [x] P2 FLINT Groebner engine (v0.6078; Lex over Q zero-dim by DRL -> certified modular FGLM) (fmpz_mpoly / nmod_mpoly, Gebauer-Moeller + sugar), gb_compute dispatcher used
      by GroebnerBasis, Eliminate, Solve, NSolve, solvemod; Modulus with elimination variables
- [x] P3 PolynomialReduce on FLINT for Lex/DegLex/DegRevLex (v0.6078); Eliminate monomial-factor bug fixed
      (Modulus PolynomialReduce stays on gfp_divmod: no benchmark row, deferred)
- [x] P4 CAD (v0.6079): SOS witness exit, no bound Simplify, ball arithmetic before exact algebraic signs
      (the profile showed lifting, not projection, so Discriminant/native projection were not the lever;
      open: n-variable lift (umbrella 249 vs 21 ms), quartic declines at a degree-48 algebraic fibre)
- [x] P5 iterate on 122 (v0.6080): 37 AHEAD / 2 SLOWER / 1 CHECK-FAIL, ahead of Maxima and PARI on every case;
      the three remaining rows are explained in benchmarks/122-polynomial-ideals/README.md
- [x] P6 docs, changelog, docstrings, versions/tags (v0.6077-v0.6080), memory

Open (outside this task): Solve's parametric Root bodies keep unevaluated `1 x` coefficients
(Solve[w^3 + x^2 w^2 + x^4 - 4 w^2 == 0, w] -> Root[1 #1^3 + 1 #1^2 x^2 ...]).

Baseline wrong answers / gaps found by 122 (fix in the phase that owns them):
- Eliminate[{x^2+y+z==1, x+y^2+z==1, x+y+z^2==1}, {y,z}] drops the factor x^2 (gb_poly_strip_monomial): the root
  x = 0 (point (0,0,1)) is lost.  P3.
- PolynomialLCM[x^2-1, x-1, Modulus->5] -> (-1+x^2)(Modulus->5): the option is multiplied in.  P1.
- PolynomialLCM[p, 0] -> p (should be 0); PolynomialGCD[0, 0] -> 1 (should be 0).  P1.
- PolynomialGCD[E^(2x)-1, E^x-1] -> 1 (should be -1+E^x).  P1.
- PolynomialGCD / LCM with Modulus and several variables: unevaluated.  P1.
- CylindricalDecomposition[x^4+y^4-3xy+x<1, {x,y}] unevaluated (degree-12 projection).  P4.
- GroebnerBasis katsura-4 Lex does not finish in minutes (other systems 27-61 ms); cyclic-6 DRL > 60 s (0.24-2.9 s).  P2.

# Together / Cancel / Apart review, FLINT-native pipelines, 5-engine shootout (2026-10-28)

Plan: ~/.claude/plans/similar-to-our-recent-declarative-turtle.md. Trig -> True is LOW priority (last, optional).

- [x] P0 benchmarks/98-rational-shootout (Mathilda / reference / Maxima / SymPy / PARI), independent verification, baseline REPORT, per-stage profile
- [x] P1 seeded stress driver tests/test_rational_stress.c + correctness fixes (Apart Q-loss, RowReduce shape, head checks; threading asymmetry; ratcanon vs Extension)
- [x] P2a conversion layer (var hash, binary pow, Plus-as-polynomial, ctx cache, direct output)
- [x] P2b Together/Cancel single classifier dispatch, no duplicate detector walks
- [x] P2c Apart FLINT pipeline: Q (raw input), Q(params) via gr fmpz_mpoly_q, Gaussian/Extension via nf
- [x] P2d options: Modulus for all three, messages, Trig -> False accepted (v0.6072); (low) Trig -> True still open
- [x] P2e FORM PARITY (user decision 2026-10-28: match the reference): Together/Cancel denominators as a product of a
      gcd-free basis of the input denominator factors (not a full factorisation; a lone factor is left alone),
      integer content in front of the denominator, numerator = content x expanded primitive part; Apart terms with the
      content extracted ((-202 + 13 x)/(529 (5 + 3 x^2))).  Probe the rule in w3-style scripts first and re-baseline tests.
- [x] P3 perf from the report; guard Integrate/DSolve/Simplify corpora
- [x] P4 docstrings, spec, changelog, flint-offload.md, version/tags, lessons + memory
- [x] P5 factored form at user level by default (v0.6073); Maxima as a fifth weekly system + experiment 121 (stress sizes)
- [x] P6 gaps from 121: inner-sum squarefree and linear-denominator gcds (v0.6074); canonical certified output,
      dense one-variable sums (v0.6075)
- [x] P7 regressions vs v0.6072 found by the closing shootout (tog-mod, can-mod, tog-gauss): Modulus one variable on
      the gcd-free-basis form (no factorisation; matches the reference on 27/28 probes), numbers under Modulus,
      Gaussian single-piece and linear-pair shortcuts, huge-degree abort declined (v0.6076)

## Review (2026-10-28)

- Versions v0.6072-v0.6075. A factored rational-function engine in FLINT (`src/poly/ratfun.c`), with Modulus for all
  three heads, the factored output form at user level, and Apart over Q(params) natively.
- Shootout (`benchmarks/98-rational-shootout`): Mathilda fastest or within 1.2x on every Cancel/Apart family over Q,
  Extension, Gaussian and Modulus Together. Remaining gaps: small-input overhead (can-small, ap-param small,
  tog-telescope) and nested fractions (2x Maxima).
- Weekly job: Maxima column; experiment 121 has Mathilda ahead on 11 of 12, behind only on nested depth 32
  (0.15 vs 0.071 ms).
- Open:
  - Trig -> True.
  - The nested-fraction denominator form (a reference splits it against deeper denominators).
  - Dense univariate factors throughout `Rf`, to close the nested gap.
  - Cold sympy timing in the other group-A `.py` halves (flagged to the user).
  - Modulus with several variables: still factors the denominator mod p; the reference's normalisation there.

# Factor review, option parity, 4-engine shootout (2026-10-27)

Plan: ~/.claude/plans/let-s-do-a-thorough-floofy-peach.md. Trig -> True is LOW priority (last, optional).

- [ ] P0 benchmarks/99-factor-shootout (Mathilda / reference / PARI / Maxima [/ SymPy]), independent verification, baseline REPORT
- [ ] P1 reviews + seeded stress drivers: (a+b) FLINT uni/multi bridges, (c) pre-stages, (d) Trager Extension, (e+f) no-FLINT fallback + FactorList/SquareFree/Terms family
- [ ] P2 correctness fixes (each with a regression test first)
- [ ] P3 options: single-pass parse, Factor::optx, Modulus (0 = Q, ::modp, any position, multivariate nmod/fmpz_mod), GaussianIntegers, multivariate Extension, forwarding in FactorList/FactorSquareFree/IrreduciblePolynomialQ
- [ ] P4 performance from the shootout report (fast entry, single-pass conversion, direct result build, FLINT squarefree, Trager norm via resultant)
- [ ] P5 docstrings, spec, changelog, version/tags, lessons + memory
- [ ] (low) Trig -> True

# Symbol table: outstanding items after v0.6057 (2026-10-27)

Baseline (v0.6057, same machine; reference = local 14.0 kernel):
| | Mathilda | reference |
|---|---|---|
| #112 generation n=14 | 2.42 s | 6.52 s |
| Cases+Union | 0.39 s | 0.68 s |
| peak (MaxMemoryUsed) | 125.7 MB | 101.8 MB |
| skeleton n=14 | 0.51 s | 0.63 s |
| Module loop 1e5 / rule+Module 1e5 | 0.147 / 0.169 s | 0.208 / 0.222 s |
| 200k escaped-then-dropped locals inside ONE Do | +49.6 MB | +0.002 MB |
| 200k closures (Module + Function) inside ONE Do | +41.4 MB | +0.002 MB |
| 400k non-escaping locals (graveyard names) | +9.4 MB | ~0 |

- [x] A. Free reclaimed locals' names (v0.6058): audit found CompiledFunction.arg_names (now pinned by a params
      List), two globals a TimeConstrained unwind abandoned (restored in tc_run_guarded), the pack-diag table (copies)
- [x] B. Nested sweeps (v0.6058): evaluate() sweeps once the queue reaches 4096 / 2x survivors; 200k in-loop
      escapes +49.6 -> +0.9 MB.  548/548 unit suites; 34 suites clean under ASan with MATHILDA_TEMP_SWEEP_EAGER=1
- [x] C. #112 peak memory (v0.6059): result held 5.4 nodes/term (private Times head symbol + coefficient per term);
      one canonical node per SymbolDef + shared small ints [-256,1024]: peak 120 -> 94 MB (ref 102), gen 2.50 -> 2.29 s
- [ ] dsolve_tests t_system_varcoeff flaked once under ctest -j10 (PossibleZeroQ got List[True, True]: DSolve answer shape
      changed under load; 40/40 standalone OK). Not investigated -- timing/state dependent
- [x] Head[CompiledFunction[...]] unevaluated, [[0]] and _CompiledFunction too: expr_head + blank_head_matches ; v0.6060
- [x] Tests (test_modular.c reclamation cases, test_expr_sharing canonical nodes), docs, v0.6058 + v0.6059 tags

# GitHub issues #112 / #113 / #114 (2026-10-27)

- [x] #113 Together on a denominator-free sum: expand (and re-evaluate) the combined numerator when lcm_den = 1 ;
      Cancel's embedded-rational gate no longer walks into exponents (y^(1/(a-1)) is one generator) ; v0.6050
- [x] #114 Factor over Q(I): factor_degree_one concludes irreducible when gcd(L, C) = 1 (indeterminate vars only) ;
      Gaussian FLINT zero test (I -> tau, deg_tau <= 1) ; 63 s -> 26 ms ; v0.6051
- [x] #112 rb_rec shares unchanged subtrees ; Expand estimates once (child bounded) + hash-screened kernel table ;
      Union hash-dedupe then sort ; Cases `_h` head test ; Echo + EchoTiming ; v0.6052
      n=14: generation 24.1 -> 18.8 s, Cases+Union 3.62 -> 0.88 s, MaxMemoryUsed 462 -> 224 MB (same machine, A/B)
- [x] #112 round 2 (v0.6053): INERT fixed points (symbolic stamps survive unrelated bindings), Plus/Times
      self-certified fixed points, adaptive run-merge Orderless sort, Times opaque fast path, Expand no-op for
      expanded input, empty-env replace_bindings, Cases `_h` match, intern pointer cache, Block restore fix.
      n=14 generation 18.5 -> 4.5 s, Cases+Union 0.84 -> 0.33 s (same machine A/B)
- [x] #112 follow-up: FLINT expansion with opaque kernels as generators (v0.6054; gated on real multiplication)
- [x] #112 per-call overhead (v0.6055): Sum enumerates ranges <= 10^6 (closed form on a symbolic index was WRONG for
      structural bodies), probe reused as first term (body evaluated once per index), integer index step,
      symtab stored hashes + doubling table. Generation 4.7 -> ~3.2 s.
- [x] Reclaim Module locals (v0.6056): +111,005 -> +9 symbols per trace, 548 -> 50 B per Module call
- [x] Rest of docs/notes/SYMBOL_TABLE_IMPROVEMENTS.md (v0.6057): ordered iteration, single copy per changed path + stack buffers, single interning per local; symbol ids declined on measurements
- [ ] Factor over Q(I), univariate: Factor[x^2 + 2 I x - 1] should be (I + x)^2 (pre-existing gap)
- [ ] Cancel numerator form: Mathilda expands the quotient numerator ((1 + a) y -> y + a y); the reference keeps it

# Test-suite drift: 42 failing registered tests (2026-10-26)

Inventory (ctest on main at -O2; assert binaries re-run with NDEBUG so every FAIL prints): ~37 individual
failures behind 42 tests. Root cause of the drift: CI only compiles; nothing runs the suite. Rule for every
item: confirm the failure, then decide CODE BUG vs STALE EXPECTATION from first principles (derivation +
numeric check), never from another system's output; a stale expectation is fixed only after the new form
is shown equal (or the new behaviour shown right).

## W1 Wrong answers (first)
- [x] integer_assume_simplify: unsound for real n -- three root causes (prov_pos real exponent, Reduce complex
      inequalities, mod-2 even-coefficient early return) ; v0.6027
- [x] sf_machine_error: Erf of a huge complex argument dropped its imaginary part; erf = 1 - erfc via arb ; v0.6029
- [x] numeric: N[expr, p] adapts its working precision until p digits settle ($MaxExtraPrecision, N::meprec) ; v0.6031
- [x] fit: NumericSpec.contagion was stack garbage at 30+ sites; moved to file-static state ; v0.6028
## W2 Capability regressions
- [x] integrate_diffunderint: the inner integral is now ConditionalExpression[1/(1+a), Re[1+a] > 0]; the stage carries the
      condition (plain under a > -1, decided by Refine) -- 3 expectations were missing the convergence condition
- [ ] DSolve: 31 corpus cases over baseline in 17 sections (per-case diff vs ce03e99d 2026-10-18), dsolve_tests
      t_m37_abel_air, dsolve_m61 (Bessel family), dsolve_m34 (cubic series), dsolve_m63 + corpus 2.1.2 timeouts
## W3 Behaviour changes to judge
- [x] packed_list (5): stale -- explicit ToNDArray/ToPackedArray widens by design (c47229b8, 2026-09-03); automatic
      packing stays lossless (now asserted); Apart print form
- [x] ndarray_linalg (4): Chop's exact 0 is its documented semantics (92453b91); Eigenvalues[NDArray] -> NDArray (convention);
      REAL BUG fixed: Equal/Unequal left a visible NDArray vs List unevaluated
- [ ] plot3d (4): polygon counts (18 vs 9) for PlotPoints -> 4, MaxRecursion -> 0
- [x] expint_machine: numeric value verified by quadrature (v0.477 mixed-argument numericalize); risch_transcendental (2):
      inexact-base antiderivative verified by differentiation, inexact-out as the test intended
- [x] flint_bridge: Factor print form under canonical order
## W4 Print-form drift (verify equality numerically, then update)
- [x] normalize, legendreq (2), stats, simp (2), trigexpand, ml_dist: forms proven equal (symbolic + numeric); LegendreQ[1/2, 5.]
      verified against Legendre's ODE and 30 digits; FullForm now prints every digit of a machine real
## W5 Timeouts (profile, as NMinimize)
- [ ] simp_sf_airy (Struve, known 22 s), mathar_regression, mathar_definite_corpus, dsolve_m63, corpus 2.1.2
## Found along the way (not test drift; separate work)
- [ ] NOT IMPLEMENTED (unevaluated even on exact matrices): Eigensystem, MatrixExp, SingularValueList,
      CholeskyDecomposition, MatrixFunction, HessenbergDecomposition
- [ ] Simplify[Re[1 + a] > 0, a > -1] stays unevaluated while Refine gives True (Simplify should subsume Refine)
- [ ] Refine[ConditionalExpression[v, 1 + a > 0], a > -1] does not reduce to v; ConditionalExpression[v, c] - y does
      not combine to ConditionalExpression[v - y, c]
- [x] Integrate[E^Sin[x], x] spent the full 10 s ParallelMixedSpecial budget: not a missing certificate but a
      wrong tower (E^(r Tan[x]), r free -- trigStep's parity tests read a transcendental r as a coefficient) and
      ClassifyPrime on the free r never returning.  Both fixed (v0.6038): 10.2 s -> 0.66 s, dsolve_m34 passes
- [ ] Integrate[(1 + Cos[4] Csc[2]^2)/(1 + Cos[2 + Y] Csc[2] Csc[Y]), Y] at HEAD ignores TimeConstrained[_, 20]
      (> 10 min); with this batch it aborts at 20 s -- find which stage was not unwinding
- [ ] A rule whose RHS has a scoping construct binding a pattern variable is alpha-renamed:
      Function[{x}, b] /. Function[v_, w_] :> Function[v, w] gives Function[v$1, ..] (the danger set holds the
      pattern names for the Module-shadow case, A11).  Decide the intended substitution semantics.
- [ ] The message funnel has the General::stop template but does not cap repeats: a DSolve IVP fit prints
      Solve::ifun four times
- [ ] DSolve[{y'' + y'^2 + y' == 0, y[0] == 0, y'[0] == 1}] is -I Pi + Log[-2 + E^-x] (correct, principal
      branch) where Log[2 - E^-x] is the natural form: the general solution carries Log[E^-x - C[1]]
- [ ] Integrate`RischTranscendental[6^x/(6^x + 1), x] returned an UNVERIFIED 0 when Solve[{eq}, {x, y}] (one
      equation, several unknowns) solved instead of declining -- the stage trusts a Solve shape; it should
      diff-back verify (seen while trying the {eqn} unwrap for several variables; the unwrap is one-variable only)
- [ ] DSolve corpus 2.2.33-3275 (y y'' - y^2 y' == y'^2, y[0]==2, y'[0]==1) and -3277 (y y'' == y^3 + y'^2,
      y[0]==1, y'[0]==2) regressed at b3c5eb34 (v0.267): the general solution became the Tanh family
      -C[2]/2 + C[2]/2 Tanh[C[2](C[1]-x)/2], which meets these conditions only at complex constants (the old
      E^(C[2] x)/(C[1] - E^(C[2] x)/C[2]) fits with real ones); 3277's implicit form leaves C[2] unfitted
- [ ] The corpus prelude accepted a WRONG implicit answer for 2.2.37-3683 before v0.355 (residual 1.14 by the
      implicit-function rule) -- its implicit verify can be fooled; and 2.2.26-2524..2526's old implicit
      Riccati answers were wrong too (now declined, correctly)
- [ ] Exact's potential integral pays ParallelMixedSpecial's 10 s cascade budget to decline (2.2.8-753 at 10.6 s
      against the corpus' 8 s); recheck after the v0.6038 tower fix
- [ ] Element[x, dom] decides from $Assumptions during plain evaluation (simp_builtins.c builtin_element), so
      Block[{$Assumptions = Element[a, Integers]}, Refine[e, Element[a, Integers]]] passes Refine a True and
      $Assumptions itself evaluates to True.  Only assumption-aware functions should consult it; check what relies
      on the current behaviour before changing it (Limit works around it, v0.6044)
- [ ] Gruntz on Exp[Exp[PolyGamma[0, PolyGamma[0, x]]]]/x declines; with Series' coefficient-expansion threshold
      at 64 leaves (series.c so_normalize_radical_coefs) it answers 1/Sqrt[E] (correct) but in 109 s of Expand; the
      threshold is 128 for now.  A cheaper normal form than Expand for those coefficients would get both
- [ ] Refine[Tan[x] > 0, 0 < x < Pi/2] stays unevaluated though Refine decides Sin[x] > 0 and Cos[x] > 0 there
      (so Simplify[Log[Sin[x]] - Log[Tan[x]] - Log[Cos[x]], 0 < x < Pi/2] does not fuse); prov_pos derives Tan's
      sign from Sin/Cos facts but not from a box on x
- [ ] Root is HoldAll (src/root.c root_init), so its INDEX is never evaluated: k = 2; Root[#^2 - 2 &, k] and
      Table[Root[f, k], {k, 4}] stay symbolic.  Function holds its own body, so HoldFirst (or none) suffices
- [ ] ClassifyPrime[T, c] does not terminate for a constant c (TowerSpecialPlaces now filters constants out;
      the function itself should decline)
- [ ] dsolve_tests / m61 / m63 run 150-200 s (alarm raised to 400 s in each main); the costs are no-crash /
      implicit checks that were as slow at v0.350, and F1's Bessel quadratures now closing (2-15 s each)

## W6 Keep it from drifting again
- [x] CI job: .github/workflows/tests.yml -- `ctest -LE nightly` (unit suites) on push/PR, the full suite daily

# Dialog[] (2026-10-26)

Ground truth from the local 14.0 kernel (interactive mode, piped and pty): numbering continues from
the global $Line (+1), Dialog[expr] shows expr as the first dialog Out line, only an INPUT of the
form Return[] / Return[x] exits (a Return produced by evaluation does not), Return[] gives the last
line's value (a `;` line's value counts: `8;` then Return[] is 8), dialog lines are top-level
(Throw::nocatch, Break::nofwd, Abort[] -> $Aborted stay inside), prompts indent one space per level,
on exit numbering reverts, Options[Dialog] = {DialogProlog -> Null, DialogSymbols -> {}},
Dialog::optx / ::nonopt / ::lvlist leave the call unevaluated, $IgnoreEOF = False.

- [x] modular.c: C-side Block API (mth_block_enter/mth_block_leave), builtin_block refactored onto it
- [x] eval.c: exported eval_report_uncaught() (the top-level sentinel chain, used twice in evaluate_core)
- [x] src/dialog.{c,h}: session helpers ($Line / In / Out) moved from repl.c; reader/writer hooks;
      the dialog loop; builtin_dialog (HoldAll, Protected), options + messages, $IgnoreEOF,
      DialogSymbols/DialogProlog (Protected), docstrings, sym_names entries
- [x] repl.c: interactive reader (readline, add_history) and a stdin reader for -file; pipe mode none
- [x] tests/test_dialog.c (scripted reader + capturing writer): numbering, %, Return forms, EOF,
      $IgnoreEOF, $Epilog, $MessageList, DialogSymbols/Prolog order, nesting, Throw/Break/Abort,
      $Pre/$Post, Quit, bad options, Protected, no reader; valgrind/leaks clean
- [x] docs: control-flow.md, changelog, version bump + tag

Review: 26/26 dialog_tests; leaks --atExit 0 for dialog_tests and a -file script with nested dialogs;
evaluator / Block / Return / Catch-Throw / Abort / Goto / TimeConstrained / message suites pass;
interactive pty transcript matches the reference session line for line; check-c99, check-messages OK.
- [ ] Follow-up (pre-existing): the top-level REPL prints `Out[n]= Null` for `x = 5;` / `Print[..]` and stores
      Out[n] = Null for a `;` line (the dialog loop does neither); repl.c process_input could use the dialog's
      line handling
- [x] repl_hooks_tests: test_exit_requests_quit_with_status expected the input after Exit[3] to run (value is Quit[3]);
      the $PreRead test's stderr -> /dev/null -> /dev/tty restore hid every later FAIL line off a terminal
- [x] 14 more test files restored stderr with freopen("/dev/tty"): now test_stderr_redirect/restore (dup2) in
      test_utils.h, and test_symbol.c mutes via mth_msg_suppress_push/pop
- [x] 35 test binaries pointed stderr at /dev/null for the rest of the run: now test_mute_messages() (a whole-run
      Quiet frame), stderr stays open
- [x] nminimize_tests killed by the alarm(120) watchdog (171 s): the penalty path interpreted every constraint and
      its gradient (393f2767's active-set BFGS multiplied the iterations), and Abs' sent gradients to finite
      differences. Constraint programs registered (fm_progs_register), Abs' -> Sign in the solver gradient; 54 s ; v0.6026
- [x] test_fixed_charge_flow 33 s in nminimize_tests vs ~12 s standalone: the CMake test build had no -O (no build
      type); now -O2 -g by default (asserts kept). Full-suite A/B: 0 new failures, 4 fixed
- [ ] Follow-up (pre-existing, identical at -O0 and -O2): 42 of 417 registered tests fail on main -- 20 DSolve
      corpus/stress, 16 "Subprocess aborted" (numeric, simp, stats, fit, packed_list, ndarray_linalg, ml_dist,
      normalize, risch_transcendental, diffunderint, plot3d, sf_machine_error, legendreq, expint_machine,
      integer_assume_simplify, trigexpand), 3 SIGALRM (dsolve_m63_stress, simp_sf_airy, mathar_regression),
      2 Timeout (mathar_definite_corpus, dsolve_corpus_2_1_2), flint_bridge

# FeynCalc H1: named characters, test harness, first baseline (2026-10-24)

Plan: ~/.claude/plans/let-s-continue-our-feyncalc-humble-karp.md. Gap doc: packages/feyncalc/FEYNCALC_GAPS.md (local).

- [x] P1 named characters in source: WL classification table, lexer (letters/operators/constants/space), tests, docs ; v0.587 (full-suite A/B vs HEAD clean; match_stress_corpus 2000-deep case overflows only the unoptimised main test build)
- [x] P2 harness: tools/feyncalc/run_fc_tests.m, tools/check_feyncalc.py, make check-feyncalc
- [x] T2 Check/Quiet nesting, Quiet/Check message specs, named emission of legacy helpers ; v0.588 (A/B Release vs Release: one regression found and fixed -- Inverse::sing via an anonymous helper)
- [x] T3 assignment LHS pattern-free parts evaluated (h[y_ - 4]) ; v0.589 (Contract 67 -> 92 / 119)
- [x] T7 Cycles/PermutationGroup/GroupElements/Permute over a group ; v0.590
- [x] P3 first baseline (11362 pass, ~85% of runnable) + three triage agents (NCA+Tables, Dirac+Feynman, Shared+LoopIntegrals) -> FEYNCALC_GAPS.md §7
- [x] Batch 1 ; v0.591: -(x+y) sign distribution + flat-chain parse provenance (+ Root sign, negative coefficient under radicals, FactorList sign, series/FunctionExpand cancelling pairs), HoldPattern LHS, ::write, N-values, indexed iterators, valued pattern heads, OptionsPattern guard/binding, 2ep lexer. Full A/B Release vs Release clean except pre-existing Fit WorkingPrecision sequence bug
- [ ] Follow-up (pre-existing): Fit[..., WorkingPrecision -> 24] returns the exact fit after a WorkingPrecision -> Infinity call in the same session
- [ ] Follow-up: positive non-reducing coefficient under radicals (Sqrt[12 a] -> 2 Sqrt[3] Sqrt[a] in the reference)
- [x] WL reference run of the harness (--per-dir --binary wolframscript): 13412/13412 runnable pass, 235 s -- the harness is faithful; per-dir is the comparison mode (FeynCalc caches across tests)
- [x] Batch 2 ; v0.592: Sequence splice only at pattern-variable sites, nested option lists (OptionValue/OptionsPattern; fixed the LeviCivita loop v0.591 exposed), Coefficient negative orders / denominator terms / expand only in the form; TimeConstrained unwinds only at a safe point (grace for the evaluator poll, then a frame-pointer check that no frame is inside libc) -- the Lorentz directory died with SIGILL after a timeout, losing 880 tests (now 1738/1812)
- [ ] O1 canonical order (src/sort.c expr_compare rewrite, worktree o1wt): numbers by value (complex: re, |im|, sign), Times > Power > Plus compared as polynomial terms from the last factor, atoms before compounds, strings before symbols, symbol collation, head/length/args. Random corpus 400/400, held-out 567/569 (2 misses are evaluation differences). Full-suite A/B and expectation triage against the reference pending
- [ ] Follow-up: DummyIndexFreeQ on a 50000-term sum: 25x the reference -- re-evaluation churn (Block bumps the eval clock, the big Plus is re-walked per call) and matcher speed
- [ ] Follow-up: Infinity / ComplexInfinity / -Infinity are symbols / Times[-1, Infinity] here, DirectedInfinity[...] in the reference
- [ ] Follow-up (pre-existing, wrong answer): PossibleZeroQ[D[F[x+1],x] - Cosh[x+1] Cosh[2x+2] + 10^-2] is True (reference False; F the Weierstrass antiderivative of Cosh[x] Cosh[2x]) -- every sample point lands at large |x|, where the operand scale (~10^40) hides the constant
- [x] v0.592 per-dir FeynCalc run: 12248/13412 (91.3%; Lorentz 877 -> 1738 after the TimeConstrained fix); baseline rewritten from it (per-dir mode)
- [ ] Follow-up: 32 ids pass per-file but fail per-dir (state left by earlier tests; reference passes them per-dir) -- mostly DiracTrick*FreeIndices-ID8, GluonVertex-ID1, FCFeynmanPrepare-ID51/53 (list: scratchpad fc_lost_v592.txt)
- [x] v0.593 PossibleZeroQ exact-point ladder; v0.594 Times numeric-base power merge; v0.595 assumption-rule pattern capture (Refine[1/Sqrt[r^2], r > 0] was -2)
- [x] v0.596 PossibleZeroQ climbs every rung on exact huge-scale points; v0.597 canonical order (O1) + Solve/SolveAlways ordering + Meijer-G convolution side + Simplify power-distribution polish
- [ ] O2 power normal forms (prototype in scratchpad o1fix: power.c numeric-factor extraction for any non-integer exponent, (1/d)^e -> d^-e; times.c rational-coefficient absorption 12 2^n -> 3 2^(2+n)). All match the reference, but they break recognisers that read the old forms: Simplify half-angle on Sqrt[(1 - Cos[x])/2] (now Sqrt[2 - 2 Cos[x]]/2), infinite Product families on (1/3)^(2^k), Gamma duplication Simplify, a Mellin Bessel integral, Goursat/Residue integrals. Needs those recognisers updated with it
- [ ] Follow-up: the symbolic-n Struve ODE identity takes ~22 s (1.2 s before the canonical order) -- the special-function engine families churn Factor/Together over an expansion
- [ ] Follow-up: numeric_sign (simp_assume.c) reads only Integer/BigInt/Real, so x > 1/2 does not imply x > 0 in fact_implies_strict_positive
- [ ] Follow-up: Sqrt[(2 I) x] -> Sqrt[2] Sqrt[I x] in the reference (complex coefficient)
- [x] v0.597 per-dir FeynCalc run: 12644/13412 (94.3%); 405 gained, 9 lost (DiracOrder-ID54 is dummy-index relabelling of an equal answer; FCFeynmanParametrize-ID28..30 and FCToTeXPreviewTermOrder-ID3 need the Times-level power split below; TID-ID42/44/65 and SetMandelstam-ID7 not yet traced)
- [x] v0.598 Multinomial + Binomial negative/half-integer closed forms (PaVeUVPart 137/139); v0.599 FormatValues (SharedObjectsUpValues 92/97); v0.600 Return from a definition (13 FCLoopBasis aborts); v0.6001 linear box syntax (48 ToExpression::sntx)
- [x] v0.6002 Infinity ordering, v0.6003 PolynomialGCD/Cancel over compound variables (x[i], Pair[..]), v0.6004 NValues + N[x, {p, a}] + N not Listable
- [x] v0.6004 per-dir FeynCalc run: 12951/13412 (96.6%); 317 gained, 10 lost vs v0.597: 6 FCLoopIntegralToGraph (per-dir state; per-file v0.6004 fails 8 where v0.597 fails 12), 4 timeouts under A/B load (DiracTrace-ID29, TdecGrc-ID3, TID-ID13/14 -- being timed)
- [x] v0.6005 O3 power-of-product split + Mellin fullsimp2 shifts only partnered Gammas
- [ ] v0.6006 exponent generators (Variables/Factor/FactorList over b^t, expsplit.c), Simplify exponent/argument normal form + cross-term coprime basis, O4 Times coprime basis for numeric powers -- combined full A/B running
- [ ] v0.6009 number atoms (rules/Cases/Count/Position/DeleteCases never enter Rational/Complex; FreeQ/MemberQ see their head) -- A/B running; v0.6010 DownValue order ignores conditions (FromGFAD 28), v0.6011 OptionsPattern in any position (SUNSimplify 10), v0.6012 SyntacticNegativeQ sums/complexes (SUNSimplify 2) -- written, A/B next
- [ ] Matcher: in a Flat+Orderless head, rest_ s[a_] binds the specific pattern to the EARLIEST factors in the reference (Replace[s[1] s[2] s[3], rest_ s[a_] :> {rest, a}] is {s[2] s[3], 1}; here {s[1] s[2], 3}) -- SUNSimplify ID83-86
- [x] DownValue order: a later repeated-name pattern p[a_, a_] goes ahead of p[a1_, b1_] defined earlier, as in the reference (v0.6022, generality order)
- [ ] Simplify[(18 - 6 n^2 + n^4)/(96 n^2)] keeps the input in the reference (here the expanded 1/96 (-6 + 18/n^2 + n^2)) -- SUNSimplify ID23
- [ ] Risch real-trig reconstruction (rt_realify): 1/(3 + 5 Cos[x]) (real zeros) comes back as complex-argument logs; the test's FreeQ[g, Complex] was vacuous until v0.6009
- [x] v0.6013 D of a compound head + Factor num/den via FLINT; v0.6014 ReleaseHold in heads; v0.6015 ToExpression non-string; v0.6016 Times O4 loop fix + merged-exponent sign (last term, complex) + numeric opposite exponents; v0.6017 Together/Cancel indexed variables (29 s -> 0.02 s); v0.6018 Factor Gaussian content (DiracTrace 53/61); v0.6019 contagion only numeric summands; v0.6020 pattern-headed call single operand in Flat -- all pushed 2026-10-26; FeynCalc per-dir measure pending
- [x] v0.6020 per-dir FeynCalc run: 13132/13412 (97.9%); 165 gained, 12 lost -- 9 FCFeynmanPrepare (v0.6009 FreeQ-on-Rational let Factor2 send 2 - D/2 inside Gamma to Together, which rewrote kernel arguments: fixed v0.6021, Together kernels opaque as in the reference), DiracTrickInsideDiracTrace-ID17/29 (v0.6010 exposed score-based DownValue order: rule order by pattern generality, v0.6022), FCFeynmanFindDivergences-ID24 timeout (not yet traced)
- [x] Perf: memoised definitions were quadratic (f[n_] := f[n] = n over 20000 values ~38 s); per-symbol hash index of the pattern-free prefix -> 0.024 s ; v0.6023 (targeted rule suites + reference diff of Unset/Clear/Block/DownValues=/conditional literals; full-suite A/B still to run)
- [ ] v0.6022 per-dir FeynCalc run in progress (scratchpad fc6022): Dirac 1881 (+7), Shared 4621 (+5), SUN 162 (+2), Pauli 225 (+1) so far; then --update-baseline (local)
- [ ] FCFeynmanFindDivergences-ID24: 15 s alone on v0.6008 and v0.6022 (passes); 120 s TIMEOUT only in the per-dir run -- per-dir state or load, not a kernel regression
- [ ] DownValue redefinition: the reference also replaces equivalent-but-not-alpha-equal patterns (x_ then _ ; x_h then h[___]; x:(_Integer|_Real) then x:(_Real|_Integer)) -- one rule there, two here
- [ ] Apart: splits over a linear factor with radical or I constants only when they are symbolic -- Apart[1/((x - Sqrt[2])(x - Sqrt[3])), x], Apart[1/((x - (-1)^(1/3))(x + (-1)^(1/3))), x], Apart[1/((-i Sqrt[3] + x)(i Sqrt[3] + x)), x] are recombined here, split in the reference; Apart[g[x]/((x - 1)(x + 1)), x] (non-polynomial numerator) declines here; Apart[(z^2 - 1)/z, y] (variable absent) splits here, the reference returns the factored Together -- ToFCPartialFractionForm 6
- [ ] Solve[1 + z^2 + z^4 == 0, z] gives nested Sqrt[1/2 (-1 - I Sqrt[3])] here, (-1)^(1/3)-style roots in the reference; Roots[...] unevaluated here
- [ ] (2026-10-26 finding) the reference's multivariate factor sign is NOT a monomial order: 2-variable cases are lex (x > y), but b - a c / a - b c / a b - c and b c - a d all keep the term without the alphabetically largest variable positive, and no lex/grlex/grevlex permutation or numeric-leading-coefficient heuristic beats plain lex on 559 random factors (493/559) -- an artifact of its factoring algorithm; deprioritised (3 DiracTrace tests). Unit content does differ simply: Factor[4 I a - 4 I b] is (4 I) (a - b), Factor[-I a + I b] is -I (a - b)
- [ ] Factor over Z[i] when the input has complex coefficients: Factor[I x^2 + I] is (1 + I x) (I + x) in the reference (I (1 + x^2) here); unit choice for general complex content ((-4 + 2 I) a + (6 - 2 I) b is (2 I) ((1 + 2 I) a - (1 + 3 I) b) there); Factor[(I x + I)^2] gives -(-1 - x)^2 (pre-existing); Factor of a fraction with Gaussian content (num factored internally, not via builtin)
- [ ] Simplify's RadicalCanon re-dispatch (simp_search.c ~1308, `<=` gate) recurses without bound whenever radical_canon and Times disagree on a canonical form (M11338 stack overflow, fixed at the Times side in v0.6016) -- consider a cycle guard
- [ ] Plus/Times: none open from contagion; Flat ReplaceList singleton wraps (fl[a]) still differ (documented)
- [ ] Simplify tie: (-3 + d)(-2 + d) and 6 - 5 d + d^2 both score 9; update_best keeps the first (strict <); the reference returns the factored one (LeafCount 7 vs 8?) -- unverified tiebreak, wide blast radius (DiracTrace ID71)
- [ ] Factor sign/unit convention: Factor[8 b c - 8 a d] is -8 (-b c + a d) here, 8 (b c - a d) in the reference (FactorList already agrees); I-content not pulled: Factor[4 I a - 4 I b] is (4 I) (a - b) there -- DiracTrace ID45/53/61
- [ ] Factor with Gaussian-integer coefficients is very slow multivariate (18 vars: 8 s vs 0 s) -- DiracTrace ID46/29 timeouts
- [ ] (f -> Times sub-product part fixed v0.6020; remaining: the reference enumerates subject SUBSETS in order and assigns pattern elements within each, Mathilda is pattern-element-major, so q[x_, z_] orientation is retried per f there, not here) Matcher: an UpValue on a Times pattern q[z_, x_] f_[y__] (q Orderless) binds f_ to a sub-product of the Flat Times and tries orientations in a different order (q[a, b] e[a] g[b]: e[a] r[b, g] here, g[b] r[a, e] there); condition evaluations 2/16/48/320 for n = 1/2/3/5 vs 1/4/9/25, plus a second full pass -- DiracTrace ID70 and ID29 (PairContract)
- [ ] Together of a sum with kernel factors does not combine: Together[f[a] (-p - q) + f[a] (p + q)] is 0 in the reference (unchanged here; with a symbol it works); the symbol path in turn expands a single product (Together[x (y + 1)] is x + x y here, x (1 + y) there) -- see feedback_together_no_expand: wide blast radius -- DiracTrace ID61#2
- [ ] Simplify should prefer the factored form: Simplify[(2 - d)(3 - d)] is (-3 + d)(-2 + d) in the reference (6 - 5 d + d^2 here) -- DiracTrace ID71
- [ ] AbsoluteTime[] has 1 s resolution; Unique[{"a","b","c"}] gives all three the same number
- [ ] Factor of a large rational expression is slow (FactorList2-ID10: Factor[q] 51 s while Together 1 s + numerator/denominator factoring 0.4 s): the classical exact_poly_div path in builtin_factor's rational branch; FactorList2 ID5/ID10 time out
- [ ] Coefficient[x^(1+e) + 3 x^(2+e), x] is 0 in the reference (x^e here)
- [ ] FactorList order with a radical generator: {Sqrt[x+y], 1} before {(x+y)^e, 1} in the reference
- [x] O3 Times-level power split (v0.6005), measured against the reference: in a product, Power[Times[c, rest], e] with c a positive number splits into c^e rest^e (y (2x)^b is 2^b x^b y; Pi^a (2 Pi)^b is 2^b Pi^(a+b); (-2x)^b x^a is 2^b (-x)^b x^a), while a lone (2x)^b stays. Needed by FCFeynmanParametrize (Pi^(D/2) (2 Pi)^-D) and FCFeynmanFindDivergences
- [x] Variables/Factor/FactorList split b^(n + k e) into b^n (b^e)^k (FactorList[(x+y)^(-3+3e)] is {{1,1},{(x+y)^e,3},{x+y,-3}} in the reference) -- FactorList2 ID7/8/11, FCFeynmanFindDivergences
- [x] Simplify[a^(3(-1+e)) - a^(-3+3e)] is 0 in the reference (exponent normal form); Simplify[la^2 (la^2)^(e-1), la > 0] leaves la^(2 + 2(-1+e)) (assumption power rule result not re-simplified) -- FactorList2 aborts, FCFeynmanProjectiveQ
- [ ] Together does not pull out integer content: Together[-2 a x6 + 2 b x8] is -2 (a x6 - b x8) in the reference
- [ ] Indeterminate < 1 is False here (IEEE-unordered, deliberate) but unevaluated in the reference
- [ ] Parser: Hold[_Symbol - 1] keeps Times[-1, 1] (the reference folds the literal to -1)
- [ ] P4 acceptance loop: fix largest root-cause groups in the kernel, ratchet

# SparseArray through the functional and procedural builtins (2026-10-24)

Plan: ~/.claude/plans/let-s-continue-our-implementation-clever-frog.md. Spec: docs/spec/builtins/data-structures.md ("Functional and procedural heads").

- [x] Reference ground truth (wolframscript) for ~120 functional/procedural call shapes
- [x] W1 level specs: native Map (leaf-level O(nnz), shape-preserving resparse), MapIndexed/Apply/Scan with a level spec (src/sparse/sparse_func.c)
- [x] W2 iterators: Table/Do/Sum/Product over {i, s} (iter.c, one fix point)
- [x] W3 ordering/set/fold heads answered on Normal (Sort/SortBy kept sparse; Fold keeps a sparse initial value); atom ::normal for Sort/SortBy/Union/Intersection/Complement/Split/DeleteDuplicates (builtin_normal_error)
- [x] W4 AppendTo/PrependTo on sparse; s[[All,k]] = sparse; shape-breaking assignment leaves a List
- [x] W5 Thread over sparse (incl. Thread[s]); Flatten DEEP; Outer sparse only for Times/one-arg; Times with dense arrays sparse; lower-rank sparse broadcast (m + v)
- [x] W6 Clip/Chop element-wise; resparse background = commonest element when the leading one vanished
- [x] listable_sparse: densify cap + Thread::tdlen settle (Range[0] bigv no longer allocates 10^9)
- [x] Dense fixes found: Map/Apply/Scan atomic Rational/Complex; Map All; packed w[[All]] = 0
- [x] tools/sparse_head_sweep.py (sparse vs dense vs reference, 958 heads x 10 shapes); fuzzer func/proc families
- [x] tests/test_sparse_func.c
- [ ] Follow-up (pre-existing, dense): Order/Sort compare Lists from the LAST element (Plus/Times canonical order) -- Sort[{{0,x},{1,0}}] gives {{1,0},{0,x}}; reference {{0,x},{1,0}}. Core expr_compare change, large blast radius.
- [ ] Follow-up (pre-existing, dense): ~700 dense-side differences in the head sweep (DENSE class), e.g. Floor[0, 0] -> 0 (ref Indeterminate), AnyTrue with a non-Boolean test stays unevaluated (ref builds Or), Thread with unequal lengths gives no Thread::tdlen. Re-run: python3 tools/sparse_head_sweep.py
- [ ] Follow-up: Map[f, 1/2, {0, Infinity}] returns 1/2 (ref f[1/2])
- [ ] Follow-up: First[x, s] / Last[x, s] with a sparse default densify the default (head only)
- [ ] Follow-up: Table[i, {i, v}] is a List; the reference re-sparsifies a one-dimensional Table (documented)

# SparseArray stress test, fixes, and beating scipy.sparse (2026-10-24)

Plan: ~/.claude/plans/let-s-thoroughly-stress-test-indexed-hopcroft.md.

- [x] Phase 0: head sweep (1247 heads x 22 call shapes) -> gate modes; `tools/sparse_stress.py` differential fuzzer
- [x] Phase 1: transparency gate (`src/sparse/sparse_gate.c`, `SymbolDef.sparse_mode`, eval pre/post hooks, 2^27 cap + SparseArray::dense, MATHILDA_SPARSE_DIAG=gate)
- [x] Fuzz-found fixes: packed 0^0, `{0.,1.}^0` -> 1, `2.5^0` -> 1., -0.0 printing, symmetric Eigenvalues tie order, linalg FullForm messages, double message on declining builtins (evaluator settled step)
- [x] Full suite A/B vs HEAD worktree (core evaluator change) before push -- clean after fixing the settled-step rule (pointer snapshot) and the eigen test; pushed v0.574-v0.581
- [x] Phase 2a: native reductions (v0.577), structural ops (v0.578), Dot/SpGEMM/MatrixPower/Inner (v0.579)
- [ ] Phase 2b: sparse LinearSolve direct+Krylov, Eigenvalues[s, k] Lanczos/Arnoldi (worker running)
- [x] Phase 3: benchmark 120 broadened, 3-way, every row <= 1.0x
- [x] Phase 4: docs, valgrind on the fuzz corpus, `make check-sparse`
- [x] `x = rhs` re-evaluated its value after the assignment (double messages, Quiet leak) -- v0.581
- [x] Fuzz/worker-found fixes: constructor NULL-side crash (v0.580), Unitize (v0.582), dense KroneckerProduct/triangularize/Transpose/Reverse/IdentityMatrix/Flatten (v0.583), MatrixRank generic rank (v0.575)
- [ ] Follow-up (pre-existing): `Eigensystem` does not exist
- [ ] Follow-up (pre-existing): non-symmetric machine Eigenvalues tie order: {{1.,2.},{0.,-1.}} -> {1.,-1.}, expected {-1.,1.}
- [ ] Follow-up (pre-existing): `ReplaceAll[s, 1]` (a non-rule) returns s; it should stay unevaluated
- [x] `MatrixRank` of a 5x5 mixing machine reals, a symbol and 25-digit integers ran > 30 s -- v0.575
- [ ] Follow-up: `a = Range[6]; Dot[a, {1,2}]` prints Dot::dotsh twice (builtin declines a packed arg -> post-gate unpacks -> re-run)
- [ ] Follow-up: `Take[{1,...,6}, 10]` emits no Take::take
- [ ] Follow-up: dense exact matrix 2-norm `Norm[{{1,0,0},{0,0,2},{0,-1,4}}]` unevaluated
- [ ] Follow-up: `nd_parallel_for` spawns threads per call (no pool); limits small-op speed (Max, A.M10)

# SparseArray (2026-10-23)

Plan: ~/.claude/plans/pasted-content-id-bdcb-let-s-plan-smooth-shore.md. Spec: docs/spec/builtins/data-structures.md.

- [x] Semantics probes (canonical form, raw-object matching, construction, merge rule, Band, Part, properties, messages)
- [x] EXPR_SPARSE node type (src/sparse/): payload, shared pattern, canonical kinds; expr/eval/match/print/sort/precision plumbing
- [x] Constructor (rules, {pos}->{vals}, patterns, Band incl. blocks, dense/packed, existing arrays, standard form) + messages
- [x] COO -> CSR: composite-key LSD radix (12-bit digits), first-wins dedup, background drop; packed-rule direct keys
- [x] Normal (packed when one machine type), Length/Dimensions/ArrayDepth/ArrayQ/VectorQ/MatrixQ, Part read + assign, properties, ArrayRules
- [x] Listable threading (union rule, Times intersection), typed Plus/Times/Power incl. fused merges, ND-kernel reuse, generic tier, N
- [x] Tests: tests/test_sparsearray.c (constructors, messages, raw object, structure, Part, arithmetic, functions, fast-path cross-checks); valgrind clean (startup baseline only)
- [x] Audits: check-c99, check-messages, check-packed-aware, check-compile-coverage, check-nd-surfaces, check-array-exactness
- [x] Docs: data-structures.md, changelog v0.572, docstrings (SparseArray, ArrayRules, ArrayQ, ArrayDepth, Band, Normal)
- [x] benchmarks/120-sparse-arrays vs scipy.sparse: 13 cases, all at parity or ahead (v0.573)
- [ ] Follow-up session: linear algebra (Dot/SpMV/SpGEMM, Transpose, Total, LinearSolve, Det, Inverse, eigen), Map/Flatten, graph heads returning SparseArray
- [ ] Follow-up (pre-existing, not sparse): Sqrt[-4.] gives 1.22e-16 + 2. I (scalar and packed); the exact answer is 0. + 2. I
- [ ] Follow-up (pre-existing): standalone `r = pos -> vals` unpacks both sides (Rule is not packed-aware by the no-nesting invariant)

## Review
- Element-wise results equal the operation on Normal[...] element by element (no whole-operand Real coercion) -- documented in data-structures.md.
- Performance at or ahead of scipy.sparse on every measured op (build 40 vs 41 ms; A+B 5.8 vs 10; A*B 3.0 vs 7.2; Sin 2.8 vs 16).
- Pre-existing failures seen in neighbouring suites (packed_list ToNDArray mixed, numeric N[E^(I Pi/4)-(-1)^(1/4),32], series Integrate[Series]) reproduce identically on a HEAD build.

# Special-function Simplify engine (2026-10-22)

Plan: ~/.claude/plans/similar-to-our-simp-bessel-merry-pudding.md. Derivations: docs/design/simp_special_functions.md.

- [x] W0 engine (`simp_sf_engine.{c,h}`) + Bessel port; Casoratians, K evenness, integer reflection
- [x] W1 Gamma / incomplete Gamma / LogGamma / BarnesG / Pochhammer / Beta / Binomial / Factorial; reflection, multiplication
- [x] W2 PolyGamma / HarmonicNumber / HurwitzZeta / Zeta[s,a] / LerchPhi / Riemann FE / PolyLog
- [x] W3 Erf family, InverseErf, ExpIntegralE, LogIntegral/Ei, Si/Ci/Shi/Chi, Fresnel
- [x] W4 orthogonal polynomials, Airy, StruveH, Whittaker, Fibonacci/LucasL
- [x] W5 hypergeometric (PFQ canonical order, contiguous lattices, Kummer/Euler)
- [x] W6 elliptic K/E (Legendre relation), ProductLog, QPochhammer
- [ ] Follow-up: N[LaguerreL[nu, -1, x]] stays unevaluated (numeric generalized Laguerre at a = -1)
- [ ] Follow-up: Simplify bottom-up spends minutes in TrigReduce on `E^-w (poly in w^(2a+j))` (no special function involved); see test_simp_sf_gamma.c GammaIncomplete note

# Group F crypto benchmarks (2026-10-21)

Plan: ~/.claude/plans/let-s-create-20-additional-quiet-clover.md. Roadmap: benchmarks/CRYPTO_GAPS.md.

- [x] Harness: PARI/GP system in run_all.py, 3-digit experiments, `--only` ranges, harness.gp, crypto_data.{m,py,gp}
- [x] 20 experiments 100-119 (.m/.py/.gp + README each), labels and checks agree across 4 systems except real Mathilda bugs
- [x] Full four-system run; CRYPTO_GAPS.md ranked roadmap; README / requirements / changelog
- [ ] NEXT (from CRYPTO_GAPS.md, in order):
  - [ ] C1 Return inside While/For exits only the loop
  - [ ] C2 AppendTo/PrependTo no-op on packed lists
  - [ ] C3 MapThread unevaluated on substituted packed arrays
  - [ ] C4 Partition[l,n,d,kL,{}] unevaluated
  - [ ] C5 CoefficientList[p,x,n]
  - [ ] C6 Quiet leaks Inverse::method / Thread::tdlen (stdout)
  - [ ] C7 Module temporaries never removed (leak + progressive slowdown)
  - [ ] C8/P3 MultiplicativeOrder discrete log: Pohlig-Hellman + BSGS/rho, no cap
  - [ ] A1 Bit* family (ND + Compile surfaces) — [x] BitXor v0.508 (all surfaces, bignums, ::int); [x] BitAnd v0.510 (all surfaces, bignums, WL-verified simplification, no message); [x] BitOr v0.511; [x] BitNot v0.512 (all surfaces incl. autocompile, bignums, involution, ::argx; NumPy-level); BitShift* open; A2 Hash SHA-2/MD5; A3 Modulus-> linear algebra
  - [ ] P1 integer ListConvolve; P4 packing through Quotient/Map/Dot; P2 Association insertion
  - [ ] P5 fp-LLL; P6 machine-int Sqrt; P8 Times product tree; P9 InterpolatingPolynomial Modulus; P10 PrimitiveRoot; P11 PolynomialRemainder x^n+1

## Review
Results: see benchmarks/CRYPTO_GAPS.md headline. All headline bugs re-verified by the lead on the current binary; the 110 Fermat root cause was corrected (machine-int `Sqrt`, not bignum sqrt).

# Kernel fixes from the book pass (2026-10-20)

Plan: ~/.claude/plans/let-s-update-the-mathilda-cryptic-bee.md (W1-W7). Parallel worktree workers; lead integrates serially with version bumps v0.482+.

- [x] W1 Integrate: spurious Infinity::indet (Check), Exp[-xt]/t^2 hang, zeta(3), PV symbolic, PV Log[b/a], E1 family
- [x] W2 MeijerG/PolyLog reductions (Erf, BesselK, Pi/2, Pi/4, Catalan, -7 zeta(3)/16)
- [x] W3 DSolve vpS leak + answer simplification; Limit DirectedInfinity
- [x] W4 PolyLog[0,2,2.] real; D StruveH; D PolyLog[n,p,z]; StruveH half-integer; HypergeometricU[-n,..]; DirichletL exact Interval
- [x] W5 Minimize: TimeConstraint overrun + separable const factor; simplified optima; MZ_MAX_INEQ cap; open-boundary pole; slow LPs; docstring
- [x] W6 StringCases $n; ToCharacterCode enc; IntegerString Roman; TakeLargest ties; Definition ^=; Dot printing; ::tagpos; Plus/Times honour cleared Orderless
- [x] W7 Compile: Boole; exact Chop tolerance
- [x] integrate serially (v0.482-v0.507, tagged), full-suite A/B vs fe2fcb54: 35 failing suites identical to baseline (same FAIL lines)
- [x] book follow-through (a4bb2d19; 564 pp, check-links OK)

# Book update v0.290 -> v0.480 (2026-10-20)

Plan: ~/.claude/plans/let-s-update-the-mathilda-cryptic-bee.md. Backlog + status: book/BOOK_UPDATE_REVIEW.md.

- [x] Step 0 backlog (BOOK_UPDATE_REVIEW.md)
- [x] Wave 1 REVIEW fixes R1-R5, R7, R8, R10 (R6 -> A, R9 -> B)
- [ ] Wave 2/3 writers A-G (parallel, disjoint chapters)
- [ ] spurious Integrate::nonelem kernel fix (worktree agent)
- [ ] site index regen, check-links, usage, pdf, index terms, changelog, commit

# RandomPrime (v0.468, 2026-10-19)

Plan: ~/.claude/plans/pasted-content-id-4ce0-let-s-implement-steady-widget.md

- [x] shared PrimeQ test mth_is_prime_u64/_mpz (facint.c); PrimeQ routed through it
- [x] RandomPrime in random.c (table / u64 rejection / mpz rejection + existence scan), packed output
- [x] SYM_RandomPrime, Protected, docstring, SKIP_DIMS entry
- [x] tests (test_random.c, 10 cases incl. primality cross-check)
- [x] docs (random-number-generation.md, changelog), version bump
- [x] valgrind (per-call clean), check-c99/messages/packed-aware, PrimeQ test binaries (primenu stale expectation is pre-existing)

# Mathar integral corpus + Integrate extensions (campaign, 2026-10-18)

Plan: ~/.claude/plans/i-would-like-to-cozy-swan.md.  Source: arXiv:1207.5845 TeX
(scratchpad, sha256 in corpus header).  Tooling: tools/mathar/ (extract.py,
tex2m.py, verify.wls [Mathematica oracle], instance.wls, build_corpus.py,
check_override.py, overrides.json).

- [x] TeX->InputForm translator + extractor (1007 integral sites; 593 auto-translated)
- [ ] oracle verification of auto translations
- [ ] hand repair (sharded agents) of SKIP-PARSE / oracle-rejected entries -> overrides.json
- [ ] corpus files + runner test_mathar_corpus.c (S/N modes, parallel pool, JSON) + CMake
- [ ] baseline run; fix every WRONG / CRASH first
- [ ] triage doc docs/design/mathar_integrate_triage.md (classes ranked)
- [ ] extension stages, one per class (commit + version + tag each)
- [ ] backlog of remaining classes

# BlakeRationalBaseDescent rebuild (v0.456)

Plan: ~/.claude/plans/in-factorisation-temp-we-have-velvety-phoenix.md.
Sources: factorisation_temp/rational_base_descent_revised.tex, lattice_rbd.tex.
(v0.455 was taken by a concurrent strings commit.)

- [x] fi_rbd.h core API; fi_rbd_descent.c (Classic, Exact scaled + jump + batched gcd, AllHits, multiplier)
- [x] fi_rbd_bases.c (Stern-Brocot BFS, Sum order, perfect-power filter, simplest rational)
- [x] fi_rbd.c glue + hit cache; FiRbdOpts in FiCtx; suboptions + alias heads in fi_builtin.c
- [x] threaded deterministic base search (FiClaim cap, lowest index wins)
- [x] fi_rbd_lll.c (native integral LLL; FLINT fmpz_lll dispatch)
- [x] fi_rbd_lattice.c (HG detector, square-free part + p-adic integer roots, (t,d) selection); Lattice/TwoPass/Coppersmith
- [x] tests (paper examples, comb, multiplier, lattice 197-bit, edge cases, docstrings)
- [x] verify_methods.py --rbd-planted (975 calls, 0 wrong); perf old vs new; leaks 0; no-FLINT compile
- [x] extensive parent docstring (user request), docs, changelog, version 0.456

## Review (v0.456)
- All the papers' claims reproduced (independent Python reference): Ex.1 (512, -2/0), 2^1041+1 four hits,
  R_2224 14 hits, (153^153+1)/154 7 hits, comb sweep 533 bases -> 22/7 and 11/7 only.
- Exact kernel 37-100x faster per base than v0.454 (single thread); threads on top.
- The paper's 2048-bit challenge modulus factored with no base given: 106 s wall, 16 threads.
- Bugs found on the way: HG rows not square-free (p-adic finder needs h/gcd(h,h'));
  lattice regions too coarse for TimeConstrained; cofactor fallback cost on AllHits leftovers;
  verify generator looped forever when (a/b)^k/q was too small to hold a prime.
- Open: native LLL slow (0.2-0.7 s per call at 197 bits vs ms with FLINT).

# DirichletCharacter + DirichletL (v0.450)

Plan: ~/.claude/plans/pasted-content-id-6d9a-let-s-implement-stateful-locket.md.
Ground truth (wolframscript 14): DirichletCharacter attrs {Protected} (threads
over lists internally); DirichletL {Listable, NHoldAll, Protected} -- user asked
Listable|NumericFunction|Protected. 2^a (a>=3) digit order: (-1)-digit most
significant, then the 5-digit; prime-power digits smallest prime first. Value =
unit * E^(I Pi r), r in (-1,1] (Mathilda does not canonicalise E^(I Pi r), so the
engine folds). Messages: intp / invl / int. Exact: s<=0 always (B_{n,chi}),
s=1 primitive non-principal; s>=2 only k=1. Mathematica MR ~300-500 us, AP 8 ms.

- [x] character engine `dirichlet.{c,h}` (factor, PH discrete log, conductor, cache)
- [x] DirichletCharacter builtin + messages + list threading + packed path
- [x] DirichletL exact (k=1 Zeta, s<=0 gen. Bernoulli, s=1 primitive)
- [x] machine kernel `sf_machine_dirichlet.c` (E-M with expm1 tail, functional eq.)
- [x] FLINT acb_dirichlet_l bridge; HurwitzZeta template fallback; precision tracking
- [x] N keeps k, j exact (numericalize special case)
- [x] NDArray/packed path, AWARE, REG_N, Compile infer/emit (MeijerG-style)
- [x] Interval (arb ball + derivative-sign certification), real characters
- [x] info.c, sym_names, docs, changelog, version 0.450, audit tool probes
- [x] tests (test_dirichletcharacter.c, test_dirichletl.c) + valgrind
- [x] benchmark 97-special-functions

## Review (v0.450)
- Mathematica parity: 148,138 character values; 1152 exact L-values; all spec examples.
- Machine kernel <= 6e-14 rel (3000 random cases vs arb); five surfaces bit-identical.
- Throughput run caught a silent slow path: complex chi at real s < 1/2 declined in
  sfg_loggamma_c (negative real axis) -> arb fallback, 368 us/call. Fixed (real lgamma
  + sign phase); now 3.5 us.
- valgrind: totals identical for trivial / 1x / 4x workloads (startup baseline only).
- Found, not fixed (pre-existing, core): (1) SameQ of a packed complex array vs its
  evaluated List fails when an element has im == 0 (Complex[x,0.] vs x); (2) a
  message from an unevaluated call inside Quiet re-fires when the enclosing
  expression re-evaluates it (LerchPhi/HurwitzZeta/JacobiSymbol too); (3)
  test_numeric N[E^(I Pi/4)-(-1)^(1/4),32] and a stale AiryBi stress gap fail.

# MeijerG machine kernel for the log case (v0.446)

Goal: `MeijerGlog` (b's an integer apart) ran 119 s per 10^6 buffer / 456 us scalar,
because every point fell to the arb perturbation engine. Replace with an exact
machine residue-series kernel.

Method: G = -sum of residues of F(s) w^s at the poles s = beta + k of each
integer-spaced class of b_1..b_m. Laurent expansion of F at s = beta + eps
(arb, once per plan), then the exact rational step F(s+1) = F(s) R(s) in
double-double to tabulate c[k][j] (the coefficient of w^(beta+k) log^j w). Per
point: ascending sum in double, measured cancellation; re-sum in double-double
when > 6 bits are lost (log w, w carried as dd); arb only past 45 bits / table end.

- [x] `sf_dd.h`: dd + complex dd arithmetic, exp/log/sincos in dd
- [x] `flint_mg_laurent`: Laurent coefficients of prod Gamma(x_i + s_i eps)^(+-1)
- [x] log table build in `sf_mg_plan`; `sf_mg_plan_free`; memo owns tables
- [x] per-point evaluator (real + complex lanes, double then dd)
- [x] scalar path memoises log-case plans; ndarray/flat free correctly
- [x] tests vs mpmath (real, complex, dual p>q, multiplicity 3, r != 1, shifts), valgrind
- [x] surfaces bit-for-bit (scalar / List / packed / NDArray / Compile)
- [x] benchmark re-run, report + tex, docs/spec + changelog, v0.446

## Review
Pushed `bf002df8` (tag v0.446) and `f3c6185b` (benchmark). MeijerGlog 456 -> 2.99 us
scalar, 119 s -> 0.32 s per 10^6 buffer; 36 references <= 4.5e-16. Found on the way:
the dual's `1 - b` rounding was amplified ~100x (fixed in the log table; the non-log
Slater dual still forms it in double).

# Special-Functions Shootout benchmark

Plan: `~/.claude/plans/let-s-write-a-benchmark-optimized-hammock.md`

Self-contained 4-engine benchmark (Mathilda vs Python scipy+mpmath, Maxima, wolframscript)
over machine-real / complex / arbitrary-precision, measuring speed (scalar latency +
machine-real vectorized throughput) and accuracy (cross-engine agreement, no oracle).
Standalone folder `benchmarks/97-special-functions/` (numbered to shelve with the rest, but
skipped by `run_all.py` discovery via its `STANDALONE` set). No version bump.

## Build
- [x] Confirm engines present: Mathilda, venv scipy/mpmath/numpy, Maxima 5.49, wolframscript 14.0
- [x] Probe engine I/O formats (Mathilda InputForm/SetPrecision; Maxima interval/bigfloat/complex)
- [x] Verify scipy/mpmath callable names + sympy availability
- [x] `registry.py` — 46 numeric heads + 7 exact heads, per-engine templates + domains + callables
- [x] `run_shootout.py` — generate `.m`/`.mac`, run 4 engines, parse, agreement+timing → results JSON
- [x] Smoke subset across 4 engines × 3 regimes; fixed 3 parser bugs (WL backtick, dps-before-parse, Maxima complex)
- [x] Full numeric sweep (53 cases); disagreements all real (Maxima branch/symbolic), not parser bugs
- [x] Vectorized machine-real throughput (Mathilda packed vs scipy)
- [x] Exact-arithmetic section (Fibonacci/LucasL/BernoulliB/EulerE/HarmonicNumber/BarnesG/Hyperfactorial)
- [x] `report.py` — results JSON → LaTeX longtable/booktabs fragments + macros
- [x] `special_function_performance.tex` — standalone article; rendered PDF (10 pp)
- [x] `README.md`
- [x] changelog note `docs/spec/changelog/2026-10-12.md`
- [x] Send PDF to user + headline summary

## Review

Delivered `benchmarks/97-special-functions/` (numbered for shelving; `run_all.py` skips it via
its `STANDALONE` set): `registry.py`
(spine), `run_shootout.py` (4-engine driver), `report.py` (JSON→LaTeX), `special_function_performance.tex`
+ rendered PDF, `README.md`, `results/2026-10-17.json`. Contributor-facing tooling, no kernel change,
**no version bump**.

Results (53 numeric heads × 3 regimes × 4 engines):
- **Coverage:** Mathilda 152 head×regime cells = Mathematica 152 = Python 152 (mpmath fallback); Maxima 124
  (lacks HurwitzZeta/PolyLog/HypergeometricU/Whittaker numerics + orthopoly-at-complex).
- **Accuracy (agreement, no oracle):** MR/C ~13–16 digits, AP ~50; exactly **1** genuine disagreement
  (Maxima `expintegral_chi` adds a spurious iπ on the real axis). All 7 exact heads byte-identical.
- **Speed:** scalar latency — Mathilda ~12× (MR) / ~9× (C) slower than SciPy's compiled C per call, but
  **~2× faster in arbitrary precision** (often the fastest engine). Vectorized throughput (10⁶ buffer):
  Mathilda's packed kernels beat SciPy on most heads, up to **209× (StruveH), 136× (BesselY), 46× (BesselJ)**.

Key gotchas fixed: wolframscript `InputForm` backtick precision suffix; parsing CAS AP strings before raising
mpmath dps (truncated 50-digit values to 15); Maxima complex `a*%i+b` / nounform / `%pi` partial-eval; scipy
0-d-ndarray + complex-declining funcs → mpmath fallback; sympy-import charged to first exact timing.

---

# DSolve — the next 500 (§2.2.46–§2.2.50, Problems 4501–5000, M74–M78)

Plan: `~/.claude/plans/let-s-continue-our-implementation-twinkling-reef.md`. Executed as ONE
batch (fixes are shared across sections) in two commits: core CAS fixes (v0.434) + DSolve
batch with the five corpora and gates (v0.435).

- [x] Fetch + convert all five (§2.2.46 has 55 systems); indvar green; spot-audit →
      converter fix `\operatorname{arccot}` → ArcCot (4924, §2.1.2-756)
- [x] Baseline 89/83/76/85/90 = 423/500, 0 FAIL
- [x] General fixes (see DSOLVE_PLAN M74–M78): Together/Cancel binder + Log semantics,
      Factor SIGFPE, ∫Secⁿ, quiet Risch probes, inert integrating factor, Riccati
      particular, ExactODE inert, Kummer U, separable split + branch check, UC screen,
      system trig 2nd pass
- [x] Full per-record A/B vs HEAD (51 sections) → found + fixed 8 regressions incl. 2 FAIL,
      1 wrong-as-UNEVAL, 2 crashes; final: 0 lost, +41 old; new 94/89/85/88/92 = 448/500
- [x] Unit suite: only pre-existing HEAD failures; ratcanon test updated; c99/messages green
- [x] Gates 6/11/16/12/8; reports; STATUS/README/PLAN/changelog/calculus.md/algebra.md
- [x] valgrind on the new paths
- [x] commit + tag v0.434, v0.435

## Review — next 500 (§2.2.46–50)
423 → 448 / 500, 0 FAIL, 0 crash; older sections +41 with zero per-record losses. The
biggest wins were not DSolve methods but core correctness: Together/Cancel were silently
value-changing (binder substitution, Log expansion), Factor could crash the kernel, and
∫Secⁿ hung. Ten older gates are red on HEAD already (drift since M73) — flagged, not
rebaselined.

---

# MeijerG (+ minimal FunctionExpand) — v0.445

Plan: `~/.claude/plans/pasted-content-id-088c-let-s-implement-warm-map.md`.

- [x] Mathematica calibration probes (reduction policy, hdiv, z = 0, Gamma-pair move)
- [x] Machine Slater engine `src/special_functions/sf_machine_meijerg.{c,h}` (plan + eval, arb fallback)
- [x] arb engine in `src/flint_num_bridge.c` (regularized pFq, log case by adaptive perturbation)
- [x] Builtin `src/special_functions/meijer_g.{c,h}`: validation, sort, hdiv, threading, z=0, reductions
- [x] NDArray/packed (AWARE) + ndkernels nary kernel + Compile (cap 8 -> 24, real-arm imag check, VKERNN tile op)
- [x] D (deriv.c), Series (series.c; + Puiseux O-term fix), FunctionExpand
- [x] tests/test_meijerg.c + CMake; valgrind (REPL 1x/4x + test binary == Print[1] baseline; fixed times.c uninit)
- [x] audits (c99, messages, packed-aware, compile-coverage [HEAD_SHAPES], nd-surfaces, fastpath gate)
- [x] benchmark registry + run + report/tex (results/2026-10-17-final-v0445.json)
- [x] docs (special-functions.md, changelog), version 0.445, commit + tag

## Review — MeijerG (v0.445)
MeijerG with exact/symbolic, machine, arb and complex numerics; buffers; Compile at scalar
and rank-1 (new VKERNN tile op); D, Series (log case by exact eps-Laurent), FunctionExpand.
Fastest engine in all 9 shootout cells (mpmath 2-234x, Mathematica 11-670x slower).
Four core bugs surfaced and fixed: Series Puiseux O-term (x^(p/q), q>2), a Times<->Power
radical-fusion cycle ($RecursionLimit), a Times uninitialised read, NumericQ on list args;
plus MPFR 2F1 via arb (60 ms -> <1 ms near |z|=1) and a Compile real-arm wrong answer
(imaginary part dropped). Full ctest: 20 failures, all reproduced on HEAD (per-record A/B:
0 regressions, 11 corpus records U->P).

---

# DSolve — the next 500 (§2.2.51–§2.2.55, Problems 5001–5500, M79–M83)

Plan: `~/.claude/plans/let-s-continue-our-dsolve-snoopy-biscuit.md`. One batch, as M74–M78.

- [x] Fetch + convert all five (all 100 scalar, 0 IVP, 0 systems; Kamke-ordered first-order:
      Abel 2nd type / homogeneous / dAlembert / Clairaut / quadrature-in-y'); indvar green
- [x] Baseline (HEAD v0.445 binary snapshot) → reports
- [x] Triage + general fixes
- [x] Full per-record A/B vs HEAD on all 56 sections; 0 FAIL / 0 crash / 0 lost
- [ ] Unit suites, c99/messages, valgrind
- [x] Gates, reports, STATUS/README/PLAN/changelog/builtins docs
- [ ] Version bump(s), commit, tag, push

## Review — next 500 (§2.2.51–55)

- **Result:** 431 → **462 / 500**, 0 FAIL. Gates 18 / 4 / 5 / 9 / 7 (gap 17 / 4 / 4 / 8 / 5 + a named load margin). Three of HEAD's 431 were
  wrong answers (5005/5007/5008).
- **A/B** vs HEAD v0.449, all 51 older sections, same load: **+22** (passed by neither HEAD
  run). The 4 losses are HEAD wrong answers (2.1.2-777, 2.2.48-4777/4779/4792; residuals
  checked), plus 2 flutterers that pass on rerun. §2.2.48 gate 16 → 20 (+2 for the slow new gains 4774/4727), §2.2.35 2 → 1,
  §2.2.36 7 → 5.
- **The late-slot lesson:** the first full A/B showed 14 real losses. Bernoulli-on-`Expand`
  and Exact-by-μ(u) sat *inside* their early methods and claimed equations the
  specialists closed. Moving both to late slots, plus the linear integrating-factor
  constant drop, cleared all 14.
- **Open:** intermittent SIGILL under load. A TimeConstrained async jump lands inside a
  non-GMP malloc (`symmemo_lookup`, Expr nodes). It is pre-existing and occurs at the
  same rate on HEAD. The fix belongs to a TC design pass (all-cooperative jumps, or
  guarded allocation of Expr nodes).
- **Residue 38:** Abel 2nd type, Riccati/Chini with symbolic exponents, elliptic
  separables, symmetry classes, `y = G(x,y')`, lincoeff `F(L1/L2)`, the `u = E^(2y)`
  substitution.

## Integrate: Bessel hypergeometric / MeijerG stage (v0.544, 2026-10-23)

- [x] Derive and numerically verify D1-D7: series (Chu-Vandermonde 2F3), Sin reflection, term-by-term with Log split, Mellin-Barnes reps (Y, K, KK, IK, JY, J^2+Y^2), recurrence/Wronskian reduction, telescoping for Int D_n/x
- [x] `src/internal/BesselHyper.m` + `try_besselhyper` after Holonomic, `Method -> "BesselHypergeometric"`
- [x] `tests/test_integrate_bessel_hyper.c` (diff-back at 3 points, 1e-20)
- [x] Holonomic / integrate suites; check-c99, check-messages; valgrind (stage pinned: no new leak)
- [x] Docs (calculus.md), changelog, version 0.544

### Review
All 9 originally failing probes now close, plus integer/symbolic Y/K orders.
Boundary (declines): >2 factors, different arguments, mixed ordinary/modified
families, Mellin log-degenerate cases with nothing to telescope.
Pre-existing, not from this change: `integrate_risch_transcendental_tests`
fails on `x*2.71828^(-x)` (closed by ParallelMixedTower); a 64-byte
Module-scoping leak on `Integrate[BesselJ[0,x]^3, x]` in a later cascade stage
(both Bessel stages pinned are clean).

## FeynCalc remaining correctness gaps (2026-11-02)
Plan: packages/feyncalc/FEYNCALC_REMAINING_GAPS.md (195 tests, clusters G1-G13).
- [x] G1 Solve empty system / momentum shifts ; v0.6127
- [x] G6 PaVeReduce ; v0.6128
- [x] G5 GA::noint ; v0.6129
- [x] G4 PowerExpand/FunctionExpand, G9 Series var, G11 GroebnerBasis option ; v0.6130-31
- [x] G2 Factor sign convention, G3 InputForm printer ; v0.6132-33
- [x] G10 FeynRule/FunctionalD ; v0.6134; G8 mostly closed (v0.6136-39)
- [x] v0.6145 round 12: FindDivergences ID2/3 (Simplify numeric-basis pre-pass kept (la^2)^ep abstracted), SetOptions on symbol lists (PaVe BReduce), Together compound-arg atoms (NPointTo4Point), pure-function body sharing (Isolate 0.11x -> 0.52x); 13368, superset
- [ ] G7 Feynman parametrisation -- bisected 2026-11-03 (scratchpad g7a/g7b), causes:
  - [x] F Power: perfect-power integer base under a symbolic exponent (4^(D/2) -> 2^D, 12^(D/2) -> 2^D 3^(D/2), (1/4)^e -> 4^-e) -- Parametrize ID10/11 (power.c) -- v0.6146; ID11 also needed Simplify's exponent-content fold (4^(-2+2e) -> 16^(-1+e)) and D; ID11's prefactor now matches, its remaining diff is B inside the power base
  - [x] C Together splits the integer part of a symbolic exponent on a non-sign-normal base ((b - a)^(-1 + n) -> -(-a + b)^n/(a - b)) -- Parametrize ID26 (rat.c) -- v0.6146, negative integer parts only; (b - a)^(1 + n) -> -((a - b) (-a + b)^n) still open (the product re-merges into (-a + b)^(1 + n))
  - [ ] A Simplify[e, sum == 1]: per-subterm/coefficient reduction, compared before the stronger plain Simplify wins -- Prepare ID11/13/14 (simp_reduce_by_equations)
  - [ ] B nested Collect with simplified coefficients (x1 x2 + ... -> x4 (x3 + x5) + x1 (...) + x2 (...); m + s x - s x^2 -> m - s (-1 + x) x) -- Parametrize ID6/20 (simp_search.c try_collect_per_variable)
  - [ ] H factor a common base with integer-differing exponents (x^(a - 1) + x^a -> x^(-1 + a) (1 + x)), also Together -- Parametrize ID44/10
  - [ ] E fold a monomial into a power of the same product -- Parametrize ID9 (simp_power.c)
  - [x] D Legendre duplication on shifted Gamma ratios over-applied (Simplify[Gamma[-2+2e]/Gamma[-2+e]] stays in the reference) -- Parametrize ID11 -- v0.6146 (simp_sf_engine.c sf_improves: a score rise needs a whole family eliminated)
  - [x] DECIDED (2026-11-03, user): `3e-3` is the Real 0.003 (Mathilda's e-notation extension stays) and `3 e-3` is Plus[Times[3, e], -3], which the parser already does. Write a space when e is a symbol.
  - [ ] NEW: OutputForm prints -x^(-2) as -1/x^2 (reference -x^(-2)); InputForm is right
  - [ ] NEW: a bignum base under a symbolic rational exponent that is not an exact root is not split per prime ((3^50 2)^(x/2); the reference gives 2^(x/2) 847288609443^x)
  - [ ] AddScalingParameter-ID6: Simplify finds a strictly simpler form than the reference (simp_split_multiplicative) -- goes with the tie-breaking item
- [ ] With does not rename a Function's bound variable captured from its value: With[{v = 1 + y}, Function[y, {y, v}]][2] is {2, 3} (reference {2, 1 + y})
- [ ] InputForm of a pure function: the reference prints "a + #1 & " (trailing space)
- [ ] purefunc: application still walks the whole body (~2 ms for a 120k-node bound value vs the reference's 12 us); a cached "contains Slot/names" summary would make it O(1)
- [ ] Simplify tie-breaking (FCDiffEqChangeVariables x3, Factor3-ID10)
- [ ] G12 SU(N) (SUNSimplify ID23/83/84, SharedObjectsBehavior-ID201)
- [ ] G8 leftovers: Collect2-ID30, FCSchoutenBruteForce-ID5
- [ ] G13 singletons (see FEYNCALC_REMAINING_GAPS.md)
- [ ] FortranForm/CForm printer (Write2 ID1-3, 5-7)
- [ ] v0.6144 per-dir corpus + timing refresh; raise BASELINE_PASS.txt

# Repository cleanup + README (2026-11-02, docs only, no bump)

Plan: ~/.claude/plans/snappy-sprouting-pie.md

- [x] Move root plans/notes, plans/, reports/, thoughts/ -> docs/notes/; blogs/ -> docs/blog/; rewrite refs (149 files)
- [x] Untrack generic AI-tool scaffolding + worktree gitlinks; remove dead makefile_tests, run_lldb.sh, load.m, *.bak, repo-map.md
- [x] .gitignore reorganised; local PDFs -> ignored papers/; git status clean
- [x] README rewrite (verified example session), CONTRIBUTING.md, docs/building.md restructure, docs/notes/README.md
- [x] Changelog entry
- [x] Verify: build, check-c99/messages/compile-coverage, link check clean; ctest -LE nightly: all pass except dsolve_tests (EigenvalueProblem y'[0]==0,y[Pi]==0 check) -- pre-existing, the diff changes no code

## SeparablePolynomialDecomposition + Simplify stage (v0.6177, 2026-11-05)
- [x] Package `src/internal/simp/SeparablePolynomialDecomposition.m` (user's code; `EqualTo[0]` → `# === 0 &`, as EqualTo is missing; per-row CoefficientRules for 88× on the flattening; zero/decline/inexact rules)
- [x] Inexact: `Rationalize[p, 0]`, else return p (user's rule)
- [x] Simplify stage `try_separable_decomposition` (simp_search.c), eager load in core_init
- [x] Tuned on a 63-case corpus: 1 orientation, ≤5 vars, no product gate (LeafCount −14%, +22% time)
- [x] tests/test_separable_decomposition.c (10 groups), tests/bench_separable_decomposition.c (nightly; verified to fail on the old code)
- [x] Docs (simplification.md), changelog, v0.6177; check-messages, check-c99, leaks clean
- Review: Rationalize defects (drops 1e-300 x beside 1e300; 30-digit Pi → ~1e-17) and Simplify[Sin[p]] > 60 s on the 251-leaf example are existing issues, recorded as follow-ups.

## v0.6178: Rationalize exact core; SeparablePolynomialDecomposition moved to FullSimplify — review
- [x] Rationalize: every mode in exact GMP arithmetic (tiny/huge/subnormal reals, MPFR precision, bignum tolerances, exact bignums). 4 new test groups in test_rationalize.c.
- [x] Re-measured the Simplify stage: +26% on the 63-case polynomial corpus (interleaved, 3x each). Removed the C stage; simp_search.c is back to v0.6176.
- [x] FullSimplify `Separable` transform, memoised per call: 38 leaves on the example, +45% FullSimplify time on polynomials, no measurable cost elsewhere.
- [x] ctest -LE nightly: 591/593 passed in the loaded run. image_tests and parallelmixedtower_tests hit SIGALRM under machine load (load average ~44); both pass alone, and interleaved old/new image_tests timings show no regression. check-messages, check-c99 and leaks: clean.
- Open: N of a near-midpoint subnormal rational; SetPrecision[x, Infinity] gives a Real; [Full]Simplify[Sin[p]] slow.

## Experimental`FeynCalc`DiracTr (v0.6189) — review

- [x] `src/experimental/feyncalc/diractr.c`: native twin of Rolf Mertig's `diracTr` (Woxi #2142), SameQ to it. No evaluator calls: shared `g` nodes, integer matching rows, a row-moving LSD radix sort (2 digits per pass), a bottom-up stamp with `eval_mark_fixed_point_like`. Repeated arguments merge natively. Rules on `g`, arguments or arithmetic fall back to `Expand[...]`, unstamped.
- [x] `eval.c`: `x = rhs` skips the mention walk when rhs is a compound, still-reusable fixed point. It was about 40% of `r = DiracTr @@ Array[x, 16]`. Diffed against the old binary on the self-referential Set cases: identical.
- [x] Measured at n = 14 on a loaded machine (load average 9-21): C 26-40 ms; `.m` 2.9-3.9 s, 9.0 s in the reference kernel; FORM 0.58 s (per the issue). n = 16: C 0.5-0.7 s. Measured before the radix rewrite: 103 ms (n = 14) and 2.6 s (n = 16).
- [x] `tests/test_experimental_diractr.c` (8 groups) passes. `ctest -LE nightly`: 601/601. valgrind leak summary identical with and without the calls. check-c99, check-messages and check-packed-aware are clean.
- [x] `benchmarks/123-native-vs-m` (group G). Found while building it: a HoldRest `benchIf` given the head name through a variable is silently left unevaluated. Also, parsing a statement interns the symbol, so `Names[]` presence must be decided before the symbol is parsed.
- Open: the remaining C-path cost is one Times node and one args array per term, which is irreducible as a boxed Expr; `Cases[_, _g, Infinity]` + `Union` on the result costs 3-5× the build.

## Dense linear algebra: remaining slow cases vs scipy (v0.6221)

Source: the open items of the sparse linear algebra review (complex Cholesky 76x, MatrixFunction[Cos] 23x, Hermitian Eigensystem 20x, exact CharacteristicPolynomial 29x, Adjugate 13x). Re-measured first; each gap traced to a root cause before any change.

- [x] Complex Cholesky: ~95% of the time was boxing the factor into Complex[...] nodes (`na_finish`, a pre-v0.402 restriction). Complex results stay packed now. 0.10 → 0.016 s (scipy 0.006).
  - [x] Found while doing it: packed complex vs boxed SameQ was False when an element has a zero imaginary part. SameQ/UnsameQ were not packed-aware, so the gate materialised Complex[1.,0.]. Put them on AWARE; expr_eq and the hash treat a zero-imaginary element as Real. A first attempt that changed ndarray_buffer_element_to_expr broke packed-complex display (`2. + 0. I`) and was reverted.
  - [x] Found while doing it: `(I/2) * packed` stayed unevaluated (scalar_value rejected Complex with Rational parts).
- [x] MatrixFunction[Cos] 100x100: actually a WRONG ANSWER, unevaluated (benchmark was timing the failure). Higham–Davies truncation bound + unbounded incremental derivative cache. 0.015 s, err ~1e-15.
- [x] Hermitian eigen: packed path used general zgeev (slow, complex eigenvalues with roundoff imaginary parts). zheevd + O(n^2) self-check on packed and boxed paths. Eigenvalues 0.086 → 0.011 s (numpy 0.015); Eigenvectors 0.11 → 0.038 s (numpy 0.042).
- [x] Adjugate: LU det*inv when rcond >= sqrt(u), SVD kept for near singular. 0.0185 → 0.0017 s (numpy 0.0026).
- [x] CharacteristicPolynomial all-rational: FLINT fmpq_mat_charpoly. 40x40: 1.25 → 0.002 s (sympy 0.22).
- [x] Tests in test_packed_list, test_matfunc, test_eigensystem, test_minors_adjugate, test_characteristicpolynomial, test_cholesky.
- Open: CharacteristicPolynomial[m, non-symbol] evaluates where the reference gives General::ivar (pre-existing). Complex Cholesky still 2.7x scipy (load + Hermitian check + copy).

## REPL modernisation PR 1: Ctrl-C to $Aborted, persistent history, Tab completion, colour (started 2026-11-11)

Plan: thoughts/shared/plans/2026-11-11-repl-modernize-pr1.md. Research: thoughts/shared/research/2026-11-11-repl-modernize-pr1.md.
Branch feat/repl-modernize (worktree mathilda-repl). Rules for this PR:
- No version bump, and no version in commit titles. Versions are assigned at rebase.
- Shared files get append-only edits: one appended block in tests/CMakeLists.txt, its own REPL changelog heading, no sym_names edits.
- In core.c, only tc_deadline_passed is touched.
- The src/term/term_interrupt.h API (mth_interrupt_pending / mth_interrupt_request) is fixed for other branches.

- [x] Phase 0: baseline worktree /tmp/mathilda-pr1-base (7e7bc9fb) built; tests/pty/repl_pty_test.py; goldens recorded from baseline; appended CMake block with repl_pty_test (skip 77); harness_can_fail seen failing
- [x] Phase 1: makefile SRC += src/term, LIB_OBJ filters term_readline.o; term_lex (mth_input_complete moved, lexer context); term_lex_tests; nm libmathilda.a has no rl_ refs
- [x] Phase 2: term_interrupt.{c,h} flag API (§4.1); eval poll (loop + stamp fast path); CheckAbort clears; AbortProtect honours; top-level clear; tc_deadline_passed interrupt-aware after the 44-site audit (A/B table in commit body); message gate; Pause sliced; Raylib loops; interrupt_tests (13 cases); kernel suites; perf gate ≤ +1% or noise floor
- [x] Phase 3: term_interrupt handler (SA_RESTART, 3 states, double Ctrl-C within 1 s gives _Exit(130)); term_input seam + stdio backend; term_readline (GNU signal hook, libedit degraded); repl.c single loop and Dialog reader; banner; term_interrupt_tests; pty ctrlc_* cases; USE_READLINE=0 builds; manual libedit check
- [x] Phase 4: term_history (own escaped format, 0600, truncate at load, MATHILDA_HISTFILE / MATHILDA_HISTSIZE); readline feed; Dialog entries; term_history_tests; pty history_* cases
- [x] Phase 5: names_complete; term_complete (symbols, context-aware, options of enclosing head, no evaluate); readline adapter; term_complete_tests; pty complete_* cases
- [x] Phase 6: term_color decision (NO_COLOR, CLICOLOR_FORCE, TERM, per-stream isatty; interactive only); mth_message_set_style; coloured In/Out (\001/\002 on GNU); "took" line ≥ 1 s on a tty; term_color_tests; pty color_* and took_line; check-pipe-protocol
- [x] Phase 7: docs/spec/repl.md + Mathilda_spec.md row; control-flow.md and time-and-date.md (Pause) updates; Pause docstring; changelog 2026-11-09.md REPL heading; SPEC.md §2/§11; building.md; tutorial + book banner; build.yml pty step
- [x] Verify: make clean && make; check-c99/messages/packed-aware/menu-ids/pipe-protocol; USE_READLINE=0 build; ctest -LE nightly (0 FAIL, grep ^FAIL); pty test; valgrind (docker) + leaks vs baseline; reference-name grep on the diff
- [x] Review section: audit table, benchmark medians and noise floor, libedit findings, leak numbers, ctest count, follow-ups for PRs 2-4

### Review (REPL PR 1, 2026-11-11)
Commits: cd687569 (phase 0), eaeb8ecd (1), 6e296319 (2), 5b0e2e5a (3), 86282ac6 (4), 534efef3 (5), f77c2351 (6), a2b899b9 (7), on the 3022bb76 checkpoint.
- **tc_deadline_passed audit:** 36 call sites in 10 files, all class A (clean decline): numloop, reduce_cad (CAD_POLL x17), ludecomp_mpfr, nt_dlog, fi_par, nsum, sparse_la_mpsolve, gb_flint, gb_f4, facpoly_bz_uni. Guarded anyway: root_numeric's session cache (no memoising a declined Root), ludecomp_mpfr's once-per-session warning flags. Added: FactorInteger's split_auto stops climbing tiers on an interrupt (ECM stage 2 never polls, so a Ctrl-C took ~3.5 s; now < 1 s).
- **Performance** (median of 9 interleaved -file runs, baseline 7e7bc9fb vs final branch; A/A noise in brackets): fib[28] -0.35% (-0.21%), Do[f[i], {i, 4*10^6}] +0.44% (+0.50%), Expand[(1+x+y+z)^60] -0.94% (+0.04%), 12 Integrate calls -1.51% (-0.19%). At phase 2: fib +0.91% (A/A +0.68%). Gate (<= +1%) holds; the poll is not throttled.
- **ctest -LE nightly:** branch 645/645 pass, 0 `^FAIL`; baseline 638/638 (the 7 extra are the new suites and repl_pty_test). No pre-existing failures.
- **pty test:** 32/32 on the branch. Invariance cases (pipe, pipe with tty stderr, -file, -file on a tty) byte-identical to goldens recorded from the baseline.
- **Byte-for-byte corpus:** 4086 `In[n]:=` examples from docs/spec/builtins through pipe mode: stdout identical (859 MB); stderr identical except 17 lines of `[x.xx ms]` timing output. 3953 of them as one -file script: stdout (284 MB) and stderr identical after the same timing normalisation.
- **Build:** make clean && make: the same 6 warnings as the baseline clean build (FLINT headers, reduce_companions), none new. check-c99, check-messages, check-packed-aware, check-menu-ids, check-pipe-protocol pass. libmathilda.a: 0 undefined rl_ symbols; T mth_interrupt_pending / mth_interrupt_request present. USE_READLINE=0 (scratch copy): warning-free, no readline linked, baseline fgets banner, Ctrl-C ends the process as before.
- **Leaks (macOS leaks; Docker daemon not running, so no valgrind):** term_lex, term_history, term_complete, term_color, term_interrupt, interrupt suites: 0 leaks each. Interactive session (common script) baseline 0 / branch 0; full script with two Ctrl-C aborts, a prompt Ctrl-C, completion, history recall and the took line: 0 leaks.
- **libedit** (make HOMEBREW= MACPORTS=, links /usr/lib/libedit.3.dylib): Ctrl-C during Do, While, Pause, FactorInteger, CheckAbort, AbortProtect, TimeConstrained and Dialog gives $Aborted; double Ctrl-C exits 130; Ctrl-C at the prompt is ignored (no exit, no submit, line kept); Symbol and unique option completion work; Tab-Tab does not list; prompts are uncoloured; Out and messages are coloured. Multi-line input does not work on libedit (pre-existing: the RET rebinding is ignored there), so the multi-line history and continuation cases fail there as before. The "bind: Invalid command enable-bracketed-paste" line on libedit is pre-existing.
- **Follow-ups:** classic Buchberger fallback after FLINT declines, and GMP-ECM stage 2, ignore the interrupt (finish their step; double Ctrl-C works). Ctrl-C during `$Post`/printing aborts the print evaluation only. PR 2: in-tree editor (would give libedit users line clearing and multi-line), token lexer/highlighting, replace find_unterminated / mathilda_ffi_is_complete with term_lex. PR 3: status line/spinner. PR 4: init.m, themes ($-symbols), slash commands. Pipe-mode (notebook) interrupt is still out of scope.

## REPL PR 2: the input editor (feat/repl-editor, stacked on feat/repl-modernize)
Plan: thoughts/shared/plans/2026-11-11-repl-modernize-pr2.md. Research: thoughts/shared/research/2026-11-11-repl-modernize-pr2.md.
Rules: no version bump; append-only edits to tests/CMakeLists.txt, the changelog and docs/spec/repl.md; term_color.{c,h} untouched; repl.c edits limited to the reader selection, the banner and one result-print line.

- [x] Phase 1: term_width, term_highlight (incremental per-line lexer, symbol classes, styles); unit tests; CMake block
- [x] Phase 2: term_render (grid diff, one write), term_editor (multi-line, keymap, history edges, smart/forced Return, kill/yank, undo, paste, brackets), term_editor_tty (raw mode, restore paths, SIGWINCH, DSR, MATHILDA_EDITOR); repl.c wiring; unit + pty tests; milestone 1
- [x] Phase 3: completion popup (Tab, filter, docstring row, options first, placement, degradation); tests
- [x] Phase 4: Ctrl-R history list; tests
- [x] Phase 5: result colouring (expr_print_to_string, safe subset); tests
- [x] Phase 6: latency gate test (200 lines, < 2 ms, one write, 0 allocs); resize / narrow / UTF-8 / restore pty cases
- [x] Phase 7: docs/spec/repl.md, changelog, SPEC.md
- [x] Verify: clean build, checks, USE_READLINE=0, ctest -LE nightly, pty suite, byte-identical pipe/-file, leaks, latency numbers, readline fallback; /tmp/pr2-body.md

### Review (REPL PR 2, 2026-11-11)
Commits: c30f8472 (phase 1), 45bb2573 (2), 52e9fe72 (3), 6fdc5a2a (4), 226139ee (5), 7d04dbee (6), f1aff3e4 (7), on PR 1's a0c27991.
- **Build:** make clean && make: only the pre-existing FLINT acb.h -Wstringop-overread warnings; none from the new code. check-c99, check-messages, check-packed-aware, check-menu-ids, check-pipe-protocol pass. USE_READLINE=0 (scratch copy, -j4): warning-free, no readline linked, invariance pty cases pass, stdio banner as PR 1.
- **ctest -LE nightly:** 648/648 pass (PR 1: 645 + term_highlight, term_editor, term_editor_latency).
- **pty:** 48/48 (PR 1's 32, two adjusted only for readline specifics: no_sgr checks SGR only; the second-Tab listing is the popup). Readline fallback (MATHILDA_EDITOR=readline): smart Return, Tab, ^C as PR 1.
- **Byte-identical vs PR 1 (a0c27991):** pipe mode over 5157 doc examples: 5151 responses identical, stderr identical; the 6 that differ are RandomReal examples that also differ PR 1 vs PR 1. -file over 4989 parsable examples wrapped in Print (5026 stdout lines, 152451 stderr lines): stderr identical, stdout identical but for one RandomReal line that also differs PR 1 vs PR 1. Goldens (pipe, pipe+tty stderr, -file, -file on a tty) identical.
- **Latency gate** (200-line buffer, 120x40, colour, 3000 keys): mean 22.6 us, p50 25 us, p99 29 us, max 34 us; whole-buffer relex 69 us; 0 allocations; never more than one write per key. End-to-end through a pty: p50 0.05-0.07 ms, p95 0.10 ms, max 0.22 ms.
- **Leaks (macOS leaks):** term_highlight, term_editor, term_editor_latency suites 0; interactive session (multi-line, popup, Ctrl-R, undo, paste, Ctrl-C, resize, Dialog, coloured results) 0.
- **Follow-ups:** Up/Down move by logical line, not by wrapped row; a resize on terminals that reflow may leave a stale row above the area; Dialog reads during an evaluation need PR 3's status line paused (their side of the hand-over).

## REPL PR 3: spinner, status line, progress (branch feat/repl-status)
Plan: thoughts/shared/plans/2026-11-11-repl-modernize-pr3.md. Research: thoughts/shared/research/2026-11-11-repl-modernize-pr3.md.
Stacked on feat/repl-modernize (a0c27991). Shared files append-only: own CMake block, own changelog heading, own repl.md sections; term_color append-only; term_interrupt one counter.
- [x] Phase 1: term_progress stack; hooks in Do, Table, Sum, Product, iter_run (Nest/NestList/Fold/FoldList/NestWhile/FixedPoint); worker counter in five parallel drivers; term_progress_tests; perf check
- [x] Phase 2: term_spinner layout (fit/drop order, colour/plain/ASCII, bar, ETA, took and status lines); TermTone appended to term_color; term_spinner_tests
- [x] Phase 3: term_status forwarder thread (pipe for fd 1/2, ECHO off, 200 ms delay, 10 Hz, col-0 rule), Ctrl-C ^C/quiet period, double-tap exit fn, fatal-signal and atexit cleanup, Dialog suspend/resume, took/status lines in repl.c; term_status_tests; 15 new pty cases
- [x] Phase 4: docs/spec/repl.md sections, changelog heading, Mathilda_spec.md row
- [x] Verify: make clean && make; check-*; USE_READLINE=0; ctest -LE nightly; pty; pipe/-file invariance vs PR 1; leaks; perf gate; TSan driver

### Review (REPL PR 3, 2026-11-11)
Commits: b6abee2e (progress stack), e5ff0bb9 (layout), 985c4f3e (forwarder + REPL), a4697067 (docs), 72db8cb6 (check-c99 guard).
- **Build:** make clean && make: the 6 FLINT header warnings only, none new. USE_READLINE=0 scratch build: same warnings, no readline linked, fgets loop unchanged, no spinner.
- **check-*:** check-c99 (after 72db8cb6), check-messages, check-packed-aware, check-menu-ids, check-pipe-protocol pass.
- **ctest -LE nightly:** 648/648 (PR 1: 645 + term_progress/term_spinner/term_status), 0 FAIL. Built with -DCMAKE_C_FLAGS=-O2 -j4 for disk space.
- **pty:** 47/47 (32 PR 1 + 15 new).
- **Invariance vs PR 1 (a0c27991):** 500-request NDJSON and 600-statement -file corpora: stdout and stderr byte-identical, no ESC bytes; pty goldens from 7e7bc9fb pass.
- **Leaks (macOS leaks):** interactive session with spinner/Print/progress/Dialog/Ctrl-C: 0 (PR 1: 0). term_spinner_tests 0, term_status_tests 0, term_progress_tests 7 nodes / 480 B: the deliberate TimeConstrained timeout through nested loops (documented siglongjmp leak) plus one 64 B node in a Nest case that does not reproduce in isolation or in the binary.
- **Perf (median of 9 interleaved -file runs, PR 1 vs branch; A/A in brackets):** fib +0.87% (+0.32%), Do 4M -0.33% (+0.59%), Expand -0.08% (+1.69%), Integrate x12 +1.04% (-0.93%). Loop extras: FoldList 1.5M -0.34% (+1.22%), NestList 1.5M +1.07% (+1.06%), 1000 Nest calls -0.04% (-1.71%).
- **TSan** (Apple clang, driver over term_status/spinner/progress with stubs): only the 4 by-design racy reads (progress depth/slots, worker count, dispatch head); 0 with those suppressed.
- **Deviations:** no MemoryConstrained exists (memory share is of RAM; hook provided); Map not tracked (does not evaluate its elements); status line off by default; NO_READLINE builds get no spinner.
- **Rough edges / follow-ups:** Run[] children see a non-tty stdout during evaluation; Ctrl-Z mid-evaluation leaves the line and ECHO-off until resume/shell reset; SIGKILL cannot restore ECHO; the live head can flicker between builtins; Linux not run locally (CI runs the pty suite); FoldList[g, x, Range[300000]] segfaults (deep nesting) in PR 1 and here alike.

## REPL PR 4: slash commands, themes, mods, Monitor; PR 2 + PR 3 merged (branch feat/repl-menus)
Stacked on feat/repl-editor (PR 2) with feat/repl-status (PR 3, 8c9ac092) merged in.
- [x] Merge PR 3 into PR 2: conflicts in Mathilda_spec.md, tests/CMakeLists.txt, changelog, todo, repl.md, pty harness; editor/spinner hand-over verified; 6 combo_* pty cases
- [x] Slash-command menu (leading '/' on an empty buffer only; never evaluated): /help /theme /history /timing /statusline /clear /mods /quit + mod commands
- [x] Themes dark/light/high-contrast/mono/ansi; 24-bit/256/16/none from COLORTERM/TERM; NO_COLOR wins; /theme picker with live preview (Return applies and saves, Esc reverts); ~/.mathilda/config; MATHILDA_THEME
- [x] Mods: ~/.mathilda/init.m + mods/*.m; SlashCommand[]; $PreRead/$Pre/$Post/$PrePrint; /mods; failures reported once
- [x] Monitor[expr, mon] in the status line
- [x] Polish: Ctrl-Z at the prompt and under the spinner; no cursor-query wait; Up/Down by wrapped row; redo (Alt-_); banner; docs/spec/repl.md rewritten
- [x] Verify (below)

### Review (REPL PR 4, 2026-11-12)
Commits: 097dd288 (merge), f9073acb (PR 4 core), bb8c6ba4 (tests), ad060249 (docs/polish), 2eacd69b (dark palette below 24-bit, latency gate with the thread).
- **Build:** make clean && make: only the pre-existing FLINT acb.h warnings. USE_READLINE=0 scratch build: same, no readline linked, stdio loop unchanged ("/help" is a parse error there, as before).
- **check-*:** c99, messages, packed-aware, menu-ids, pipe-protocol pass.
- **ctest -LE nightly:** 652/652 (PR 1 645 + PR 2's 3 + PR 3's 3 + term_menus_tests), 0 failures.
- **pty:** 83/83 (PR 1 32, PR 2 16, PR 3 15, 6 combo, 14 PR 4).
- **Invariance vs PR 1 (a0c27991), doc-example corpus, streamed through sha256:** pipe stdout identical (8251 lines, c94f4ce5...); -file stdout identical (3999 lines, 7aebd3f0...); no ESC bytes. Stderr (157602 / 157554 lines) differs on 11 lines in each mode, all `[x ms]` timing lines; PR 1 vs PR 1 differs on 18 / 13.
- **Leaks (macOS):** term_menus, term_editor, term_highlight, term_spinner, term_status: 0. term_progress: 7 nodes / 416 B, the documented TimeConstrained siglongjmp case from PR 3. Interactive session (slash menu, theme picker apply/cancel, mods incl. failing ones, mod commands, Monitor, Dialog, Ctrl-C, resize): 0, via `leaks <pid>` on the live process.
- **Latency:** term_editor_latency_tests p50 26 us, p99 35 us, max 82 us, 0 allocations. Pty e2e with the status thread started: p50 0.04 ms, p95 0.05 ms, max 0.12 ms.
- **Perf (-file, medians, PR 1 vs branch, A/A in brackets):** 9 rounds: fib +1.22% (+0.71%), Do 4M +0.45% (-0.54%), Expand +1.80% (+3.75%), Integrate -1.82% (-0.49%). 15 rounds: fib +0.02% (-0.54%), Do +1.70% (+0.00%), Expand -1.33% (-2.45%), Integrate -3.35% (-2.29%). Focused Do 4M, 15 rounds: branch +0.86%, PR 3 +0.12%, A/A +0.74%.
- **Deviations:** Ctrl-Z pty case runs under a job-control wrapper (a pty child's group is orphaned); dark uses the terminal palette below 24-bit; commands go into history; Monitor updates on loop steps and the 128-step evaluator tick; key-hint footer and input box not done.
- **Rough edges / follow-ups:** no OSC 11 background detection; Monitor text is shown in InputForm (strings quoted); a late cursor-position reply after very fast typeahead might echo; the spinner clears with EL before each frame.

## Bronstein completion campaign (2026-11-14) — one simple radical, full decision procedure

Plan: `src/internal/bronstein/BRONSTEIN_PLAN.md` § "2026-11-14 — Completion campaign".

- [x] Phase 0: soundness — (a) rt_field_integrate certification scopes + Cos[E^x] surface-arg bug fixed, (b) structure theorem tri-state, (c) refuted (Def 2.21 allows logs of any k-element) — v0.6545
- [x] Phase 1: failure taxonomy, Integrate::incomplete (Method + cascade), completeness gate — v0.6545
- [ ] Phase 2: C substrate
  - [x] WeakNormalizer over k(t); ParametricLogarithmicDerivative (exact); RdeSpecialDenomExp / RdeBoundDegreeExp / RdeBoundDegreePrim(db=da) cancellation; PolyRischDECancel log-derivative branch; dependent-tower detection — v0.6546
  - [ ] rewrite dependent Log/Exp monomials into an independent tower
  - [ ] ParamRischDE (Ch. 7) + LimitedIntegrate with new logarithms (removes "limited integration undecided", "log-derivative branch integral undecided")
  - [ ] rt_field_integrate completeness (m > 1 new logs, mixed resultants); constant fields in C
- [ ] Phase 3: Shape A completion (.m)
  - [x] front end: cyclic compositum as one radical (E6b), E^(r Log a) + log-linear exponent split + monomial pruning (E6c), reducible logs via kappa (E6a); fix C dependent-log retry wrong answer (Log[x^2]/x) — v0.6547
  - [x] compositum up to a locally constant sign (E6d) — v0.6549
  - [ ] multi-kappa / kappa-in-denominator certificates; units of order > 2 (Phase 4)
  - [ ] C: Log[x^2] as its own monomial (RischTranscendental leaves Int Log[x^2]/x unevaluated now)
  - [ ] PRE-EXISTING (found by the E6a review, also on main): PossibleZeroQ[Log[E^x] - x] is True (false zero; -2 Pi I for Im x in (Pi, 3 Pi]) so Integrate[(Log[E^x]-x) Sqrt[Log[x]] + 1/x, x] = Log[x] (wrong); Integrate[(Log[x^2]-2Log[x]) Exp[x^2], x] = Indeterminate (Automatic cascade; Method->Bronstein certifies); Integrate[Exp[x + Log[x^2]/4], x] -> uncaught Throw[..., "PIM"]
  - [x] #9 genus >= 1 special places (E3: specialization + two primes with special places; goodQ hook) — v0.6548
  - [ ] perf: three-radical compositum D[x^2 Sqrt[x] (x+1)^(1/5) (x+2)^(1/3)] (Kummer Z/30) > 300 s; Sqrt[x](x+1)^(1/5)(x+2)^(1/3) > 120 s
  - [ ] perf: D[Log[y+u] + 2 Log[y^2+u y+u^2]] on y^3 = u^3+u+1 takes 225 s (profile the torsion/principal steps)
  - [x] #12 weights: theorem E6e (they vanish), Internal invariant; exp special places keep the site — v0.6550
  - [ ] perf: 13 of 120 random m = 3,4 components time out at 30 s
  - [ ] brParamIntegrate -> LimitedIntegrate; LogPart:171; Divisors multi-ramified (E5); Tower:logrel/exprel robustness
- **Known open gaps: the machine-checked register `tests/bronstein_paper/open_gaps.m`** (Part D of
  docs/notes/bronstein_extensions.tex; nightly ctest bronstein_open_gaps_stress). Every newly found gap
  goes there with an example the same day; closing one moves it to examples.m.
- [x] Adversarial reviews of Bronstein / Trager / PMT (2026-11-14): false certificates in all three fixed (v0.6555); open findings in open_gaps.m
- [ ] Phase 4: constant fields (NF, Q(i)+trig, parameters)
  - [x] 4.1 Q(i) with a literal I: tragerFactorList, residue fields containing i (fieldSubI), mod-p refusal, Q(i) coefficient tests — v0.6551
  - [x] 4.1b I mod p (gfMod, goodPrime consistent root), Gaussian residue polynomials, Gaussian kappa — v0.6554
  - [x] 4.2 trig front end (TrigToExp) — v0.6554
  - [ ] presentation: Root objects in Sqrt[Tan[x]]-type answers -> radicals / real form; Trager:residue-field-QI (tower method over Q(i))
  - [ ] 4.3 general number fields (fieldFn[vars, {alpha, mu, theta}]); 4.4 symbolic parameters (E4)
  - [ ] flaky under full parallel ctest load (seen 3x since v0.6548, never alone or in a fresh loaded process): integrate_cascade_trager_bronstein_tests:89, Check[Integrate[5/((x - 3) Sqrt[x^3 - 2]), x]; "ok", "msg"] gets "msg" in 1.9 s -- an order/state-dependent message in the same process, not a timeout; find the message name (wrap with a $MessageList print in the suite)
- [ ] Phase 5: radical anywhere (brDField, transcendental Risch over algebraic base, RDE on a curve over a tower)
- [ ] Theory note docs/notes/bronstein_extensions.tex (E1–E6; LaTeX per user 2026-11-14), each adversarially reviewed — E6a–c written v0.6547, review pending
- [ ] Close-out
