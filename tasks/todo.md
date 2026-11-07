# SIAM 100-Digit Challenge benchmark + holdouts + fixes (started 2026-11-07)

Plan: ~/.claude/plans/pasted-content-id-9fcd-similar-to-fluttering-bentley.md
Register: benchmarks/SIAM100.md. Scorer: benchmarks/siam100_table.py.
Target: on every case, Mathilda's correct digits >= the best of WL / scipy+mpmath /
Maxima / FriCAS, and time < 1.5x the fastest correct system; machine + wp50 tiers.

- [x] Phase 0: harness -- FriCAS column in run_all.py (.input halves, VALUE lines kept),
      siam100_common.{m,py,mac,input}, siam100_cap.lisp, siam100_table.py
- [x] C1 (exp 166) encoded in all five systems, reference by two methods
- [x] Phase 1-2: encode C2-C10, D1-D20 (exps 167-195) with references (weak: C3 40d, D8 5d, D14 12d, D15 8d, D18 7d, D20 7d)
- [ ] Baseline run of 166-195, status table in SIAM100.md
- [ ] Phase 3: holdout A, 10 per problem (exps 196-225), run once BEFORE fixes
- [ ] Phase 4: sealed holdout B, 10 per problem, branch siam100-holdout-b, baselines only
- [ ] Phase 5: fix waves (W wrong / X absent / D fewer digits / S slower)
- [ ] Phase 6: final run, holdout B once, verdict

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
- [ ] Phase 3: benchmark 120 broadened, 3-way, every row <= 1.0x
- [ ] Phase 4: docs, valgrind on the fuzz corpus, `make check-sparse`
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
