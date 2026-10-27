# GitHub issues #112 / #113 / #114 (2026-10-27)

- [x] #113 Together on a denominator-free sum: expand (and re-evaluate) the combined numerator when lcm_den = 1 ;
      Cancel's embedded-rational gate no longer walks into exponents (y^(1/(a-1)) is one generator) ; v0.6050
- [x] #114 Factor over Q(I): factor_degree_one concludes irreducible when gcd(L, C) = 1 (indeterminate vars only) ;
      Gaussian FLINT zero test (I -> tau, deg_tau <= 1) ; 63 s -> 26 ms ; v0.6051
- [x] #112 rb_rec shares unchanged subtrees ; Expand estimates once (child bounded) + hash-screened kernel table ;
      Union hash-dedupe then sort ; Cases `_h` head test ; Echo + EchoTiming ; v0.6052
      n=14: generation 24.1 -> 18.8 s, Cases+Union 3.62 -> 0.88 s, MaxMemoryUsed 462 -> 224 MB (same machine, A/B)
- [ ] #112 follow-up: FLINT expansion with opaque kernels as generators (exclude List/Interval/SeriesData/...)
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
