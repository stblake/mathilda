# Loop performance: Do / For / Nest → close gap to Mathematica

Plan: /Users/user/.claude/plans/here-are-updated-timings-golden-meadow.md

## Stage 1 — Targeted cleanups (safe, shippable) — DONE ✓
- [x] 1a. Scalar Nest/Fold/FixedPoint: bounded history window (NestWhile excluded — uses buf->count)
- [x] 1b-i. Bucketed args-array free-list in expr.c (expr_new_function/expr_free/expr_unshare)
- [x] 1b-ii. Plus/Times final_args stack buffer (plus.c, times.c)
- [x] 1b-iii. numeric_contagion_args caller-buffer (numeric.c/.h + plus/times callers)
- [x] Stage 1 verify: Do 1.16→0.78, Nest 1.48→0.88, For 2.34→1.68; all Out identical;
      tests pass (expr_pool/expr/eval/funcprog/fold/purefunc/numeric/expand); bench_assoc PASS;
      valgrind clean (only macOS dyld/objc baseline, 0 Mathilda frames)

## Stage 2 — Automatic numeric loop fast-path — DONE ✓
- [x] 2.0 numloop.{c,h}: RPN compiler + double VM; makefile auto-discovers; CMake COMMON_SRC
- [x] 2a. Wire Nest (nest_impl) + Do (builtin_do count form)
- [x] 2b. Wire For (builtin_for) + While (builtin_while)
- [x] 2c. Extend to Fold, FixedPoint, NestWhile (+ bare-head f like Nest[Cos,...]) per user request
- [x] 2d. Compound multi-statement bodies (CompoundExpression of Sets) via NumBlock — shared
      register file, sequential statements, temp/multi-var; + Do integer range form. Per user
      follow-up (`Do[a; b; c, {n}]` was falling back). `Do[…;…;…, {10^6}]` 1.85s→0.027s.
- [x] tests/test_numloop.c: differential (fast vs interp, 1e-9 rel) + fallback correctness (33 tests
      incl. compound/multivar/range/temp-var)
- [x] Stage 2 verify: Do/Nest/For ~0.013s (≈100x over Mathematica); compound Do/For/While fast;
      exact/symbolic/read-before-assign/non-Set untouched; valgrind clean (0 numloop frames);
      bench_assoc PASS; funcprog/fold/eval/cond suites PASS; Table/Sum/nested-Do sane

## Docs
- [x] docs/spec/changelog/2026-07-20.md perf entry (both stages)

## Result summary
| Benchmark (10^6) | before | after | Mathematica |
|---|---|---|---|
| Do logistic  | 1.16s | 0.013s | 0.78s |
| Nest logistic| 1.48s | 0.014s | 0.029s |
| For logistic | 2.34s | 0.014s | 1.30s |
Fast path is transparent + gated on machine-real + inexact-result; falls back to
interpreter otherwise. Not bit-identical (Orderless runtime-value sort) — agrees
to FP rounding.

## Review
Correctness gate is the crux: fast path only fires when the interpreter's result
is provably an inexact machine real, so exact/symbolic loops are byte-identical
to before. Non-finite intermediates bail to the interpreter. Stage 1 (pool/window
cleanups) is independently valuable and shippable even without Stage 2.
