# M29 — DSolve §2.2.9 corpus (Problems 801–900) + piecewise differentiation of rounding functions

## Plan
Two intersecting parts (they meet at problem 898 `y''+9y=2 Sec[3x]`):
- **Part A**: add §2.2.9 (Problems 801–900) to the DSolve corpus (baseline 0 — 100/100
  out of the box, measured).
- **Part B**: piecewise derivatives of Floor/Ceiling/Round/IntegerPart/FractionalPart
  (user request; upgrades 898 from UNK-trusted to genuinely verified).

## Progress
- [x] Empirical: §2.2.9 = 100 scalar problems, converter handles unchanged, 100/100 PASS.
- [x] Part B: 5 derivative handlers in `src/calculus/deriv.c` (after UnitStep block).
- [x] Part B: rebuilt; all 5 outputs + chain rule match Mathematica exactly.
- [x] Part B: tests in `tests/test_deriv.c` (`test_rounding_deriv`) — pass.
- [x] Part A: generated `DSolve_test_status/DE_examples_229.m` (curl real page + --url).
- [x] Part A: registered `dsolve_corpus_2_2_9_tests` in `tests/CMakeLists.txt` (baseline 0).
- [x] Part A: reports/2.2.9.{tsv,md}, STATUS.md block + wave bullet, README row.
- [x] Part A: regression — regenerated §2.2.1–§2.2.8, records byte-for-byte identical.
- [x] `t_m29_sec_floor_verifies` in `tests/test_dsolve.c` (898 residual numericizes) — pass.
- [x] Bookkeeping: DSOLVE_PLAN.md M29 entry, docs/spec/builtins/calculus.md, changelog 2026-09-07.md.
- [x] Verify: §2.2.1–§2.2.9 gates pass; §2.1.2 gate passes; deriv/piecewise/nderiv units pass;
      dsolve_tests pass; check-c99 clean; leaks 0 on deriv path.
- [~] dsolve stress suites (running).

## Review
**Outcome.** §2.2.9 (Problems 801–900) added to the DSolve corpus and solves **100/100** out of
the box (baseline 0) — dominated by 2nd-order linear const-coeff (homogeneous + nonhomogeneous),
Euler/Emden–Fowler, complex-coefficient, and x(t)-dependent-variable forms, all already covered
by existing specialists. No ODE-solver fix was needed.

**Engine feature (user request).** Added piecewise differentiation of the integer-rounding
functions in `src/calculus/deriv.c` (5 handlers after the UnitStep block), returning the exact
Mathematica `Piecewise[{{v, cond}}, Indeterminate]` forms and composing with the chain rule.
This upgraded corpus 898 (`y''+9y==2 Sec[3x]`, a VoP solution carrying a Floor branch-tracking
term) from a non-numericizable "trust DSolve" pass to a genuine numeric residual (~1e-38).

**Verification.** 9× §2.2.x gates pass; §2.1.2 (1204 records) passes; no regression anywhere.
deriv/deriv_array/deriv_symbolic_order/nderiv/piecewise/dsolve unit suites pass; check-c99 clean;
0 leaks on the new derivative path (macOS `leaks`).

**No new builtins/attributes/symbols** — all heads (Piecewise/NotElement/Element/Integers/Re/Im)
and SYM_ constants already existed. Change set: deriv.c + 2 test files + CMakeLists + corpus file
+ 2 reports + STATUS/README/PLAN/spec/changelog docs.
