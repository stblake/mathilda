# Route every user-facing message through the Quiet/Check funnel

Plan: `/Users/user/.claude/plans/we-need-to-do-misty-boot.md`

## Commit 0 — tooling (no version bump)  ✅ DONE (22bbd3da)
- [x] `tools/check_message_routing.py` — multi-line-aware detection, EXEMPT + BASELINE, ratchet
- [x] Seed BASELINE with the full current backlog (294)
- [x] `make check-messages` target + `.PHONY`
- [x] Wire into `.github/workflows/build.yml`
- [x] Verify: gate green on unchanged tree

## Commit 1 — funnel + highest leverage (bump+tag v0.222)  ✅ DONE
- [x] `mth_message` / `mth_message_gated` / `mth_message_v` / `mth_message_cont` in message.h/.c (guarded printf attr)
- [x] Migrate `common.c:builtin_arg_error` (4 fprintf) → dropped common.c from BASELINE (290 left)
- [x] Build clean, check-messages ratchets, check-c99, behavior verified (Check[Fourier[],CAUGHT]→CAUGHT, Quiet suppresses)

## Commit 2 — confirmed correctness bugs (bump+tag v0.223)  ✅ DONE
- [x] power.c (6), plus.c (2), times.c (1) → mth_message_gated(g_arith_warnings_muted,...)
- [x] linalg/matpow.c (3) → expr_to_string + mth_message
- [x] solve/solvenlsys.c warn_nsdim, solve/solveinv.c emit_ifun (note; ifun-suppress stays moot)
- [x] findmin_common.c fm_warn → mth_message_v(g_fm_quiet,...) (adds note)
- [x] Verified: Check[MatrixPower[..,1/2]]→CAUGHT, Check[Power[0,-2]]→CAUGHT, Quiet suppresses, Quiet[Check]→FAILED. BASELINE 290→279.

## Commit 3 — subsystem helpers → wrappers (bump+tag)
- [ ] fit_warn, fm_warn, DRY re-body fs_msg/dt_msg/ops_msg/root_warn/inv_warn/purefunc

## Commits 4..N — inline sites by module (bump+tag each)
- [ ] solve/
- [ ] linalg/
- [ ] poly/
- [ ] calculus/ (incl. Integrate::nonelem + RischTranscendental routed)
- [ ] numerical_calculus/
- [ ] numerical_roots/
- [ ] strings/ + strings/regex/
- [ ] product/
- [ ] top-level A: int.c, real.c, special_functions/, numbertheory/
- [ ] top-level B: eval/core/context/refine/rootreduce/numberform/funcprog/names/interp/complex_expand/options_builtin/vectors/vectoranal/partitions/bitwise/list/precision/numeric/random/nc_accuracy/expand*

## Final commit — assert-empty + docs (bump+tag)
- [ ] BASELINE == {}; flip gate to assert-empty
- [ ] tests: test_messages.c (or extend test_parallelmixedtower); update deliberate-raw test
- [ ] docs: CLAUDE.md, SPEC.md §9, docs/design/message_routing.md, changelog
- [ ] valgrind spot-check; rebuild code-review graph

## Review
(to be filled in)

## `Integrate` — `ParallelMixedTower`: the real form of the logarithmic part (logrewrite.m, 2026-09-28)

Port of the research campaign of 2026-09-28 (logrewrite.{py,wl,mac}: Rioboo's collapse of conjugate
pairs of logarithms to real logarithms, arctangents and hyperbolic arctangents) to the Mathilda
package, with pre/post measurements. `.m` files, one C test and docs only; no C change.

### Plan
- [x] 1. Probe every kernel construct of logrewrite.wl in Mathilda first (three probe scripts): FreeQ on
  complex atoms, ComplexExpand on nested radicals, CountRoots, $InputFileName, nested TimeConstrained,
  Do-iterator capture, PolynomialExtendedGCD / PolynomialRemainder over algebraic constants, D and N of
  ArcTan/ArcTanh forms for the verify gate
- [x] 2. Baseline (pre) on the unchanged package: Charlwood 50 (charlwood_wl.py mathilda), review corpus
  (review_wl.py, SYSTEM=mathilda), tests/build/parallelmixedtower_tests
- [x] 3. `src/internal/mixed/logrewrite.m`: one function per function of logrewrite.wl, loaded with
  `LoadModule["mixed/logrewrite.m"]` inside the private context; host adaptations (`_Complex`,
  `RectPow`, `SturmCount`, `RootRadicals`, `CanRaw`, lr-prefixed iterators); default `LogToReal` in
  ParallelMixed.m when the file is missing; hook at the end of iPIM (`lrTerms` -> `LogToReal`)
- [x] 4. Develop against a copy of src/internal selected with MATHILDA_HOME; compare 42 integrands with
  Mathematica (trace lines and derivative checks at points including x = 2, 3, 5/2)
- [x] 5. Fix found in every port (a defect of the campaign, not of the port): the arctangent argument over
  the radical is a polynomial in y (YQuot / _y_quot / pm_y_quot) -- a radical in a denominator was
  read on the principal branch by the surface layer (Charlwood P5 wrong for Tan[x] < 0 in Mathematica
  and Mathilda); re-verified in the four ports
- [x] 6. Install, C test `test_method_real_form`, MATHILDA_DIVERGENCES.md A19-A24 / B7 / D, changelog
- [x] 7. Post measurements on the installed package (Charlwood, review corpus, C tests) and the review
  section below

### Review
- Both measurements on build 0.222 on a quiet machine, the baseline through `MATHILDA_HOME` pointing at
  a copy of the previous module tree (the binary was rebuilt from 0.221 to 0.222 during the session):
  Charlwood 48/50 -> 49/50 (A19 passes the gate in its real form), kernel 20.4 s -> 23.6 s, median
  0.209 -> 0.204 s, results with I 16 -> 0; review 303/68/0 -> 303/68/0, 693 s -> 602 s, 44 -> 0
  returned integrals with I; `parallelmixedtower_tests` 11 tests pass incl. `test_method_real_form`.
- The port exposed a defect of the campaign itself (the arctangent argument with the radical in a
  denominator, P5 wrong for Tan[x] < 0 in Mathematica too), fixed in all four ports and re-verified.
- Six new host divergences (A19-A24) and one behavioural one (B7), each with a one-line repro and a
  workaround in `logrewrite.m`; A24 (`Can`'s field detour returning `Dot[{}, Inverse[{}], {}]`) has
  no standalone repro yet -- it appeared on A40 after P4 in one kernel and is worth a core look.
- Not done: the `E^x` / `E^(2 x)` independent-generator defect of BuildTower (both ports certify
  `1/(1 + E^x + E^(2 x))` as not elementary) was noticed and left alone; the rewrite's 30 s budget is
  not interruptible by the package budget (B7). Nothing committed; the user was committing their own
  work in this repository at the time.
