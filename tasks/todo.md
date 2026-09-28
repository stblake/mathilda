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

## Commit 2 — confirmed correctness bugs (bump+tag)
- [ ] power.c (6), plus.c (2), times.c (1) → mth_message_gated(g_arith_warnings_muted,...)
- [ ] linalg/matpow.c (3) → expr_to_string + mth_message
- [ ] solve/solvenlsys.c:139, solve/solveinv.c:149
- [ ] findmin_common.c fm_warn note

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
