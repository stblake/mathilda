---
schema_version: 2
title: UnitBox builtin
slug: unitbox
status: in-review
source: direct-user-request
owner: Michael Sollami
issue: pending
pull_request: https://github.com/stblake/mathilda/pull/67
started: 2026-08-18
last_updated: 2026-08-18
blocked_by: none
goal_lock:
  status: active
  stamped: 2026-08-18 14:51
  scope:
    - "src/piecewise.c"
    - "src/piecewise.h"
    - "src/info.c"
    - "docs/spec/builtins/elementary-functions.md"
    - "docs/spec/changelog/2026-08-17.md"
    - "tests/test_unitbox.c"
    - "tests/CMakeLists.txt"
  success_criteria:
    - "AC-1 -- UnitBox[x] returns 1 for -1/2 <= x <= 1/2 (closed) and 0 outside, for exact/Real/MPFR/symbolic-real arguments"
    - "AC-2 -- UnitBox is registered with Listable, NumericFunction, Protected (NOT Orderless)"
    - "AC-3 -- UnitBox has a symtab_set_docstring entry"
    - "AC-4 -- docs/spec/builtins/elementary-functions.md has a UnitBox section before UnitStep"
    - "AC-5 -- current week's changelog documents the addition, including the two explicitly-deferred/accepted tradeoffs"
    - "AC-6 -- tests/test_unitbox.c exists and is registered with both add_executable AND add_test()"
---

# UnitBox builtin

This file is the lifecycle record for a single feature. It is created by `/tracking` and updated through the spine (`/dev-research`, `/dev-plan`, `/dev-implement`, `/dev-validate`, `/dev-reflect`).

**This file is sealed once `closeout complete` is stamped. Do not edit it after that point.**

## Tracking Metadata

- Title: `UnitBox builtin`
- Slug: `unitbox`
- Status: `in-review`
- Source: `direct-user-request`
- Owner: `Michael Sollami`
- Issue / ticket: `pending`
- Pull request: `#67 Add UnitBox: the rectangular pulse (box) function`
- Started: `2026-08-18`
- Last updated: `2026-08-18`
- Blocked by: `none`

## Feature Definition

- One-line goal: Add `UnitBox[x]` as a new Mathilda built-in: `1` for `-1/2 <= x <= 1/2`, `0` otherwise.
- Problem: Mathilda lacks the rectangular-pulse function; `Ramp`/`UnitStep` are the closest existing siblings and establish the pattern to follow.
- Requested by: direct user request, on branch `feature/unitbox` (off `main` at `dfaa8aa1df0518dbf7613151196a476411bd2baa`).
- Related links:
  - `thoughts/shared/research/2026-08-18-unitbox-builtin.md` (research doc)

## Working Description

`UnitBox` is a unary, two-sided threshold (box/rectangular-pulse) function, unlike its siblings `UnitStep`/`Ramp` which are one-sided. It will live in `src/piecewise.c` next to `UnitStep`/`Ramp`, reuse `ustep_class()` applied twice (once on `x + 1/2`, once on `1/2 - x`) to answer the two-sided test, and follow the same registration/attribute/docstring/spec/changelog/test conventions as `Ramp`.

## Current State Study

- Relevant existing files:
  - `src/piecewise.c` — `Ramp`/`UnitStep` declaration, registration (`piecewise_init()` lines ~22-45), `ustep_class` machinery (lines ~402-511), `builtin_unitstep`/`builtin_ramp` (~523-593), `piecewise_interval()` helper (~596-601).
  - `src/info.c` — docstrings for `UnitStep`/`Ramp` at lines ~2791-2809.
  - `src/pack.c` — packed-array `AWARE` list, lines ~1000-1029; `UnitStep` present, `Ramp` absent.
  - `src/sym_names.c` — `SYM_UnitStep` interned at line 1429; no `SYM_Ramp`.
  - `src/interval.c` / `.h` — `interval_apply_function()` (`interval.h:73`, `interval.c:906-995`), a `strcmp`-based dispatch (not symbol-pointer based) covering ~30 monotone/classified functions. `Floor`/`Ceiling` call it via `piecewise_interval()`; `UnitStep`/`Ramp` do not, and neither appears in the dispatch table.
  - `docs/spec/builtins/elementary-functions.md` — `UnitStep` entry lines 575-621, `Ramp` entry lines 623-660; insertion point for `UnitBox` is immediately before line 575.
  - `docs/spec/changelog/2026-08-17.md` — current week's changelog file.
  - `tests/test_clip.c` + `tests/test_utils.h` — structural exemplar for a new `tests/test_unitbox.c`.
  - `tests/CMakeLists.txt` — `clip_tests` target (lines ~1276-1277) lacks `add_test()`; header comment (lines 4-7) warns this makes a target invisible to `ctest`.
- Adjacent modules touched: none beyond the above (no `.m` bootstrap changes anticipated).
- Existing behavior to preserve: `UnitStep[0] = 1` (closed boundary precedent, spec line 578); `Ramp`'s precision/numeric-head preservation pattern.
- Constraints from the current codebase: C99-strict, POSIX-guard rules, builtin ownership contract (§4 of SPEC.md).
- Existing tests / validation paths: `tests/test_clip.c` pattern; no `add_test()` precedent to copy — must add it explicitly per user ruling.

## Implementation Spec

- New files to create:
  - `tests/test_unitbox.c`
- Existing files to update:
  - `src/piecewise.c` — `builtin_unitbox()` (reuses `ustep_class()` twice, on `x + 1/2` and `1/2 - x`), registration + attributes in `piecewise_init()`.
  - `src/piecewise.h` — declaration.
  - `src/info.c` — `symtab_set_docstring("UnitBox", ...)`.
  - `docs/spec/builtins/elementary-functions.md` — new `## UnitBox` section before `## UnitStep`.
  - `docs/spec/changelog/2026-08-17.md` — new entry.
  - `tests/CMakeLists.txt` — `add_executable(unitbox_tests ...)` + `add_test(NAME unitbox_tests COMMAND unitbox_tests)`.
- Data models / contracts: none new — reuses the existing `Expr*` builtin-ownership contract and `ustep_class()`'s `USTEP_NEG/NONNEG/UNKNOWN` enum.
- Import directions: `piecewise.c` already includes `sym_names.h`/`eval.h`; no new includes needed.
- Execution flow: `builtin_unitbox` builds `x + 1/2`, evaluates, classifies; if not NEG builds `1/2 - x`, evaluates, classifies; combines per the plan's decision table (NEG on either -> 0; both NONNEG -> 1; else NULL/unevaluated).

## Scope

### In Scope

- `UnitBox[x]` builtin in `src/piecewise.c`/`.h`, registered with `Listable | NumericFunction | Protected`.
- Docstring (`src/info.c`), spec doc entry, changelog entry.
- `tests/test_unitbox.c` + CMake registration with `add_test()`.

### Out Of Scope

- `pack.c` `AWARE` list addition — deferred as unmeasured (explicit user ruling; decision 2).
- `SYM_UnitBox` interned symbol — not added; confirmed not needed since `UnitBox` does not thread through `piecewise_interval()`/`interval_apply_function()` (decision 1).
- `Compile[]`, `D`, `Series`, `Limit`, `Integrate` special-casing of `UnitBox` (mirrors `UnitStep`'s extra call sites) — out of scope for first pass.
- Optimizing the per-element two-`evaluate()`-call cost inside `builtin_unitbox` — accepted, unmeasured, deferred (decision 5).

## Success Criteria

- `AC-1` — `UnitBox[x]` returns the exact integer `1` for `-1/2 <= x <= 1/2` (closed boundary) and `0` outside, for exact/Real/MPFR/certifiable-symbolic-real arguments; left unevaluated for non-real or uncertifiable arguments.
- `AC-2` — `UnitBox` has attributes `Listable`, `NumericFunction`, `Protected`, and explicitly NOT `Orderless`.
- `AC-3` — `UnitBox` has a docstring via `symtab_set_docstring`.
- `AC-4` — `docs/spec/builtins/elementary-functions.md` has a `## UnitBox` section, positioned before `## UnitStep`.
- `AC-5` — the current week's changelog (`docs/spec/changelog/2026-08-17.md`) documents the addition and both explicitly-deferred/accepted tradeoffs (no `pack.c` AWARE entry; accepted per-element evaluate() cost).
- `AC-6` — `tests/test_unitbox.c` exists, is registered in `tests/CMakeLists.txt` with both `add_executable` AND `add_test(NAME ... COMMAND ...)`, and `ctest -R unitbox_tests` finds and runs it.

### Non-functional

- `NFR-1` — no new build warnings under the project's `-Wall -Wextra -Werror=...` flags; `make check-c99` continues to pass (no new POSIX symbols introduced).

## Tests

- Must-have tests: interior points; closed boundary at ±1/2 (rational and Real forms); outside the box; ±Infinity; an exact symbolic real (`Pi`) via certification; symbolic pass-through; complex rejection + zero-imaginary-part resolution; `Listable` threading over a `List`; arity errors (`UnitBox[]`, `UnitBox[1,2]`); attributes (including the *absence* of `Orderless`); a memory-hygiene loop.
- Areas requiring full coverage: `builtin_unitbox`'s three-way branch (NEG-on-lower, NEG-on-upper, both-NONNEG, else-NULL).
- Out of scope for testing: `Interval` threading (feature doesn't have it), `pack.c` buffer fast path (feature doesn't have it).
- Test framework / runner: project's CMake/ctest suite, `tests/test_utils.h` macros (`assert_eval_eq`, `ASSERT`, `ASSERT_MSG`, `TEST`).

## Task List

- [x] Draft implementation plan (files, scope, success criteria) | independent | done
- [x] Phase 1: `builtin_unitbox` in `src/piecewise.c`/`.h`, registration + attributes | independent | done
- [x] Phase 2: docstring (`src/info.c`), spec doc entry, changelog entry | depends on: Phase 1 | done
- [x] Phase 3: `tests/test_unitbox.c` + CMake registration (`add_executable` + `add_test`) | depends on: Phase 1 | done

## Test Results

- Command: `cd tests/build && cmake .. && make unitbox_tests && ./unitbox_tests`
- Outcome: passed
- Summary: 14 test functions run (`Running test: ...` x14), `All UnitBox tests passed.`, exit 0.
- Command: `ctest -R unitbox_tests`
- Outcome: passed — `1/1 Test #61: unitbox_tests ... Passed 0.02 sec`, confirming the `add_test()` wiring works (unlike `clip_tests`).
- Failures fixed: none — clean on first run.
- Known exceptions:
  - `valgrind --leak-check=full` — not run. No `valgrind` on this machine (Apple Silicon macOS). macOS `leaks --atExit` reported "Process ... is not debuggable" (sandboxing/entitlement restriction, not a code finding); its leak count under that condition is not trustworthy and was not treated as a result. Manually traced `builtin_unitbox`'s allocation/free/consume paths (`half`, `lower`, `upper`, `neg_x`) instead — each is freed or consumed exactly once on every branch, matching the codebase's builtin-ownership contract. Flagged for a real valgrind run on Linux CI or a machine that has it installed.

## Checkpoints

- [x] start | completed: `2026-08-18 14:25`
- [x] spec / plan created | completed: `2026-08-18 14:51` (backfilled — plan was approved earlier in-session at `/Users/67840/.claude/plans/zany-pondering-pelican.md`; this stamp records that approval)
- [ ] threat-model stamped | completed: `pending` (not applicable — pure numeric builtin, no security-relevant surface; see Risk Register)
- [x] implementation started | completed: `2026-08-18 14:39` (backfilled — Phase 1 build timestamp)
- [x] implementation complete | completed: `2026-08-18 14:51`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [x] feature validated | completed: `2026-08-18 15:10` (`/ais:validate_plan`: implementation matches the plan, no deviations of substance)
- [x] PR created | completed: `2026-08-18 15:39`
- [ ] closeout complete | completed: `pending`

## PR Updates

- `2026-08-18 15:39` PR opened: https://github.com/stblake/mathilda/pull/67 (branch `feature/unitbox` -> `main`, commit `8b066dca`).

## Decisions

- No `SYM_UnitBox` interned symbol, conditional on the interval-threading decision below. Investigation of `src/interval.c:906-995` shows `interval_apply_function()` dispatches by `strcmp` on the head name, not by symbol pointer — so even routing through it would not itself require an interned symbol. Separately, neither `UnitStep` nor `Ramp` (UnitBox's direct siblings) currently thread through `piecewise_interval()`/`interval_apply_function()` — only `Floor`/`Ceiling` do, specifically because they are monotone. `UnitBox` is a non-monotone two-sided box function and does not fit that mechanism's monotone/region model without new code. Decision: `UnitBox` will NOT thread through `piecewise_interval()` in this first pass, matching its siblings — therefore no `SYM_UnitBox` is needed for this or any other reason currently in scope.
- No `pack.c` `AWARE` list addition — every existing entry is justified by a measured workload per the file's own comment convention; `UnitBox` has none yet. Shipped without it, explicitly deferred rather than silently omitted.
- Closed-interval boundary confirmed: `UnitBox[-1/2] = UnitBox[1/2] = 1`, consistent with `UnitStep[0] = 1` (spec doc line 578).
- New `tests/test_unitbox.c` CMake target will include `add_test()` (not just `add_executable`), correcting the `clip_tests` precedent per the CMakeLists.txt header's own stated intent.

## Risks And Unknowns

- None outstanding — the four open questions from the research doc were resolved by explicit user ruling (see Decisions).

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | N/A — pure numeric builtin, no security-relevant surface | low | low | n/a | n/a | Michael Sollami |

## Dependencies / Blockers

- Dependencies: none
- Blockers: none

## Proofs Of Completion

*(pending)*

## Tech Debt Review

*(pending)*

## Activity Log

- `2026-08-18 14:15` Research completed via `/ais:research_codebase`; doc written to `thoughts/shared/research/2026-08-18-unitbox-builtin.md`.
- `2026-08-18 14:25` Tracking file created; `start` checkpoint stamped. User rulings on the four open questions captured in Decisions. Next: draft and confirm implementation plan.
- `2026-08-18 14:30` Decisions section appended to the research doc (user request), mirroring the four rulings with self-contained reasoning.
- `2026-08-18 14:39` Plan drafted via `EnterPlanMode`/`/ais:create_plan`, refined with a fifth decision (accepted per-element `evaluate()` cost), approved by user. Phase 1 implemented: `builtin_unitbox` in `src/piecewise.c`, declaration in `src/piecewise.h`, registration + attributes in `piecewise_init()`. Build clean, `make check-c99` clean. Smoke-tested via `./Mathilda -file` against all planned manual-verification cases -- all matched.
- `2026-08-18 14:40` `/ais:verify-implementation` run: build/runs PASS, `piecewise_tests` (pre-existing) PASS with no regression, no debug residue; working tree correctly flagged as uncommitted/in-progress (Phases 2-3 pending at that point).
- `2026-08-18 14:51` Phases 2 and 3 completed: docstring added to `src/info.c`; `## UnitBox` spec section added to `docs/spec/builtins/elementary-functions.md` before `## UnitStep` (doc's own List example verified against the binary before writing); changelog entry added to `docs/spec/changelog/2026-08-17.md` covering both deferred/accepted tradeoffs; `tests/test_unitbox.c` written (14 test functions) and registered in `tests/CMakeLists.txt` with both `add_executable` and `add_test()`. `./unitbox_tests` and `ctest -R unitbox_tests` both pass. `valgrind` unavailable on this machine; documented the fallback (manual ownership trace) in Test Results. Tracking file backfilled: Implementation Spec/Scope/Success Criteria/Tests/Task List/Test Results populated; `goal_lock` activated; checkpoints `spec / plan created`, `implementation started`, `implementation complete` stamped.
- `2026-08-18 15:10` `/ais:verify-implementation` re-run post-Phase-3: build/tests/debug-residue all PASS; working tree correctly flagged as the only remaining "not done" item (nothing committed yet). `/ais:validate_plan` run: implementation matches the plan phase-by-phase with no deviations of substance; `feature validated` checkpoint stamped.
- `2026-08-18 15:39` Committed (`8b066dca`, identity set to Michael Sollami / michaelsollami@gmail.com per project convention) and pushed `feature/unitbox` to `origin` (`stblake/mathilda`). Opened PR #67: https://github.com/stblake/mathilda/pull/67. `PR created` checkpoint stamped.

## Reflection

*(pending — filled at closeout via `dev-reflect`)*

## Follow-Up

- `none yet`

## Team Addendum

- `none`
